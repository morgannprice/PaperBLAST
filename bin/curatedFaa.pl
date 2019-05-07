#!/usr/bin/perl -w
use strict;

use DBI;
use FindBin qw{$Bin};
use lib "$Bin/../lib";
use pbutils; # for ReadFastaEntry()
use Getopt::Long;

my $usage = <<END
Usage: curatedFaa.pl -db litsearch.db -uniq uniq.faa -out curated.faa
Filters out curated entries that are probably not actually characterized.
  -filtered filteredFile -- save the non-curated entries to filteredFile
  -showdesc -- include definitions of the curated entries
  -curatedids -- report the curated ids rather than the unique id,
     and save all lengths and descriptions in a out.info file
END
;

my ($dbfile, $uniqfile, $outfile, $filterfile, $showdesc, $curatedids);
die $usage
  unless GetOptions('db=s' => \$dbfile,
                    'uniq=s' => \$uniqfile,
                    'out=s' => \$outfile,
                    'filter=s' => \$filterfile,
                    'showdesc' => \$showdesc,
                   'curatedids' => \$curatedids)
  && defined $dbfile && defined $uniqfile && defined $outfile
  && @ARGV ==0;
die "No such file: $dbfile\n" unless -e $dbfile;
die "No such file: $uniqfile\n" unless -e $uniqfile;
my $dbh = DBI->connect("dbi:SQLite:dbname=$dbfile","","",{ RaiseError => 1 }) || die $DBI::errstr;

my %curatedIds = (); # db::protId => desc
my $genes = $dbh->selectall_arrayref("SELECT db,protId,desc FROM CuratedGene");
my %filtered = (); # db => protId => 1

foreach my $gene (@$genes) {
  my ($db, $protId, $desc) = @$gene;
  # Curated annotations for proteins that are not actually characterized often match these patterns
  # [DUPF]+ matches PF or UPF or DUF, which together with numbers are common uninformative annotations
  # descriptions like protein yyyY occur in EcoCyc
  # These uniformative patterns usually occur at the beginning of the description, but in EcoCyc
  # they may appear after another comment like "Rac prophage;"
  # In CharProtDB, "probable" and "putative" are often not informative
  # In EcoCyc, "putative" is often not informative
  # In other databases, these terms often indicate that there is genetic evidence but
  # not biochemical proof of the exact role; decided to keep those.
  my $maybe = $desc =~ m/^putative/i || $desc =~ m/^probable/i;
  if ($desc =~ m/^uncharacterized/i
      || $desc =~ m/\buncharacterized protein/i
      || $desc =~ m/^(DUF|UPF|PF)\d+ family protein/i
      || $desc =~ m/^(DUF|UPF|PF)\d+ protein/i
      || $desc =~ m/; (DUF|UPF|PF)\d+ family protein/i
      || $desc =~ m/; (DUF|UPF|PF)\d+ protein/i
      || ($desc =~ m/^protein [a-zA-Z]+$/i && $db ne "SwissProt")
      || ($desc =~ m/^putative/i && ($db eq "CharProtDB" || $db eq "ecocyc"))
      || ($desc =~ m/^probable/i && $db eq "CharProtDB")) {
    $filtered{$db}{$protId} = $desc;
  } else {
    $curatedIds{$db . "::" . $protId} = $desc;
  }
}

my %isDup = ();
my $seqdups = $dbh->selectall_arrayref(qq{SELECT sequence_id, duplicate_id FROM SeqToDuplicate WHERE duplicate_id LIKE "%::%"});
my %keepIds = (); # id => list of curated ids that are duplicates of it
foreach my $row (@$seqdups) {
  my ($seqId, $dupId) = @$row;
  if (exists $curatedIds{$dupId}) {
    push @{ $keepIds{$seqId} }, $dupId;
    $isDup{$dupId} = 1;
  }
}

my %written = ();
open(my $fhU, "<", $uniqfile) || die "Cannot read $uniqfile";
open(my $fhOut, ">", $outfile) || die "Cannot write to $outfile";
my $fhInfo;
if (defined $curatedids) {
  open($fhInfo, ">", "$outfile.info") || die "Cannot write to $outfile.info";
  print $fhInfo join("\t", "ids", "length", "descs")."\n";
}

my $state = {};
while (my ($uniqId,$sequence) = ReadFastaEntry($fhU, $state)) {
  if (exists $curatedIds{$uniqId} || exists $keepIds{$uniqId}) {
    if (defined $curatedids) {
      my $ci = FetchCuratedInfo($dbh, $uniqId);
      die "$uniqId is not curated" unless @$ci > 0;
      my @ids = map { $_->[0] . "::" . $_->[1] } @$ci;
      my @desc = map $_->[2], @$ci;
      print $fhOut ">" . join(",", @ids) . "\n" . $sequence . "\n";
      print $fhInfo join("\t", join(",",@ids), length($sequence), join(";; ", @desc))."\n";
    } else {
      my $showid = $uniqId;
      my $desc;
      my $defline = ">" . $showid;
      if (defined $showdesc) {
        my @ids = ();
        push @ids, $uniqId if exists $curatedIds{$uniqId};
        push @ids, @{ $keepIds{$uniqId} } if exists $keepIds{$uniqId};
        foreach my $id (@ids) {
          die unless exists $curatedIds{$id};
          $defline .= " ";
          $defline .= ";; " unless $id eq $ids[0];
          $defline .= "$id $curatedIds{$id}";
        }
      }
      print $fhOut "$defline\n$sequence\n";
    }
    $written{$uniqId} = 1;
  }
}
close($fhU) || die "Error reading $uniqfile";
close($fhOut) || die "Error writing to $outfile";

foreach my $curatedId (keys %curatedIds) {
  die "Did not write a sequence for $curatedId\n"
    unless exists $written{$curatedId} || exists $isDup{$curatedId};
}
print STDERR join(" ", "Found", scalar(keys %curatedIds), "characterized entries and wrote",
                  scalar(keys %written), "different", "sequences", "to", $outfile)."\n";

if (defined $fhInfo) {
  close($fhInfo) || die "Error writing to $outfile.info";
  print STDERR "Wrote $outfile.info\n"
}

if (defined $filterfile) {
  open (my $fh, ">", $filterfile) || die "Cannot write to $filterfile";
  foreach my $db (sort keys %filtered) {
    foreach my $id (sort keys %{ $filtered{$db} }) {
      my $desc = $filtered{$db}{$id};
      print $fh "$db\t$id\t$desc\n";
    }
  }
  close($fh) || die "Error writing to $filterfile";
}
foreach my $db (sort keys %filtered) {
  print STDERR "Filtered out from $db: " . scalar(keys %{ $filtered{$db} }) . " entries that are probably uncharacterized.\n";
}

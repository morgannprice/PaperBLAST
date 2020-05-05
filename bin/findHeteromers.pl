#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw{$Bin};
use lib "$Bin/../lib";
use pbutils qw{ReadTable};
sub ReportProtId($);

my $staticDir = "$Bin/../static";

my $usage = <<END
findHeteromers.pl -sprot sprot.subunits > work/hetero.tab

Compiles heteromer information from BRENDA.subunits,
metacyc.reaction_links, sprot.subunits, TCDB.curated_parsed, and
reanno.curated_parsed

BRENDA.subunits should have fields EC, organism, and one or more
protein identifiers (of the form db::protId)

metacyc.reaction_links should have fields reactionId, reaction name,
and one or more protein identifiers

sprot.subunits should include the fields db, id, isHetero, and
subunit.

Outputs a tab-delimited file with the fields db, protId, and comment
(which has the first sentence of the SUBUNIT field from SwissProt).

By default looks for all input files, except sprot.subunits, in:
  $staticDir
Use -static static to change this.
END
;

my $sprotSubunitsFile;
die $usage
  unless GetOptions('sprot=s' => \$sprotSubunitsFile,
                    'static=s' => \$staticDir)
  && defined $sprotSubunitsFile;
die "Not a directory: $staticDir\n" unless -d $staticDir;
my $brendaFile = "$staticDir/BRENDA.subunits";
my $metacycFile = "$staticDir/metacyc.reaction_links";
my $reannoFile = "$staticDir/reanno.curated_parsed";
my $tcdbFile = "$staticDir/TCDB.curated_parsed";
foreach my $file ($brendaFile, $metacycFile, $reannoFile, $tcdbFile) {
  die "No such file: $file\n" unless -e $file;
}

# To avoid duplicate output rows
my %reported = (); # db => protId => 1

print join("\t", qw{db protId comment})."\n";
open(my $fhB, "<", $brendaFile) || die "Cannot read $brendaFile";
while(my $line = <$fhB>) {
  chomp $line;
  my @F = split /\t/, $line;
  my $EC = shift @F;
  my $organism = shift @F;
  my @protIds = @F;
  if (@protIds > 1) {
    foreach my $protId (@protIds) {
      ReportProtId($protId);
    }
  }
}
close($fhB) || die "Error reading $brendaFile";

open (my $fhM, "<", $metacycFile) || die "Cannot read $metacycFile";
while (my $line = <$fhM>) {
  chomp $line;
  my @F = split /\t/, $line;
  my $reactionId = shift @F;
  my $reactionName = shift @F;
  my @protIds = @F;
  if (@protIds > 1) {
    foreach my $protId (@protIds) {
      ReportProtId($protId);
    }
  }
}
close($fhM) || die "Error reading $metacycFile";

# When this was parsing subunit itself, the logic was:
# next unless $subunit =~ m/^SUBUNIT:/; # a monomer
# $subunit =~ s/^SUBUNIT: *//;
# $subunit =~ s/[.] .*//; # keep the first sentence
# my $startHetero = $subunit =~ m/^hetero[a-z-]*mer/i
#    || $subunit =~ m/^Forms a hetero[a-z-]*mer/i;
# if ($startHetero
#      || !($subunit eq ""
#           || $subunit =~ m/monomer/i
#           || $subunit =~ m/homo[a-z]+mer/i
#           || $subunit =~ m/^interacts/i
#           || $subunit =~ m/^dimer$/i
#           || $subunit =~ m/^dimer of (identical|homo|dimer|trimer|tetramer)/i
#           || $subunit =~ m/^binds/i)) {

my @sprot = ReadTable($sprotSubunitsFile, ["db","id","subunit","isHetero"]);
foreach my $sprot (@sprot) {
  print join("\t", $sprot->{db}, $sprot->{id}, $sprot->{subunit})."\n"
    if $sprot->{isHetero};
}

open (my $fhTCDB, "<", $tcdbFile) || die "Cannot read $tcdbFile\n";
my %tcdbToAcc = (); # tcdbId => list of acc
my %accToTCDB = (); # acc => list of tcdbId
while(my $line = <$fhTCDB>) {
  chomp $line;
  my ($db, $acc, $ids, undef, $desc, $org, $seq, $comments, $pmids) = split /\t/, $line;
  die unless $db eq "TCDB";
  die "Invalid pmids in $line" unless $pmids =~ m/^[0-9,]*$/; # may be missing
  my @tcdbIds = split /,/, $ids;
  foreach my $tcdbId (@tcdbIds) {
    die "Invalid tcdbid $tcdbId" unless $tcdbId =~ m/^[0-9][.][A-Z][.][0-9.]+$/;
    push @{ $tcdbToAcc{$tcdbId} }, $acc;
    push @{ $accToTCDB{$acc} }, $tcdbId;
  }
}
close($fhTCDB) || die "Error reading $tcdbFile";
foreach my $acc (sort keys %accToTCDB) {
  my $isHetero = 0;
  foreach my $tcdbId (@{ $accToTCDB{$acc} }) {
    $isHetero = 1 if @{ $tcdbToAcc{$tcdbId} } > 1;
  }
  ReportProtId("TCDB::" . $acc) if $isHetero;
}

open (my $fhReanno, "<", $reannoFile) || die "Cannot read $reannoFile\n";
while (my $line = <$fhReanno>) {
  chomp $line;
  my ($db, $acc, undef, undef, $desc, $org) = split /\t/, $line;
  die unless $db eq "reanno";
  ReportProtId("reanno::" . $acc)
    if $desc =~ m/\b(component|subunit)s?\b/i;
}
close($fhReanno) || die "Error reading $reannoFile";

sub ReportProtId($) {
  my $protId = shift;
  my @parts = split /::/, $protId;
  die "Invalid protId $protId" unless @parts == 2;
  my ($db,$id) = @parts;
  print join("\t", @parts, "")."\n"
    unless exists $reported{$db}{$id};
  $reported{$db}{$id} = 1;
}

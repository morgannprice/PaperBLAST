#!/usr/bin/perl -w
# Build the sqlite3 database of curated proteins for GapMind and Curated Clusters
use strict;
use DBI;
use Getopt::Long;
use FindBin qw{$RealBin};
use lib "$RealBin/../lib";
use pbutils qw{ReadFastaEntry ReadTable SqliteImport};
sub DescToEc($);

my $staticDir = "$RealBin/../static";
my $usage = <<END
buildCuratedDb.pl -dir tmp/path.aa

Assumes that the directory already contains curated.faa,
curated.faa.info, curated2.faa, and hetero.tab. If you don't want to
inlude pfam hits (used by curatedClusters.cgi, but not by GapMind
itself), set pfam.hits.tab to be empty.

Also uses some files from static/ -- metacyc.reaction_compounds,
metacyc.reaction_links, TCDB.curated_parsed

Optional arguments:
-curated dir/curated.faa
-curated2 dir/curated2.faa
-hetero dir/hetero.tab
-pfamhits dir/pfam.hits.tab
-static $staticDir
-out dir/curated.db
END
;

my ($dir, $curatedFile, $curated2File, $heteroFile, $pfamHitsFile, $dbFile);
die $usage
  unless GetOptions('dir=s' => \$dir,
                    'curated=s' => \$curatedFile,
                    'curated2=s' => \$curated2File,
                    'hetero=s' => \$heteroFile,
                    'pfamhits=s' => \$pfamHitsFile,
                    'static=s' => \$staticDir,
                    'out=s' => \$dbFile)
  && @ARGV == 0
  && defined $dir;

die "No such directory: $dir\n" unless -d $dir;
die "No such directory: $staticDir\n" unless -d $staticDir;

$curatedFile = "$dir/curated.faa" unless defined $curatedFile;
my $curatedInfoFile = "$curatedFile.info";
$curated2File = "$dir/curated2.faa" unless defined $curated2File;
$heteroFile = "$dir/hetero.tab" unless defined $heteroFile;
$pfamHitsFile = "$dir/pfam.hits.tab" unless defined $pfamHitsFile;
$dbFile = "$dir/curated.db" unless defined $dbFile;
my $reactionLinksFile = "$staticDir/metacyc.reaction_links";
my $reactionCompoundsFile = "$staticDir/metacyc.reaction_compounds";
my $tcdbFile = "$staticDir/TCDB.curated_parsed";

foreach my $file ($curatedFile, $curatedInfoFile, $curated2File, $heteroFile, $pfamHitsFile,
                  $reactionLinksFile, $reactionCompoundsFile, $tcdbFile) {
  die "No such file: $file\n" unless -e $file;
}

my $tmpDir = $ENV{TMP} || "/tmp";
my $tmpDbFile = "$tmpDir/buildCuratedDb.$$.db";
print STDERR "Building temporary database $tmpDbFile\n";

my $schema = "$RealBin/../lib/curated.sql";
system("sqlite3 $tmpDbFile < $schema") == 0
  || die "Error loading schema $schema into $tmpDbFile -- $!";

# the orgs and id2s fields are optional
my @info = ReadTable($curatedInfoFile, ["ids", "length", "descs"]);
# some other inputs have just 1 id, instead of the full ids, so build
# the mapping
my %idToIds = ();
my @idToIds = ();
foreach my $info (@info) {
  my $curatedIds = $info->{ids};
  foreach my $id (split /,/, $curatedIds) {
    $idToIds{$id} = $curatedIds;
    push @idToIds, [ $id, $curatedIds ];
  }
}
my @curatedInfo = map { $_->{descs} =~ s/\r */ /g;
                        [ $_->{ids}, $_->{length}, $_->{descs},
                          $_->{id2s} || "", $_->{orgs} || "" ]
                      } @info;
SqliteImport($tmpDbFile, "CuratedInfo", \@curatedInfo);
print STDERR "Loaded CuratedInfo\n";
SqliteImport($tmpDbFile, "CuratedIdToIds", \@idToIds);

my @curatedSeq = ();
open(my $fhFaa, "<", $curatedFile) || die "Cannot read $curatedFile\n";
my $state = {};
while (my ($curatedIds, $seq) = ReadFastaEntry($fhFaa, $state)) {
  push @curatedSeq, [ $curatedIds, $seq ];
}
close($fhFaa) || die "Error reading $curatedFile";
SqliteImport($tmpDbFile, "CuratedSeq", \@curatedSeq);
print STDERR "Loaded CuratedSeq\n";

my @curated2 = ();
open($fhFaa, "<", "$curated2File") || die "Cannot read $curated2File\n";
$state = {};
while (my ($header, $seq) = ReadFastaEntry($fhFaa, $state)) {
  my @F = split / /, $header;
  my $protId = shift @F;
  my $desc = join(" ", @F);
  push @curated2, [ $protId, $desc, $seq ];
}
close($fhFaa) || die "Error reading $curated2File";
SqliteImport($tmpDbFile, "Curated2", \@curated2);

# heteroFile has one entry per id, so the same sequence could be listed
# several times (from different databases)
my @heteroIn = ReadTable($heteroFile, ["db","protId","comment"]);
my %hetero = (); # curatedIds to list of comments
my $nUnknownHetero = 0;
foreach my $row (@heteroIn) {
  my $id = $row->{db} . "::" . $row->{protId};
  if (exists $idToIds{$id}) {
    push @{ $hetero{$idToIds{$id}} }, $row->{comment};
  } else {
    $nUnknownHetero++;
  }
}
print STDERR "Warning: skipped $nUnknownHetero entries from $heteroFile with unknown protein ids\n"
  if $nUnknownHetero > 0;
my @hetero = (); # curatedIds to comment
foreach my $curatedIds (sort keys %hetero) {
  # Most comments are empty (they just show that protein is part of a complex)
  # As of February 2021, only SwissProt has comments, so there wouldn't be
  # more than one, but in case it happens, join them together
  my @comments = grep { $_ ne "" } @{ $hetero{$curatedIds} };
  push @hetero, [ $curatedIds, join(". ", @comments) ];
}
SqliteImport($tmpDbFile, "Hetero", \@hetero);
print STDERR "Loaded Hetero\n";

my @pfamHits = ();
open (my $fhHits, "<", $pfamHitsFile) || die "Cannot read $pfamHitsFile\n";
while (my $line = <$fhHits>) {
  chomp $line;
  my ($protId, $hmmName, $hmmAcc, $eval, $bits,
      $protBeg, $protEnd, $protLen,
      $hmmBeg, $hmmEnd, $hmmLen) = split /\t/, $line;
  die unless defined $hmmLen && $hmmLen =~ m/^\d+$/;
  push @pfamHits, [ $protId, $hmmName, $hmmAcc, $eval, $bits,
                    $protBeg, $protEnd, $protLen,
                    $hmmBeg, $hmmEnd, $hmmLen ];
}
SqliteImport($tmpDbFile, "CuratedPFam", \@pfamHits);
print STDERR "Loaded CuratedPFam\n";

my @compoundInReaction = ();
my %crKey = (); # rxnId::cmpId::side must be unique
open (my $fhCompounds, "<", $reactionCompoundsFile)
  || die "Cannot read $reactionCompoundsFile";
while (my $line = <$fhCompounds>) {
  chomp $line;
  my @F = split /\t/, $line;
  my $rxnId = shift @F;
  my $rxnLocation = shift @F;
  foreach my $spec (@F) {
    my ($side, $coeff, $compartment, $cmpId, $cmpDesc) = split /:/, $spec;
    die unless defined $cmpDesc;
    my $key = join("::", $rxnId, $cmpId, $side);
    push @compoundInReaction, [ $rxnId, $rxnLocation, $cmpId, $cmpDesc,
                                $side, $coeff, $compartment ]
      unless exists $crKey{$key};
    $crKey{$key} = 1;
  }
}
close($fhCompounds) || die "Error reading $reactionCompoundsFile";
SqliteImport($tmpDbFile, "CompoundInReaction", \@compoundInReaction);

my @enzymeForReaction = ();
my %erKey = (); # curatedIds ::: $rxnId should be unique
open(my $fhRxnLinks, "<", $reactionLinksFile)
  || die "Error reading $reactionLinksFile";
my $nEnzSkip = 0;
while (my $line = <$fhRxnLinks>) {
  chomp $line;
  my @F = split /\t/, $line;
  my $rxnId = shift @F;
  my $enzDesc = shift @F;
  foreach my $id (@F) {
    if (!exists $idToIds{$id}) {
      $nEnzSkip++;
    } else {
      my $curatedIds = $idToIds{$id};
      my $key = join(":::", $curatedIds, $rxnId);
      push @enzymeForReaction, [ $curatedIds, $rxnId, $enzDesc ]
        unless exists $erKey{$key};
      $erKey{$key} = 1;
    }
  }
}
close($fhRxnLinks) || die "Error reading $reactionLinksFile";
print STDERR "Warning: skipped $nEnzSkip entries from $reactionLinksFile with unknown protein ids\n"
  if $nEnzSkip > 0;
SqliteImport($tmpDbFile, "EnzymeForReaction", \@enzymeForReaction);

# Transporters to substrates
my @transporterSubstrate = ();
open(my $fhTCDB, "<", $tcdbFile) || die "Cannot read $tcdbFile";
my $nSkipTCDB = 0;
while (my $line = <$fhTCDB>) {
  chomp $line;
  my ($db, $protId, $tcdb_class, undef, $desc, $org, $seq, $comment, $pmIds) = split /\t/, $line;
  die unless defined $pmIds;
  die unless $db eq "TCDB";
  my $id = "$db" . "::" . $protId;
  if (exists $idToIds{$id}) {
    my @comments = split /_:::_/, $comment;
    foreach my $c (@comments) {
      if ($c =~ m/^SUBSTRATES: (.*)$/) {
        my $substrate = $1;
        push @transporterSubstrate, [ $idToIds{$id}, $substrate ];
      }
    }
  } else {
    $nSkipTCDB++;
  }
}
close($fhTCDB) || die "Error reading $tcdbFile";
print STDERR "Skipped $nSkipTCDB unknown identifiers in $tcdbFile\n"
  if $nSkipTCDB > 0;
SqliteImport($tmpDbFile, "TransporterSubstrate", \@transporterSubstrate);

# EC numbers for each curated item
my %ecCurated = (); # ec => curatedIds => 1
foreach my $info (@info) {
  my @descs = split /;; /, $info->{descs};
  foreach my $desc (@descs) {
    foreach my $ec (DescToEc($desc)) {
      $ecCurated{$ec}{ $info->{ids} } = 1;
    }
  }
}
my @ecCurated = ();
foreach my $ec (sort keys %ecCurated) {
  my @curatedIds = sort keys %{ $ecCurated{$ec} };
  foreach my $curatedIds (@curatedIds) {
    push @ecCurated, [ $ec, $curatedIds ];
  }
}
SqliteImport($tmpDbFile, "ECToCurated", \@ecCurated);

my %ecCurated2 = (); # ec => protId => 1
foreach my $row (@curated2) {
  my ($protId, $desc, $seq) = @$row;
  foreach my $ec (DescToEc($desc)) {
    $ecCurated2{$ec}{$protId} = 1;
  }
}
my @ecCurated2 = ();
foreach my $ec (sort keys %ecCurated2) {
  my @protIds = sort keys %{ $ecCurated2{$ec} };
  foreach my $protId (@protIds) {
    push @ecCurated2, [ $ec, $protId ];
  }
}
SqliteImport($tmpDbFile, "ECToCurated2", \@ecCurated2);

system("cp $tmpDbFile $dbFile") == 0 || die "Copying $tmpDbFile to $dbFile failed: $!";
unlink($tmpDbFile);
print STDERR "Built curated database $dbFile\n";

# Does not guarantee that results are unique
# Finds EC numbers of the form EC 1.1.1.1 or EC:1.1.1.1, delimited with
# spaces, commas, semicolons, spaces, or brackets i.e. " EC 1.1.1.1,"
# It also allows identifiers like this 3.2.1.20|3.2.1.28 (appears in CAZy)

sub DescToEc($) {
  my ($desc) = @_;
  my @ecSeen = ();
  my @words = split /[ |]/, $desc;
  foreach my $word (@words) {
    # Removing leading ( or [
    $word =~ s/^[\(\[]//;
    # Remove leading EC: (the space case does not arise)
    $word =~ s/^EC://i;
    # Remove trailing , ; ) ]
    $word =~ s/[,;\)\]]+$//;
    # allow non-specific identifiers like 1.1.1.-
    # allow BRENDA identifiers like 3.2.1.B4, sucrose-6-phosphate hydrolase
    # or metacyc identifiers like 4.4.1.m
    push @ecSeen, $word
      if $word =~ m/^\d+[.][0-9-]+[.][0-9-]+[.][0-9-]+$/
        || $word =~ m/^\d+[.][0-9]+[.][0-9]+[.][A-Z][0-9]+$/
        || $word =~ m/^d+[.]\d+[.]\d+[.][a-z]$/;
  }
  return @ecSeen;
}




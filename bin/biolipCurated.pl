#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw{$RealBin};
use lib "$RealBin/../lib";
use pbutils qw{readPdbNames};

my $usage = <<END
biolipCurated.pl -biolip BioLiP_2013-03-6_nr.txt BioLiP_UP_nr.txt
   -pdbnames pdbNames.tsv
   -pdbligands PDBLigands.tab > biolip.curated_parsed

Produces the biolip curated_parsed file.
END
;

my @biolipFiles;
my ($pdbNamesFile, $pdbLigandsFile);

die $usage
  unless GetOptions('biolip=s{1,}' => \@biolipFiles,
                    'pdbnames=s' => \$pdbNamesFile,
                    'pdbligands=s' => \$pdbLigandsFile)
  && @ARGV == 0
  && @biolipFiles > 0
  && defined $pdbNamesFile
  && defined $pdbLigandsFile;

my %pdbDesc = %{ readPdbNames($pdbNamesFile) };
print STDERR "Read " . scalar(keys %pdbDesc) . " pdb descriptions from $pdbNamesFile\n";

# id of the form 1a1cB => hash with
#    desc, uniprotId, seq, pubmedIds (comma-delimited string) and ligandIds (which is a list)
my %idInfo = ();

foreach my $file (@biolipFiles) {
  open(my $fhAnno, "<", $file) || die "Cannot read $file\n";
  while (my $line = <$fhAnno>) {
    chomp $line;
    my ($pdbId, $chain, $resolution, $siteId, $ligandId, $ligandChain, $ligandNo,
        $bindingResiduesPDB, $bindingResidues,
        $siteResiduesPDB, $siteResidues,
        $EC, $GO,
        $affinityPM, $affinityMOAD, $affinityPDBBind, $affinityBindingDB,
        $uniprotId, $pubmedIds, $ligandResidues, $seq) = split /\t/, $line;
    die "Invalid line $line" unless $seq =~ m/^[A-Za-z]+$/;
    $seq = uc($seq); # not sure why there is a lowercase character in 1aw8B
    $pubmedIds =~ m/^[0-9,]*$/ || die "Invalid pubmed ids $pubmedIds in line\n$line";
    my $id = $pdbId.$chain;
    $idInfo{$id}{uniprotId} = $uniprotId;
    $idInfo{$id}{desc} = $pdbDesc{$pdbId} || "";
    $idInfo{$id}{seq} = $seq;
    $idInfo{$id}{pubmedIds} = $pubmedIds;
    push @{ $idInfo{$id}{ligandIds} }, $ligandId;
  }
}
print STDERR "Read biolip information for " . scalar(keys %idInfo) . " chains\n";

my %ligandDesc = ();
open (my $fhLigand, "<", $pdbLigandsFile) || die "Cannot read $pdbLigandsFile\n";
while(my $line = <$fhLigand>) {
  chomp $line;
  my ($ligandId, $ligandDesc) = split /\t/, $line;
  die $line unless $ligandDesc;
  $ligandDesc{$ligandId} = $ligandDesc;
}
close($fhLigand) || die "Error reading $pdbLigandsFile";
print STDERR "Read descriptions for " . scalar(keys %ligandDesc) . " ligands\n";

# In rare cases, chain identifiers differ only by case, i.e. 1rxs:A and 1rxs:a
# These will cause problems downstreai with BLAST tools.
# So, suppress any ids already seen as lower case
my %lcidSeen = ();
my $nSkipLC = 0;
my $nWritten = 0;
foreach my $id (sort keys %idInfo) {
  if (exists $lcidSeen{lc($id)}) {
    $nSkipLC++;
    next;
  }
  $lcidSeen{lc($id)} = 1;
  my $info = $idInfo{$id};
  my %ligandSeen = ();
  my @ligandIds = grep { my $keep = !exists $ligandSeen{$_};
                         $ligandSeen{$_} = 1;
                         $keep } @{ $info->{ligandIds} };
  my @ligandDescs = map { exists $ligandDesc{$_} ? lc($ligandDesc{$_}) : $_ } @ligandIds;
  my $comment = "";
  $comment = (@ligandDescs > 1 ? "Ligands:" : "Ligand:") . " " . join("; ", @ligandDescs)
    if @ligandDescs > 0;
  print join("\t", "biolip",
             $id, $info->{uniprotId},
             "", # short_name
             $info->{desc},
             "", # organism
             $info->{seq}, $comment, $info->{pubmedIds})."\n";
  $nWritten++;
}
print STDERR "Skipped $nSkipLC entries due to non-uniqueness after matching case\n";
print STDERR "Wrote $nWritten lines in curated_parsed format\n";


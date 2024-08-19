#!/usr/bin/perl -w

# Required parameters:
# set -- which database of curated proteins (defaults to gaps2, i.e. using ../tmp/path.gaps2/curated.faa*)
# ids (multiple values) -- each value is an ids (comma-delimited curated sequence identifier)
#    or is of the form UniProt::id, where this uniprot is described in the query file for
#    this pathway
#
# Optional arguments:
# path -- which pathway is under consideration (required if any ids are from UniProt)

use strict;
use CGI qw(:standard Vars start_ul end_ul);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use Time::HiRes qw{gettimeofday};
use DBI;
use lib "../lib";
use pbutils qw{ReadTable NewerThan};
use pbweb qw{start_page AddCuratedInfo GeneToHtmlLine};
use Steps qw{ReadSteps};
use DB_File;
use List::Util qw{min max};
use IO::Handle qw{autoflush};

sub CompoundInfoToHtml($$$);

my $minCoverageFrac = 0.5;

my $set = param('set') || "gaps2";
$set =~ m/^[a-zA-Z0-9._-]+$/ || die "Invalid set $set";

my $stepPath = "../gaps/$set";
die "Invalid set $set: no $stepPath directory" unless -d $stepPath;
my $queryPath = "../tmp/path.$set"; # intermediate files
die "Invalid set $set: no $queryPath directory" unless -d $queryPath;

my @ids = param('ids');
die "Must specify at least two identifiers" unless @ids > 1;

my $curatedDb = "$queryPath/curated.db";
die "No such file: $curatedDb" unless -e $curatedDb;
my $dbhC = DBI->connect("dbi:SQLite:dbname=$curatedDb","","",{ RaiseError => 1 }) || die $DBI::errstr;
my $stepsDb = "$queryPath/steps.db";
my $dbhS = DBI->connect("dbi:SQLite:dbname=$stepsDb","","",{ RaiseError => 1 }) || die $DBI::errstr;

my $pathSpec = param('path') || "";
die "Invalid pathspec $pathSpec" unless $pathSpec =~ m/^[a-zA-Z_0-9-]*$/;
my ($banner, $bannerURL);

my $blastall = "../bin/blast/blastall";
my $formatdb = "../bin/blast/formatdb";
foreach my $x ($blastall,$formatdb) {
  die "No such executable: $x" unless -x $x;
}

my $tmpPre = "/tmp/curatedSim.$$";

my ($setDesc) = $dbhS->selectrow_array("SELECT desc FROM Pathway WHERE pathwayId = 'all'");

start_page('title' => 'Similarities of Characterized Proteins',
           'banner' => $banner, 'bannerURL' => $bannerURL);
autoflush STDOUT 1; # show preliminary results

my %seqs = ();
my %curatedInfo = ();
foreach my $ids (@ids) {
  if ($ids =~ m/^UniProt::(.*)$/) {
    my $uniprotId = $1;
    die "Cannot handle uniprot inputs unless path is set" if $pathSpec eq "";
    my $row = $dbhS->selectrow_hashref(qq{ SELECT * FROM StepQuery
                                        WHERE pathwayId = ?
                                        AND queryType = "uniprot"
                                        AND uniprotId = ? },
                                       {}, $pathSpec, $uniprotId);
    die "Unknown uniprotId $uniprotId in pathway $pathSpec" unless defined $row;
    $seqs{$ids} = $row->{seq};
    # Fake a row from CuratedInfo
    $curatedInfo{$ids} = { 'curatedIds' => $ids, 'seqLength' => length($row->{seq}),
                           'descs' => $row->{desc}, 'ids2s' => $uniprotId, 'orgs' => "" };
  } else {
    my ($seq) = $dbhC->selectrow_array("SELECT seq FROM CuratedSeq WHERE curatedIds = ?",
                                       {}, $ids);
    die "No sequence for $ids" unless defined $seq;
    $seqs{$ids} = $seq;
    my $row = $dbhC->selectrow_hashref("SELECT * FROM CuratedInfo WHERE curatedIds = ?",
                                       {}, $ids);
    die "Unknown curated ids $ids" unless defined $row;
    $curatedInfo{$ids} = $row;
  }
}

print p("Comparing", scalar(keys %seqs), "sequences");

# from multi-identifier "ids" such as "CharProtDB::CH_091614,SwissProt::P24942,TCDB::P24942"
# to a short id such as "CharProtDB::CH_091614", and back

my %shortToIds = ();
my %idsToShort = ();
foreach my $ids (@ids) {
  my $short = $ids; $short =~ s/,.*//;
  die "Duplicate id prefix $short" if exists $shortToIds{$short};
  $idsToShort{$ids} = $short;
  $shortToIds{$short} = $ids;
}

my @subjectIds = @ids;
my $queryId = shift @subjectIds || die;

my $queryFile = "$tmpPre.query";
open(my $fhQ, ">", $queryFile) || die "Cannot write to $queryFile";
print $fhQ ">", $idsToShort{$queryId}, "\n", $seqs{$queryId}, "\n";
close($fhQ) || die "Error writing $queryFile";

my $subjectFile = "$tmpPre.subject";
open(my $fhS, ">", $subjectFile) || die "Cannot write to $subjectFile";
foreach my $subjectId (@subjectIds) {
  print $fhS ">", $idsToShort{$subjectId}, "\n", $seqs{$subjectId}, "\n";
}
close($fhS) || die "Error writing $subjectFile";

system("$formatdb -p T -i $subjectFile > /dev/null 2>&1") == 0
  || die "$formatdb for $subjectFile failed: $!";
my $hitsFile = "$tmpPre.hits";
system(qq{$blastall -p blastp -d $subjectFile -i $queryFile -e 1e-3 -m 8 -F "m S" -o $hitsFile > /dev/null 2>&1}) == 0
  || die "$blastall on $subjectFile $queryFile failed: $!";

unlink($queryFile);
unlink($subjectFile);
foreach my $suffix (qw{phr pin psq}) {
  unlink("$subjectFile.$suffix");
}

my @hits = ();
my %subjectHit = ();
open(my $fhH, "<", $hitsFile) || die "Error reading $hitsFile";
while (my $line = <$fhH>) {
  chomp $line;
  my ($qShort, $sShort, $identity, $alen, $mm, $gap, $qbeg, $qend, $sbeg, $send, $eval, $bits) = split /\t/, $line;
  die unless $qShort eq $idsToShort{$queryId};
  my $subjectId = $shortToIds{$sShort} || die;
  next unless ($qend-$qbeg+1) >= $minCoverageFrac * length($seqs{$queryId})
    && ($send-$sbeg+1) >= $minCoverageFrac * length($seqs{$subjectId});
  $subjectHit{$subjectId} = { 'qbeg' => $qbeg, 'qend' => $qend, 'sbeg' => $sbeg, 'send' => $send,
                              'identity' => $identity }
    unless exists $subjectHit{$subjectId};
}
close($fhH) || die "Error reading $hitsFile";

unlink($hitsFile);

print h3("Query"),
  p(CompoundInfoToHtml($queryId, $curatedInfo{$queryId}, $seqs{$queryId})),
  "\n";

my @subjectsHit = sort { $subjectHit{$b}{identity} <=> $subjectHit{$a}{identity} }
  grep { exists $subjectHit{$_} } @subjectIds;

print h3("Other Sequences with Hits"), "\n";
print p("None") if @subjectsHit == 0;

foreach my $subjectId (@subjectsHit) {
  my $simString = "";
  if (exists $subjectHit{$subjectId}) {
    my $hit = $subjectHit{$subjectId};
    my $qLen = length($seqs{$queryId});
    my $sLen = length($seqs{$subjectId});
    my $covFracQ = ($hit->{qend} - $hit->{qbeg} + 1) / $qLen;
    my $covFracS = ($hit->{send} - $hit->{sbeg} + 1) / $sLen;

    $simString = b(int($hit->{identity}) . "%") . " identical to query, "
      . a({-title => "$hit->{sbeg}:$hit->{send}/$sLen of this sequence is similar to"
           . " $hit->{qbeg}:$hit->{qend}/$qLen of query"},
          b(int(100 * min($covFracQ,$covFracS)) . "%") . " coverage");
  } else {
    $simString = "No BLAST hit to query";
  }
  print p(CompoundInfoToHtml($subjectId, $curatedInfo{$subjectId}, $seqs{$subjectId}),
          br(), $simString), "\n";
}

my @subjectsNotHit = grep !exists $subjectHit{$_}, @subjectIds;
print h3("Other Sequences without Hits"), "\n";
print p("None") if @subjectsNotHit == 0;

foreach my $subjectId (@subjectsNotHit) {
  print p(CompoundInfoToHtml($subjectId, $curatedInfo{$subjectId}, $seqs{$subjectId})), "\n";
}

print end_html;
exit(0);


sub CompoundInfoToHtml($$$) {
  my ($compoundId, $info, $seq) = @_;
  die unless $compoundId;
  die "No info for $compoundId" unless $info;
  die "no seq for $compoundId" unless $seq;
  my @ids = split /,/, $compoundId;
  die unless @ids > 0;
  my @descs = split /;; /, $info->{descs};
  my @orgs = split /;; /, $info->{orgs} if defined $info->{orgs};
  my @id2s = split /;; /, $info->{id2s} if defined $info->{id2s};
  die "Mismatched length of ids and descs" unless scalar(@ids) == scalar(@descs);
  my $len = $info->{seqLength};
  die unless $len;
  my ($hetComment) = $dbhC->selectrow_array("SELECT comment FROM Hetero WHERE curatedIds = ?",
                                            {}, $compoundId);


  my @genes;
  for (my $i = 0; $i < @ids; $i++) {
    my $id = $ids[$i];
    my $desc = $descs[$i];
    my ($db, $protId) = split /::/, $id;
    die "Cannot parse id $id" unless $protId;
    my $gene = { 'db' => $db, 'protId' => $protId, 'desc' => $desc,
                 'protein_length' => $len,
                 'comment' => '', 'name' => '', id2 => '', organism => '' };
    $gene->{id2} = $id2s[$i] if @id2s > 0;
    $gene->{organism} = $orgs[$i] if @orgs > 0;
    AddCuratedInfo($gene);
    $gene->{HTML} = GeneToHtmlLine($gene);
    $gene->{HTML} .= " (" . a({ -title => $hetComment }, "heteromeric") . ")"
      if defined $hetComment;
    push @genes, $gene;
  }
  @genes = sort { $a->{priority} <=> $b->{priority} } @genes;
  my @pieces = map $_->{HTML}, @genes;
  my $id1 = $genes[0]->{db} . "::" . $genes[0]->{protId};
  my $newline = "%0A";
  my $query = ">$id1$newline$seq";
  my @links = ();
  push @links, a({-href => "litSearch.cgi?query=$query"}, "PaperBLAST");
  push @links, a({-href => "http://www.ncbi.nlm.nih.gov/Structure/cdd/wrpsb.cgi?seqinput=$query"}, "CDD");
  return join("<BR>", @pieces,
                small($len, "amino acids: ", join(", ", @links)));
}

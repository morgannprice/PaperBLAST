#!/usr/bin/perl -w
use strict;

# Optional parameters:
# set -- which database of curated proteins (defaults to gaps2, i.e. using ../tmp/gaps2.aa)
# query -- what term to search for
# word -- report whole word matches only (like a perl boolean)
# identity -- %identity threshold for clustering, default 30
# coverage -- %coverage threshold for clustering, default 75

use strict;
use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use Time::HiRes qw{gettimeofday};
use DBI;
use lib "../lib";
use pbutils;
use pbweb qw{start_page AddCuratedInfo GeneToHtmlLine};
use DB_File;
use URI::Escape;
use HTML::Entities;
use IO::Handle qw{autoflush};

sub MatchRows($$$);
sub FaaToDb($$);
sub CompoundInfoToHtml($$);

my $maxHits = 250;

my $set = param('set') || "gaps2";
my $query = param('query') || "";
my $word = param('word') || 0;

my $minIdentity = param('identity');
$minIdentity = 30 unless defined $minIdentity && $minIdentity =~ m/^\d+$/;
$minIdentity = 100 if $minIdentity > 100;
$minIdentity = 0 if $minIdentity < 0;

my $minCoverage = param('coverage');
$minCoverage = 75 unless defined $minCoverage && $minCoverage =~ m/^\d+$/;
$minCoverage = 100 if $minCoverage > 100;
$minCoverage = 0 if $minCoverage < 0;

$set =~ m/^[a-zA-Z0-9._-]+$/ || die "Invalid set $set";
my $queryPath = "../tmp/path.$set"; # intermediate files
die "Invalid set $set: no $queryPath directory" unless -d $queryPath;
my $curatedFaa = "$queryPath/curated.faa";
die "No such file: $curatedFaa" unless -e $curatedFaa;

my $fastacmd = "../bin/blast/fastacmd";
my $blastall = "../bin/blast/blastall";
my $formatdb = "../bin/blast/formatdb";
foreach my $x ($fastacmd,$blastall,$formatdb) {
  die "No such executable: $x" unless -x $x;
}

my $tmpPre = "/tmp/cluratedClusters.$$";

start_page('title' => 'Clusters of Characterized Proteins');
autoflush STDOUT 1; # show preliminary results

unless (NewerThan("$curatedFaa.db", $curatedFaa)) {
  print p("Reformatting the sequence database"), "\n";
  FaaToDb($curatedFaa, "$curatedFaa.db");
}

$query =~ s/^\s+//;
$query =~ s/\s+$//;
if ($query eq "") {
  print
    start_form(-method => 'get', -action => 'curatedClusters.cgi'),
    hidden(-name => 'set', -value => $set, -override => 1),
    p("Search for:", textfield(-name => 'query', -value => '', -size => 50, -maxlength => 200)),
    p(small("Use % as a wild card character")),
    p({-style => "margin-left: 3em;" },
      checkbox(-name => "word", -label => "Match whole words only?", -checked => 1)),
    p("Cluster characterized proteins at",
      textfield(-name => "identity", -value => $minIdentity, -size => 3, -maxlength => 3),
      '%identity and',
      textfield(-name => "coverage", -value => $minCoverage, -size => 3, -maxlemgth => 3),
      '%coverage'),
    p(submit(-name => 'Search')),
    end_form;
} else {
  print p("Searching for", HTML::Entities::encode($query)),
    p("Or try", a({-href => "curatedClusters.cgi?set=${set}"}, "another search")),
    "\n";
  my @curatedInfo = ReadTable("$curatedFaa.info", ["ids","length","descs"]);
  my %curatedInfo = map { $_->{ids} => $_ } @curatedInfo;
  my $hits = MatchRows(\@curatedInfo, $query, $word);
  if (@$hits > $maxHits) {
    print p("Found over $maxHits proteins with matching descriptions. Please try a different query.");
    print end_html;
    exit(0);
  } elsif (@$hits == 0) {
    print p("Sorry, no hits were found");
  } else {
    print p("Found " . scalar(@$hits) . " characterized proteins with matching descriptions.");
  }
  print "\n";
  if (@$hits > 0) {
    # Fetch their sequences
    my %seqsAll;
    tie %seqsAll, "DB_File", "$curatedFaa.db", O_RDONLY, 0666, $DB_HASH
      || die "Cannot open file $curatedFaa.db -- $!";
    my %seqs = ();
    foreach my $id (map $_->{ids}, @$hits) {
      $seqs{$id} = $seqsAll{$id} || die "No or empty sequence for $id in $curatedFaa.db";
    }
    untie %seqsAll;
    die unless scalar(keys %seqs) == scalar(@$hits);
    print p("Fetched " . scalar(keys %seqs) . " sequences");

    open(my $fhFaa, ">", "$tmpPre.faa") || die "Error writing to $tmpPre.faa";
    foreach my $id (sort keys %seqs) {
      print $fhFaa ">$id\n$seqs{$id}\n";
    }
    close($fhFaa) || die "Error writing to $tmpPre.faa";
    print p("Running BLASTp"), "\n";
    my $covFrac = $minCoverage / 100;
    my $formatCmd = "$formatdb -i $tmpPre.faa -p T";
    system($formatCmd) == 0 || die "formatdb failed -- $formatCmd -- $!";
    my $blastCmd = qq{$blastall -p blastp -i $tmpPre.faa -d $tmpPre.faa -F "m S" -e 0.001 -m 8 -a 8 -o $tmpPre.hits >& /dev/null};
    system($blastCmd) == 0 || die "blastall failed -- $blastCmd -- $!";
    unlink("$tmpPre.faa");
    foreach my $suffix (qw{phr pin psq}) {
      unlink("$tmpPre.faa.$suffix");
    }
    my %sim = (); # id => id => bits if hit (removed self hits)
    my @sims = (); # list of [id1, id2, bits]
    open (my $fhsim, "<", "$tmpPre.hits") || die "Cannot read $tmpPre.hits";
    while(my $line = <$fhsim>) {
      my ($id1, $id2, $identity, $alnlen, $mm, $gap, $qbeg, $qend, $sbeg, $send, $eval, $bits) = split /\t/, $line;
      die $id1 unless exists $seqs{$id1};
      die $id2 unless exists $seqs{$id2};
      my $cov1 = $alnlen / length($seqs{$id1});
      my $cov2 = $alnlen / length($seqs{$id2});
      if ($id1 ne $id2 && $identity >= $minIdentity && $cov1 >= $covFrac && $cov2 >= $covFrac) {
        $sim{$id1}{$id2} = $bits;
        push @sims, [$id1,$id2,$bits];
      }
    }
    close($fhsim) || die "Error reading $tmpPre.hits";
    unlink("$tmpPre.hits");
    print p("Found similarities, at above ${minIdentity}% identity and ${minCoverage}% coverage, for",
            scalar(keys %sim), "of these sequences"), "\n";

    # Try to find a small number of seed sequences that hit all of the members
    # Don't try to be optimal, just start with the seeds that have the most hits above the threshold
    my %nSims = map { $_ => scalar(keys %{ $sim{$_} }) } (keys %sim);
    my @seeds = sort { $nSims{$b} <=> $nSims{$a} || $a <=> $b } (keys %nSims);
    my %clust = (); # sequence to cluster, which is a hash of ids with the seed having a value of 1
    foreach my $seed (@seeds) {
       if (!exists $clust{$seed}) {
         my $chash = { $seed => 1 };
         $clust{$seed} = $chash;
         foreach my $id2 (keys %{ $sim{$seed} }) {
           if (!exists $clust{$id2}) {
             $clust{$id2} = $chash;
	     $chash->{$id2} = 0;
           }
         }
      }
    }

    # And, add the singletons
    foreach my $id (keys %seqs) {
      if (!exists $clust{$id}) {
        $clust{$id} = { $id => 1 };
      }
    }

    # Report the clusters
    # First identify the unique ones
    my %clustByIds = ();
    foreach my $chash (values %clust) {
      $clustByIds{ join(":::", sort keys %$chash) } = $chash;
    }
    my @clusters = values %clustByIds;
    my @singletons = grep { scalar(keys %$_) == 1 } @clusters;
    my $clustReport = "Found " . (scalar(@clusters) - scalar(@singletons)) . " clusters of similar sequences.";
    $clustReport .= " Another " . scalar(@singletons) . " sequences are not clustered" if @singletons > 0;
    print p($clustReport), "\n";

    my @clustBySize = sort { scalar(keys %$b) <=> scalar(keys %$a) } @clusters;
    my $nCluster = 0;
    foreach my $cluster (@clustBySize) {
      my @ids = sort keys %$cluster;
      if (@ids > 1) {
        $nCluster++;
	my ($seed) = grep $cluster->{$_}, @ids;
	die unless defined $seed;
        print h3("Cluster $nCluster");
        print small("The first sequence in each cluster is the seed.") if $nCluster == 1; 
	my @other = grep ! $cluster->{$_}, @ids;
        foreach my $id ($seed, @other) {
          print CompoundInfoToHtml($id, $curatedInfo{$id}), "\n";
        }
      }
    }
    if (@singletons > 0) {
      print h3("Singletons");
      my @singletonIds = map { (keys %{ $_ })[0] } @singletons;
      foreach my $id (sort @singletonIds) {
        print CompoundInfoToHtml($id, $curatedInfo{$id}), "\n";
      }
    }
  }
}
print end_html;
exit(0);

sub MatchRows($$$) {
  my ($curatedInfo, $query, $word) = @_;
  die "Searching for empty term"
    unless defined $query && $query ne "";
  my $regexp = quotemeta($query);
  $regexp =~ s/\\%/.*/g;
  my @match = grep { $_->{descs} =~ m/$regexp/i } @$curatedInfo;
  return $word ? CuratedWordMatch(\@match, $query, "descs") : \@match;
}

sub FaaToDb($$) {
  my ($faaIn, $db) = @_;
  my %seqs;
  tie %seqs, "DB_File", $db, O_CREAT|O_RDWR, 0666, $DB_HASH
    || die "Cannot write to file $db -- $!";
  open (my $fh, "<", $faaIn) || die "Cannot read $faaIn";
  my $state = {};
  while (my ($header, $seq) = ReadFastaEntry($fh, $state)) {
    $header =~ s/ .*//;
    $seqs{$header} = $seq;
  }
  close($fh) || die "Error reading $faaIn";
  untie %seqs;
  die "Error writing to file $db" unless NewerThan($db, $faaIn);
}

sub CompoundInfoToHtml($$) {
  my ($compoundId, $info) = @_;
  die unless $compoundId;
  die "No info for $compoundId" unless $info;
  my @ids = split /,/, $compoundId;
  my @descs = split /;; /, $info->{descs};
  die "Mismatched length of ids and descs" unless scalar(@ids) == scalar(@descs);
  my $len = $info->{length};
  die unless $len;
  my @pieces = ();
  for (my $i = 0; $i < @ids; $i++) {
    my $id = $ids[$i];
    my $desc = $descs[$i];
    my ($db, $protId) = split /::/, $id;
    die "Cannot parse id $id" unless $protId;
    my $gene = { 'db' => $db, 'protId' => $protId, 'desc' => $desc,
                 'protein_length' => $len,
                 'comment' => '', 'name' => '', id2 => '' };
    AddCuratedInfo($gene);
    push @pieces, GeneToHtmlLine($gene);
  }
  return p(join("<BR>", @pieces, small("$len amino acids"))); 
}


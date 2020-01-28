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
use pbweb qw{start_page};
use DB_File;
use URI::Escape;
use HTML::Entities;
use IO::Handle qw{autoflush};

sub MatchRows($$$);
sub FaaToDb($$);

my $maxHits = 1000;

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
my $usearch = "../bin/usearch";
foreach my $x ($fastacmd,$usearch) {
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
    print p("Found " . scalar(@$hits) . " characterized proteins with matching descriptions. This was reduced to the limit of $maxHits.");
    my @hits = splice(@$hits, 0, $maxHits);
    $hits = \@hits;
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
    print p("Running usearch's ublast"), "\n";
    my $idFrac = $minIdentity / 100;
    my $covFrac = $minCoverage / 100;
    my $uCmd = "$usearch -ublast $tmpPre.faa -db $tmpPre.faa -evalue 0.001 -id $idFrac -query_cov $covFrac -target_cov $covFrac -blast6out $tmpPre.hits -quiet";
    system($uCmd) == 0 || die "usearch failed -- $uCmd -- !";
    unlink("$tmpPre.faa");
    my %hit = (); # id => id => bits if hit (removed self hits)
    my @hits = (); # list of [id1, id2, bits]
    open (my $fhHit, "<", "$tmpPre.hits") || die "Cannot read $tmpPre.hits";
    while(my $line = <$fhHit>) {
      my ($id1, $id2, $identity, $alnlen, $mm, $gap, $qbeg, $qend, $sbeg, $send, $eval, $bits) = split /\t/, $line;
      unless ($id1 eq $id2 || exists $hit{$id1}{$id2}) {
        $hit{$id1}{$id2} = $bits;
        push @hits, [$id1,$id2,$bits];
      }
    }
    close($fhHit) || die "Error reading $tmpPre.hits";
    unlink("$tmpPre.hits");
    print p("Found similarities, at above ${minIdentity}% identity and ${minCoverage}% coverage, for",
            scalar(keys %hit), "of these sequences"), "\n";

    # Sort the hits by bits and the greedily form complete clusters
    my %clust = map { $_ => { $_ => 1 } } (keys %seqs); # id to hash of members
    @hits = sort { $b->[2] <=> $a->[2] } @hits;
    foreach my $row (@hits) {
      my ($id1,$id2) = @$row;
      if (!exists $clust{$id1}{$id2}) {
        # not already clustered together, should they be?
        # need to verify that all existing members of the cluster hit id2
        my $chash = $clust{$id1};
        my $ok = 1;
        foreach my $id3 (keys %$chash) {
          $ok = 0 if $id3 ne $id2 && !exists $hit{$id2}{$id3};
        }
        if ($ok) {
          $chash->{$id2} = 1;
          $clust{$id2} = $chash;
        }
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
    my $clustReport = "Found " . (scalar(@clusters) - scalar(@singletons)) . " fully-connected clusters of similar sequences.";
    $clustReport .= " Another " . scalar(@singletons) . " sequences are not clustered" if @singletons > 0;
    print p($clustReport), "\n";

    my @clustBySize = sort { scalar(keys %$b) <=> scalar(keys %$a) } @clusters;
    my $nCluster = 0;
    foreach my $cluster (@clustBySize) {
      my @ids = sort keys %$cluster;
      if (@ids > 1) {
        $nCluster++;
        print h3("Cluster $nCluster");
        foreach my $id (@ids) {
          print p($id, $curatedInfo{$id}{descs}, small("(" . $curatedInfo{$id}{length}, "a.a.)"));
        }
        print "\n";
      }
    }
    if (@singletons > 0) {
      print h3("Singletons");
      foreach my $cluster (@clustBySize) {
        my @ids = sort keys %$cluster;
        if (@ids == 1) {
          my ($id) = @ids;
          print p($id, $curatedInfo{$id}{descs}, small("(" . $curatedInfo{$id}{length}, "a.a.)"));
        }
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
  @match = splice(@match, 0, $maxHits) if @match > $maxHits;
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

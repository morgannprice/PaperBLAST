#!/usr/bin/perl -w

# All parameters are optional. Without parameters, shows a query input page.
#
# Parameters:
# set -- which database of curated proteins (defaults to gaps2, i.e. using ../tmp/path.gaps2/curated.faa*)
#
# Specify which proteins to cluster using a query:
#	query -- what term to search for
#	word -- report whole word matches only (like a perl boolean)
# or using a step definition:
#	path -- which pathway the step is in
#	step -- which step
# Alternatively, browse the pathways or steps, use
#	path=all to list the pathways, or
#	set path but not step to list the steps in a pathway
# Can also find similar but non-matching characterized proteins, using
#	path, step, and close=1
#
# Clustering options:
# identity -- %identity threshold for clustering, default 30
# coverage -- %coverage threshold for clustering, default 75

use strict;
use CGI qw(:standard Vars start_ul end_ul);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use Time::HiRes qw{gettimeofday};
use DBI;
use lib "../lib";
use pbutils qw{ReadTable NewerThan CuratedWordMatch ReadFastaEntry};
use pbweb qw{start_page AddCuratedInfo GeneToHtmlLine};
use Steps qw{ReadSteps};
use DB_File;
use URI::Escape;
use HTML::Entities;
use IO::Handle qw{autoflush};
use List::Util qw{min max};

sub IsHetero($);
sub MatchRows($$$);
sub FaaToDb($$);
# The first argument is a "compound" id which is a list of ids
sub CompoundInfoToHtml($$$);
sub TransporterSearch($$$);

my $maxHits = 250;
my $nCPU = 8;

my $minHmmCoverageGood = 0.8;

my $set = param('set') || "gaps2";
$set =~ m/^[a-zA-Z0-9._-]+$/ || die "Invalid set $set";

my $query = param('query') || "";
my $wordMode = param('word') || 0;
my $pathSpec = param('path') || "";
my $step = param('step') || "";
my $closeMode = param('close') || 0;
die "Can only use close if path and step is set"
  if $closeMode && !($pathSpec && $step);

my $minIdentity = param('identity');
$minIdentity = 30 unless defined $minIdentity && $minIdentity =~ m/^\d+$/;
$minIdentity = 100 if $minIdentity > 100;
$minIdentity = 0 if $minIdentity < 0;

my $minCoverage = param('coverage');
$minCoverage = 75 unless defined $minCoverage && $minCoverage =~ m/^\d+$/;
$minCoverage = 100 if $minCoverage > 100;
$minCoverage = 0 if $minCoverage < 0;

my $stepPath = "../gaps/$set";
die "Invalid set $set: no $stepPath directory" unless -d $stepPath;
my $queryPath = "../tmp/path.$set"; # intermediate files
die "Invalid set $set: no $queryPath directory" unless -d $queryPath;

my ($steps, $banner, $bannerURL);
my @pathInfo;
my %pathInfo;
if ($pathSpec ne "") {
  $steps = ReadSteps("../gaps/$set/$pathSpec.steps") if $pathSpec && $pathSpec ne "all";
  @pathInfo = ReadTable("$stepPath/$set.table", ["pathwayId","desc"]);
  %pathInfo = map { $_->{pathwayId} => $_ } @pathInfo;
  $banner = "GapMind for " . ($pathInfo{all}{desc} || "pathways");
  $bannerURL = "gapView.cgi";
}

my $curatedFaa = "$queryPath/curated.faa";
die "No such file: $curatedFaa" unless -e $curatedFaa;

my $hetFile = "$queryPath/hetero.tab";
die "No such file: $hetFile" unless -e $hetFile;
my @het = ReadTable($hetFile, ["db","protId","comment"]);
my %hetComment = map { $_->{db} . "::" . $_->{protId} => $_->{comment} } @het;

my $fastacmd = "../bin/blast/fastacmd";
my $blastall = "../bin/blast/blastall";
my $formatdb = "../bin/blast/formatdb";
my $usearch = "../bin/usearch";
foreach my $x ($fastacmd,$blastall,$formatdb,$usearch) {
  die "No such executable: $x" unless -x $x;
}

my $tmpPre = "/tmp/cluratedClusters.$$";

start_page('title' => $closeMode ? "Other Proteins Similar to $step"
           : "Clusters of Characterized Proteins",
           'banner' => $banner, 'bannerURL' => $bannerURL);
autoflush STDOUT 1; # show preliminary results

unless (NewerThan("$curatedFaa.db", $curatedFaa)) {
  print p("Reformatting the sequence database"), "\n";
  FaaToDb($curatedFaa, "$curatedFaa.db");
}

$query =~ s/^\s+//;
$query =~ s/\s+$//;
if ($query eq "" && $pathSpec eq "") {
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
    end_form,
    p(a("Or",
        a({-href => "curatedClusters.cgi?set=$set&path=all"}, "cluster proteins"),
        "for steps in", a({-href => "gapView.cgi?set=$set"}, "GapMind"))),
    end_html;
  exit(0);
}

# Fetch curated sequences
my @curatedInfo = ReadTable("$curatedFaa.info", ["ids","length","descs"]);
my %curatedInfo = map { $_->{ids} => $_ } @curatedInfo;
my %seqsAll;
tie %seqsAll, "DB_File", "$curatedFaa.db", O_RDONLY, 0666, $DB_HASH
  || die "Cannot open file $curatedFaa.db -- $!";
# Some sequences are in the query file, and not from curated.faa
my %uniprotSeq = ();
my @hitIds = ();

my %idToIds = ();
foreach my $ids (keys %curatedInfo) {
  foreach my $id (split /,/, $ids) {
    $idToIds{$id} = $ids;
  }
}

my %hitHmm = (); # hits from HMMs
my %hitHmmGood = (); # hits from HMMs with high coverage
# hits from other rules (except uniprots not in the curated database which are in %uniprotSeq)
my %hitRule = ();
my %ignore = (); # ids to ignore similarity to (used in $closeMode)

if ($query =~ m/^transporter:(.+)$/) {
  my $compoundSpec = $1;
  my $sqldb = "../data/litsearch.db";
  my $dbh = DBI->connect("dbi:SQLite:dbname=$sqldb","","",{ RaiseError => 1 }) || die $DBI::errstr;
  my $staticDir = "../static";
  my $chits = TransporterMatch($dbh, $staticDir, $compoundSpec);
  my %idsHit = ();
  foreach my $chit (@$chits) {
    my $id = $chit->{db} . "::" . $chit->{protId};
    if (!exists $idToIds{$id}) {
      print p("Warning: $id ($chit->{desc}) matched but is not in the curated table"), "\n";
    } else {
      $idsHit{ $idToIds{$id} } = 1;
    }
  }
  if (scalar(keys %idsHit) == 0) {
    print p("Sorry, no transporters were found. Please try a different query.");
  } else {
    print p("Found", scalar(@$chits), "transporters with", scalar(keys %idsHit), "different sequences");
  }
  @hitIds = sort keys %idsHit;
} elsif ($query ne "") { # find similar proteins
  my $wordStatement = $wordMode ? " as complete word(s)" : "";
  print p("Searching for", b(HTML::Entities::encode($query)), $wordStatement),
    p("Or try", a({-href => "curatedClusters.cgi?set=${set}"}, "another search")),
    "\n";
  my $hits = MatchRows(\@curatedInfo, $query, $wordMode);
  if (@$hits > $maxHits) {
    print p("Found over $maxHits proteins with matching descriptions. Please try a different query.");
    print end_html;
    exit(0);
  } elsif (@$hits == 0) {
    print p("Sorry, no hits were found");
  } else {
    my $nHetero = scalar(grep IsHetero($_), map { $_->{ids} } @$hits);
    print p("Found " . scalar(@$hits) . " characterized proteins with matching descriptions. $nHetero of these are heteromeric.");
  }
  @hitIds = map $_->{ids}, @$hits;
} elsif ($pathSpec eq "all") { # show list of pathways
  die if $step ne "";
  print p("Pathways for", $pathInfo{all}{desc}),
    start_ul();
  foreach my $pathObj (@pathInfo) {
    if ($pathObj->{pathwayId} ne "all") {
      print li(a({ -href => "curatedClusters.cgi?set=$set&path=".$pathObj->{pathwayId} },
                 $pathObj->{desc}));
    }
  }
  print end_ul(),
    p("Or",a({-href => "curatedClusters.cgi?set=$set"}, "search by keyword"));
} elsif ($pathSpec ne "" && $step ne "") { # find proteins for this step
  my $stepDesc = $steps->{steps}{$step}{desc} || die "Unknown step $step in gaps/$set/$pathSpec.steps";

  my @querycol = qw{step type query desc file sequence};
  my @queries = ReadTable("$queryPath/$pathSpec.query", \@querycol);
  @queries = grep { $_->{step} eq $step } @queries;
  die "Unknown step $step, not in the $pathSpec.query file" unless @queries > 0;
  print p($closeMode ? "Finding" : "Clustering",
          "the characterized proteins for", i($step), "($stepDesc) in $pathInfo{$pathSpec}{desc}");
  if ($closeMode) {
    print p("Or see",
            a({-href => "curatedClusters.cgi?set=$set&path=$pathSpec&step=$step"},
                        "clustering"),
            "for step $step");
  } else {
    print p("Or see",
            a({-href => "curatedClusters.cgi?set=$set&path=$pathSpec&step=$step&close=1"},
              "other proteins similar to $step"));
  }
  print p("Or see all steps for",
          a({-href => "curatedClusters.cgi?set=$set&path=$pathSpec"}, $pathInfo{$pathSpec}{desc}));

  foreach my $query (@queries) {
    if ($query->{type} eq "curated") {
      my $ids = $query->{query};
      die "Unknown identifier $ids in query for $step" unless exists $seqsAll{$ids};
      push @hitIds, $ids;
      $hitRule{$ids} = 1;
    } elsif ($query->{type} eq "hmm") {
      my $hmm = $query->{query};
      my $hmmFile = $queryPath. "/" . $query->{file};
      die "No such hmm file: $hmmFile" unless -e $hmmFile;
      my $hitsFile = "$hmmFile.curatedhits";
      unless (NewerThan($hitsFile, $hmmFile) && NewerThan($hitsFile, $curatedFaa)) {
        my $hmmsearch = "../bin/hmmsearch";
        die "No such executable: $hmmsearch" unless -x $hmmsearch;
        my $cmd = "$hmmsearch --cpu $nCPU --cut_tc -o /dev/null --domtblout $hitsFile.$$.tmp $hmmFile $curatedFaa";
        system($cmd) == 0 || die "Failed running $cmd -- $!";
        rename("$hitsFile.$$.tmp",$hitsFile) || die "Cannot rename $hitsFile.$$.tmp to $hitsFile";
      }
      open(my $fhHits, "<", $hitsFile) || die "Cannot read $hitsFile";
      while(my $line = <$fhHits>) {
        chomp $line;
        next if $line =~ m/^#/;
        my @F = split /\s+/, $line;
      my ($ids, undef, $hitlen, $hmmName, undef, $hmmLen, $seqeval, $seqscore, undef, undef, undef, undef, $domeval, $domscore, undef, $qbeg, $qend, $hitbeg, $hitend) = @F;
        die "Cannot parse hmmsearch line $line"
          unless defined $hitend && $hitend > 0 && $hmmLen > 0 && $ids ne "";
        die "Unknown id $ids reported by hmmsearch" unless exists $curatedInfo{$ids};
        my $highCoverage = ($qend-$qbeg+1)/$hmmLen >= $minHmmCoverageGood;
        $hitHmmGood{$ids} = 1 if $highCoverage;
        # if in close mode, ignore low-coverage HMM-only hits
        if ($highCoverage || ! $closeMode) {
          $hitHmm{$ids} = 1;
          push @hitIds, $ids;
        }
      }
      close($fhHits) || die "Error reading $hitsFile";
    } elsif ($query->{type} eq "uniprot") {
      my $id = "SwissProt::" . $query->{query};
      if (exists $idToIds{$id}) {
        print p("Warning: uniprot rule $query->{query} is already in curated"), "\n";
        push @hitIds, $idToIds{$id};
        $hitRule{ $idToIds{$id} } = 1;
      } else {
        $id =~ s/SwissProt/UniProt/;
        my $desc = $query->{desc};
        $desc =~ s/^(rec|sub)name: full=//i;
        $desc =~ s/ *[{][A-Z0-9:|.-]+[}] *//g;
        $desc =~ s/altname: full=//i;
        $desc =~ s/;? *$//;
        $curatedInfo{$id} = { 'ids' => $id, 'descs' => $desc,
                              'length' => length($query->{sequence}) };
        $uniprotSeq{$id} = $query->{sequence};
        push @hitIds, $id;
      }
    } elsif ($query->{type} eq "curated2") {
      ; # curated2 means not actually characterized, so skip that
    } elsif ($query->{type} eq "ignore") {
      $ignore{ $query->{query} } = 1;
    } else {
      die "Unknown query type $query->{type} for step $step in $queryPath/$pathSpec.query";
    }
  }
  print p("Sorry, no curated sequences were found matching any of", scalar(@queries),
          "rules for step $step") if @hitIds == 0;

  # Remove hits that match ignore (this should only happen for HMM hits)
  @hitIds = grep !exists $ignore{ $_ }, @hitIds;
} elsif ($pathSpec ne "") {
  die "Unknown pathway $pathSpec" unless exists $pathInfo{$pathSpec};
  my @steps = sort { $a->{i} <=> $b->{i} } values %{ $steps->{steps} };
  print h3("Steps in $pathInfo{$pathSpec}{desc}"),
    p("Cluster the charaterized proteins for a step:"), start_ul();
  foreach my $stepObj (@steps) {
    print li(a({ -href => "curatedClusters.cgi?set=$set&path=$pathSpec&step=$stepObj->{name}" }, 
            $stepObj->{name}) . ":", $stepObj->{desc});
  }
  print end_ul(),
    p("Or see all", a({-href => "curatedClusters.cgi?set=$set&path=all"}, "pathways"));
} else {
  die "Invalid mode for curatedClusters.cgi";
}

if (@hitIds == 0) {
  print end_html;
  exit(0);
}

# Fetch the hits' sequences
my %seqs = ();
foreach my $id (@hitIds) {
  if (exists $seqsAll{$id}) {
    $seqs{$id} = $seqsAll{$id};
  } elsif (exists $uniprotSeq{$id}) {
    $seqs{$id} = $uniprotSeq{$id};
  } else {
    die "No sequence for $id in $curatedFaa.db";
  }
}
my %hitsUniq = map { $_ => 1 } @hitIds;
die unless scalar(keys %seqs) == scalar(keys %hitsUniq);
print p("Fetched " . scalar(keys %seqs) . " sequences"), "\n";

open(my $fhFaa, ">", "$tmpPre.faa") || die "Error writing to $tmpPre.faa";
foreach my $id (sort keys %seqs) {
  print $fhFaa ">$id\n$seqs{$id}\n";
}
close($fhFaa) || die "Error writing to $tmpPre.faa";

if ($closeMode && $pathSpec && $step) {
  die "No such file: $curatedFaa.udb" unless -e "$curatedFaa.udb";
  my $minIdClose = 0.4;
  my $minCovClose = 0.7;
  print p("Running ublast to find other characterized proteins with",
          ($minIdClose*100)."%",
          "identity and",
          ($minCovClose*100)."%",
          "coverage"), "\n";
  my $cmd = "$usearch -ublast $tmpPre.faa -db $curatedFaa.udb -evalue 0.001 -blast6out $tmpPre.close -threads $nCPU -quiet -id $minIdClose -query_cov $minCovClose";
  system($cmd) == 0 || die "Error running\n$cmd\n-- $!";
  my @close = ();
  open (my $fhClose, "<", "$tmpPre.close") || die "Cannot read $tmpPre.close";
  while(my $line = <$fhClose>) {
    chomp $line;
    my ($query, $subject, $identity, $alen, $mm, $gap, $qbeg, $qend, $sbeg, $send, $eval, $bits) = split /\t/, $line;
    die "Cannot parse blast6out line $line" unless $bits > 0;
    push @close, { 'query' => $query, 'subject' => $subject, 'identity' => $identity,
                   'qbeg' => $qbeg, 'qend' => $qend, 'sbeg' => $sbeg, 'send' => $send,
                   'eval' => $eval, 'bits' => $bits };
    die "Unknown subject $subject" unless exists $curatedInfo{$subject};
  }
  close($fhClose) || die "Error reading $tmpPre.close";
  unlink("$tmpPre.faa");
  unlink("$tmpPre.close");
  @close = sort { $b->{bits} <=> $a->{bits} } @close;
  my %subjectSeen = ();
  @close = grep { my $subject = $_->{subject};
                  my $seen = exists $subjectSeen{$subject};
                  $subjectSeen{$subject} = 1;
                  ! $seen; } @close;
  my $nHitsAll = scalar(@close);
  @close = grep !exists $hitsUniq{ $_->{subject} }, @close;
  my $nPreIgnore = scalar(@close);
  @close = grep !exists $ignore{ $_->{subject} }, @close;
  print p("Found hits to $nPreIgnore  other characterized sequences.",
          "(Found $nHitsAll hits including self hits.)");
  print p("Removing ignored items reduces this to", scalar(@close), "close sequences.")
    if scalar(@close) < $nPreIgnore;
  # And also record hmm-only good hits
  my $hasSeqRule = keys(%uniprotSeq) > 0 || keys(%hitRule) > 0;
  my @hitHmmGoodOnly = ();
  if ($hasSeqRule) {
    @hitHmmGoodOnly = grep !exists $ignore{$_} && !exists $hitRule{$_}, keys %hitHmmGood;
    print p("There are also", scalar(@hitHmmGoodOnly), "HMM-only high-coverage hits.")
      if @hitHmmGoodOnly > 0;
  }
  if (@close > 0) {
    print h3("Close sequences");
    print p("(Sequences that are similar to these will not be high-confidence candidates for", i($step).".)");
  }
  foreach my $close (@close) {
    my $subject = $close->{subject};
    my $query = $close->{query};
    my $queryDesc = CompoundInfoToHtml($query, $curatedInfo{$query}, $seqs{$query});

    my $idString = int($close->{identity})."%";
    my $covString =   "Amino acids $close->{sbeg}:$close->{send}/$curatedInfo{$subject}{length}"
      . " are $idString identical to $close->{qbeg}:$close->{qend}/$curatedInfo{$query}{length}"
        . " of the $step protein";

    print p({ -style => 'margin-bottom: 0em;' },
            CompoundInfoToHtml($subject, $curatedInfo{$subject}, $seqsAll{$subject}),
            br(),
            a({-title => $covString}, "$idString identical to")),
           p({-style => 'margin-left: 5em; margin-top: 0em; font-size: 90%;'},
                CompoundInfoToHtml($query, $curatedInfo{$query}, $seqs{$query})),
           "\n";
  }

  if (@hitHmmGoodOnly > 0) {
    print h3("HMM-only sequences");
    print p("(Since these sequences' annotations are outside the definition for", i($step).",",
            "HMM hits that are over 40% similar to these sequences will be scored as moderate confidence.)");
    foreach my $hit (@hitHmmGoodOnly) {
      print p(CompoundInfoToHtml($hit, $curatedInfo{$hit}, $seqs{$hit}))."\n";
    }
  }
  print end_html;
  exit(0);
}

# else
# Clustering mode
print p("Running BLASTp"), "\n";
my $covFrac = $minCoverage / 100;
my $formatCmd = "$formatdb -i $tmpPre.faa -p T";
system($formatCmd) == 0 || die "formatdb failed -- $formatCmd -- $!";
my $blastCmd = qq{$blastall -p blastp -i $tmpPre.faa -d $tmpPre.faa -F "m S" -e 0.001 -m 8 -a $nCPU -o $tmpPre.hits >& /dev/null};
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
my @seeds = sort { $nSims{$b} <=> $nSims{$a} || $a cmp $b } (keys %nSims);
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

# Mark HMM only hits
my $hasSeqRule = keys(%uniprotSeq) > 0 || keys(%hitRule) > 0;
my $hmmColor = $hasSeqRule ? "red" : "black";
# mark hits from Hmm only
foreach my $ids (keys %hitHmm) {
  if (!exists $hitRule{$ids}) {
    $curatedInfo{$ids}{descs} .= exists $hitHmmGood{$ids} ?
      span({-style => "color: $hmmColor;"}, " (from HMM only)")
        : " (low-coverage HMM hit)";
  }
}

# First identify the unique ones
my %clustByIds = ();
foreach my $chash (values %clust) {
  $clustByIds{ join(":::", sort keys %$chash) } = $chash;
}
my @clusters = values %clustByIds;
my @singletons = grep { scalar(keys %$_) == 1 } @clusters;
my $clustReport = "Found " . (scalar(@clusters) - scalar(@singletons)) . " clusters of similar sequences.";
$clustReport .= " Another " . scalar(@singletons) . " sequences are not clustered." if @singletons > 0;
print p($clustReport), "\n";

my @clustBySize = sort { scalar(keys %$b) <=> scalar(keys %$a) } @clusters;
my $nCluster = 0;
foreach my $cluster (@clustBySize) {
  my @ids = sort keys %$cluster;
  if (@ids > 1) {
    $nCluster++;
    my ($seed) = grep $cluster->{$_}, @ids;
    die unless defined $seed;
    my $nHetero = scalar(grep IsHetero($_), @ids);
    my $sz = scalar(@ids);
    my @len = map $curatedInfo{$_}{length}, @ids;
    my @clusterHeader = ("Cluster $nCluster,", min(@len) . "-" .max(@len), "amino acids");
    if ($nHetero ==  $sz) {
      push @clusterHeader, "(heteromeric)";
    } elsif ($nHetero > 0) {
      push @clusterHeader, "($nHetero/$sz heteromeric)";
    } else {
      push @clusterHeader, "(not heteromeric)";
    }
    print h3(@clusterHeader);
    print small("The first sequence in each cluster is the seed.") if $nCluster == 1; 
    my @other = grep ! $cluster->{$_}, @ids;
    foreach my $id ($seed, @other) {
      print p(CompoundInfoToHtml($id, $curatedInfo{$id}, $seqs{$id})), "\n";
    }
  }
}
if (@singletons > 0) {
  my @singletonIds = map { (keys %{ $_ })[0] } @singletons;
  my $nHetero = scalar(grep IsHetero($_), @singletonIds);
  my $nSingle = scalar(@singletonIds);
  print h3("Singletons ($nHetero/$nSingle heteromeric)");
  foreach my $id (sort @singletonIds) {
    print p(CompoundInfoToHtml($id, $curatedInfo{$id}, $seqs{$id})), "\n";
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

sub CompoundInfoToHtml($$$) {
  my ($compoundId, $info, $seq) = @_;
  die unless $compoundId;
  die "No info for $compoundId" unless $info;
  die "no seq for $compoundId" unless $seq;
  my @ids = split /,/, $compoundId;
  die unless @ids > 0;
  my @descs = split /;; /, $info->{descs};
  die "Mismatched length of ids and descs" unless scalar(@ids) == scalar(@descs);
  my $len = $info->{length};
  die unless $len;
  my @genes;
  for (my $i = 0; $i < @ids; $i++) {
    my $id = $ids[$i];
    my $desc = $descs[$i];
    my ($db, $protId) = split /::/, $id;
    die "Cannot parse id $id" unless $protId;
    my $gene = { 'db' => $db, 'protId' => $protId, 'desc' => $desc,
                 'protein_length' => $len,
                 'comment' => '', 'name' => '', id2 => '' };
    AddCuratedInfo($gene);
    $gene->{HTML} = GeneToHtmlLine($gene);
    $gene->{HTML} .= " (" . i(a({ -title => $hetComment{$id} }, "heteromeric")) . ")"
      if exists $hetComment{$id};
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

sub IsHetero($) {
  my ($ids) = @_;
  my @ids = split /,/, $ids;
  foreach my $id (@ids) {
    return 1 if exists $hetComment{$id};
  }
  return 0;
}

sub TransporterMatch($$$) {
  my ($dbh, $staticDir, $compoundSpec) = @_;
  my @compoundList = split /:/, $compoundSpec;
  my %compoundList = map { $_ => 1, lc($_) => 1 } @compoundList;

  # metacyc files linking reactions to compounds and to protein(s)
  my $reactionCompoundsFile = "$staticDir/metacyc.reaction_compounds";
  my $reactionProteinsFile = "$staticDir/metacyc.reaction_links";

  my %rxnCompounds = (); # rxnId => compounds with side, coefficient, compartment, compoundId, and compoundName
  my %rxnLocation = ();
  open(my $fhc, "<", $reactionCompoundsFile) || die "Cannot read $reactionCompoundsFile";
  while (my $line = <$fhc>) {
    chomp $line;
    my @F = split /\t/, $line;
    my $rxnId = shift @F;
    my $rxnLoc = shift @F;
    $rxnLocation{$rxnId} = $rxnLoc;
    my @cmp = ();
    foreach my $f (@F) {
      my @pieces = split /:/, $f;
      my $side = shift @pieces;
      my $coefficient = shift @pieces;
      my $compartment = shift @pieces;
      my $compoundId = shift @pieces;
      my $compoundName = join(":", @pieces) || "";
      push @cmp, { 'side' => $side, 'coefficient' => $coefficient,
                   'compartment' => $compartment,
                   'compoundId' => $compoundId, 'compoundName' => $compoundName };
    }
    $rxnCompounds{$rxnId} = \@cmp;
  }
  close($fhc) || die "Error reading $reactionCompoundsFile";

  my %rxnProt = (); # rxnId => list of subunits of the form metacyc::id
  my %rxnName = ();
  open (my $fhp, "<", $reactionProteinsFile) || die "Cannot read $reactionProteinsFile";
  while (my $line = <$fhp>) {
    chomp $line;
    my @F = split /\t/, $line;
    my $rxnId = shift @F;
    my $rxnName = shift @F;
    push @{ $rxnProt{$rxnId} }, @F;
    $rxnName{$rxnId} = $rxnName;
  }
  close($fhp) || die "Error reading $reactionProteinsFile";

  my @out = (); # list of rows from CuratedGene
  my %metacycSaved = (); # protId => 1 for entries already in @out

  # Find the MetaCyc reactions that match the compound, either as compoundId or compoundName,
  # and determine if they are transport reactions, that is, the compound appears
  # with two different locations.
  # (For now ignore PTS)
  foreach my $rxnId (keys %rxnCompounds) {
    my $cmps = $rxnCompounds{$rxnId};
    my @cmpMatch = grep { exists $compoundList{$_->{compoundId}}
                            || exists $compoundList{$_->{compoundName}}
                              || exists $compoundList{lc($_->{compoundName})} } @$cmps;
    my %compartments = map { $_->{compartment} => 1 } @cmpMatch;
    my @compartments = grep { $_ ne "" } (keys %compartments);
    my %compartmentsAll = map { $_->{compartment} => 1 } @$cmps;
    my @compartmentsAll = grep { $_ ne "" } (keys %compartmentsAll);
    # The same compound on both sides, or, this compound on one side and something else on the other
    my $use = (@cmpMatch > 1 && @compartments > 1)
      || (@cmpMatch > 0 && @compartments > 0 && @compartmentsAll > 1);
    if ($use) {
      # Compound of interest in more than one compartment
      # So find all the proteins for this reaction
      my $protids = $rxnProt{$rxnId};
      foreach my $protid (@$protids) {
        my $id = $protid; $id =~ s/^metacyc:://;
        next if exists $metacycSaved{$id};
        $metacycSaved{$id} = 1;
        my $chit = $dbh->selectrow_hashref(qq{ SELECT * FROM CuratedGene WHERE db="metacyc" AND protId=? },
                                           {}, $id);
        if (defined $chit) {
          push @out, $chit;
        } else {
          print p("Skipped unknown metacyc id $id");
        }
      }
    }
  }

  my $tcdbs = $dbh->selectall_arrayref("SELECT * FROM CuratedGene WHERE db = 'TCDB';",
                                       { Slice => {} });
  foreach my $tcdb (@$tcdbs) {
    my @comments = split /_:::_/, $tcdb->{comment};
    my $keep = 0;
    foreach my $comment (@comments) {
      if ($comment =~ m/^SUBSTRATES: (.*)$/) {
        my @substrates = split /, /, $1;
        foreach my $substrate (@substrates) {
          $keep = 1 if exists $compoundList{lc($substrate)};
        }
      }
    }
    push @out, $tcdb if $keep;
  }

  # Search the description for words transport, transporter, exporter, porter, permease, import
  my $chits = $dbh->selectall_arrayref(
    qq{SELECT * FROM CuratedGene WHERE db IN ("SwissProt","CharProtDB","BRENDA","reanno","ecocyc")
       AND (desc LIKE "%transport%"
            OR desc LIKE "%porter%"
            OR desc LIKE "%import%"
            OR desc LIKE "%permease%"
            OR desc like "%PTS system%")},
    { Slice => {} });
  foreach my $chit (@$chits) {
    my $desc = $chit->{desc};
    my $keep = 0;
    foreach my $cmp (@compoundList) {
      my $pattern = quotemeta($cmp);
      # doing a word-based check this is tricky because substrates may be separated by "/", without spaces
      # or may appear as "glucose-binding"
      # or because another term like "6-phosphate" could be present as the next word.
      # That last issue is not handled.
      my $b = "[,/ ]"; # pattern for word boundary
      if ($desc =~ m!^$pattern$b!i # at beginning
          || $desc =~ m!$b$pattern$b!i # in middle
          || $desc =~ m!$b$pattern$!i # at end
          || $desc =~ m!^$pattern-(binding|specific)!i # at beginning
          || $desc =~ m!$b$pattern-(binding|specific)!i) { # in middle
        $keep = 1;
      }
    }
    push @out, $chit if $keep;
  }

  return \@out;
}

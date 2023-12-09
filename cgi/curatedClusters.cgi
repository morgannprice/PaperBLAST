#!/usr/bin/perl -w

# All parameters are optional. Without parameters, shows a query input page.
#
# Parameters:
# set -- which database of curated proteins (defaults to carbon, i.e. using
#	../tmp/path.carbon/curated.db and curated.faa.udb)
#
# Specify which proteins to cluster using a query:
#	query -- what term to search for, or transporter:compound:compound...
#	word -- report whole word matches only (like a perl boolean)
# or using a step definition:
#	path -- which pathway the step is in
#	step -- which step
# In cluster mode, can format the output as tab-delimited (format=tsv) or as rules (format=rules)
#	can also sort by organism (byorg)
# Alternatively, browse the pathways or steps, use
#	path=all to list the pathways, or
#       stepSpec to search for matching steps, or
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
use pbweb qw{start_page AddCuratedInfo GeneToHtmlLine DataForStepParts FormatStepPart LinkifyComment HMMToURL};
use Steps qw{ReadSteps};
use DB_File;
use URI::Escape;
use HTML::Entities;
use IO::Handle qw{autoflush};
use List::Util qw{min max};

sub GetHetero($$);
sub HeteroToText($);
sub MatchRows($$$);
sub CuratedToHtml($$); # curatedIds, sequence
sub TransporterSearch($$);
sub TSVPrint($$$);
sub TransporterMatch($$);
sub FormatClose($$$$);

my $maxHits = 250;
my $nCPU = 8;

my $minHmmCoverageGood = 0.8;

my $set = param('set') || "carbon";
$set =~ m/^[a-zA-Z0-9._-]+$/ || die "Invalid set $set";

my $query = param('query') || "";
my $wordMode = param('word') || 0;
my $pathSpec = param('path') || "";
my $step = param('step') || "";
my $closeMode = param('close') || 0;
my $format = param('format') || "";
die "Can only use close if path and step is set"
  if $closeMode && !($pathSpec && $step);
die "Format cannot be used with close mode" if $closeMode && $format ne "";
die "Unknown format spec"
  if $format ne "" && $format ne "tsv" && $format ne "rules";
my $byorg = param('byorg') || "";
die "byorg cannot be used with close mode" if $byorg && $closeMode;
die "byorg cannot be used with format" if $byorg && $format ne "";

my $minIdentity = param('identity');
$minIdentity = 30 unless defined $minIdentity && $minIdentity =~ m/^\d+$/;
$minIdentity = 100 if $minIdentity > 100;
$minIdentity = 0 if $minIdentity < 0;

my $minCoverage = param('coverage');
$minCoverage = 75 unless defined $minCoverage && $minCoverage =~ m/^\d+$/;
$minCoverage = 100 if $minCoverage > 100;
$minCoverage = 0 if $minCoverage < 0;

my $stepSpec = param('stepSpec') || "";
$stepSpec =~ s/^\s*//;
$stepSpec =~ s/\s*$//;

my $stepPath = "../gaps/$set";
die "Invalid set $set: no $stepPath directory" unless -d $stepPath;
my $queryPath = "../tmp/path.$set"; # intermediate files
die "Invalid set $set: no $queryPath directory" unless -d $queryPath;

my ($banner, $bannerURL);
my $pathInfo = (); # rows from Pathway table, in order
my %pathInfo; # pathwayId => row from Pathway table
my $dbhS; # database handle for steps.db; used only if pathSpec is set

if ($pathSpec ne "" || $stepSpec ne "") {
  my $stepsDb = "$queryPath/steps.db";
  $dbhS = DBI->connect("dbi:SQLite:dbname=$stepsDb","","",{ RaiseError => 1 }) || die $DBI::errstr;
  $pathInfo = $dbhS->selectall_arrayref("SELECT * FROM Pathway",
                                       { Slice => {} });
  %pathInfo = map { $_->{pathwayId} => $_ } @$pathInfo;
  $banner = "GapMind for " . ($pathInfo{all}{desc} || "pathways");
  $bannerURL = "gapView.cgi?set=${set}";
}

my $curatedFaa = "$queryPath/curated.faa";
my $curatedDb = "$queryPath/curated.db";
my $curatedUdb = "$queryPath/curated.faa.udb";
foreach my $file ($curatedFaa, $curatedDb, $curatedUdb) {
  die "No such file: $file" unless -e $file;
}

my $blastall = "../bin/blast/blastall";
my $formatdb = "../bin/blast/formatdb";
my $usearch = "../bin/usearch";
foreach my $x ($blastall,$formatdb,$usearch) {
  die "No such executable: $x" unless -x $x;
}

my $tmpPre = "/tmp/cluratedClusters.$$";

$query =~ s/^\s+//;
$query =~ s/\s+$//;
die "Can only use format with a query or with step" if $format ne "" && $query eq "" && ($pathSpec eq "" || $step eq "");

if ($format eq "") {
  start_page('title' => $closeMode ? "Other Characterized Proteins Similar to $step"
           : "Clusters of Characterized Proteins",
           'banner' => $banner, 'bannerURL' => $bannerURL);
} elsif ($format eq "tsv") {
  print "Content-Type:text/tab-separated-values\n";
  print "Content-Disposition: attachment; filename=clusters.tsv\n\n";
  print join("\t", qw{Sequence Cluster id id2 desc organism isHeteromeric})."\n";
} elsif ($format eq "rules") {
  print "Content-Type:text\n\n";
  print "# query: $query\n" if $query ne "";
  print "# pathway $pathSpec step $step\n" if $pathSpec ne "" && $step ne "";
} else {
  die "Unknown mode $format";
}

autoflush STDOUT 1; # show preliminary results

if ($stepSpec) {
  my $stepObjs = $dbhS->selectall_arrayref("SELECT * FROM Step WHERE stepId LIKE ?",
                                           { Slice => {} }, $stepSpec);
  if (@$stepObjs == 0) {
    print p("Sorry, no steps in $pathInfo{all}{desc} matched");
    # fall through to the all-pathways page
    $pathSpec = "all";
    $stepSpec = "";
  } elsif (@$stepObjs == 1) {
    # fall through to the 1-step page
    $pathSpec = $stepObjs->[0]{pathwayId};
    $step = $stepObjs->[0]{stepId};
  } else {
    # Show list of steps to click through to
    print p("Multiple steps matched:"), start_ul();
    foreach my $stepObj (@$stepObjs) {
      my $URL = join("&",
                     "curatedClusters.cgi?set=$set",
                     "path=$stepObj->{pathwayId}",
                     "step=$stepObj->{stepId}",
                     "identity=$minIdentity",
                     "coverage=$minCoverage");
      print li(a({ -href => $URL }, $stepObj->{stepId}),
               "from", $pathInfo{ $stepObj->{pathwayId} }{desc});
    }
    print end_ul(),
      p("Or see", a({ -href => "curatedClusters.cgi?set=$set&path=all"}, "all pathways")),
      end_html;
    exit(0);
  }
}

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
    p({-style => "margin-left: 3em;" },
      checkbox(-name => "byorg", -label => "Sort by organism?", -checked => 0)),
    p(submit(-name => 'Search')),
    end_form,
    p("Or cluster proteins for steps in GapMind for",
      a({-href => "curatedClusters.cgi?set=aa&path=all"}, "amino acid biosynthesis"),
      "or",
      a({-href => "curatedClusters.cgi?set=carbon&path=all"}, "carbon catabolism")),
    end_html;
  exit(0);
}

my $dbhC = DBI->connect("dbi:SQLite:dbname=$curatedDb","","",{ RaiseError => 1 }) || die $DBI::errstr;

# The list of matching sequences, mostly from CuratedInfo
# In step mode, some entries are from uniprot, with
#   curatedIds = UniProt::protId and the seq field set.
#	(curated2 is ignored)
my @hits = (); # Hits from matching curated description

# Remember which hits were from HMMs
my %hitHmm = (); # curatedIds => 1 if found by an HMM
my %hitHmmGood = (); # the same, but only for HMM hits with high coverage
# hits from other rules (except uniprots not in the curated database which are in %uniprotSeq)
my %uniprotSeq = ();
my %hitPart = (); # curatedIds => 1 if matches a (non-HMM) part
my %ignore = (); # ids to ignore similarity to (used in $closeMode)

if ($query =~ m/^transporter:(.+)$/) {
  my $compoundSpec = $1;
  my $queryShow = $compoundSpec; $queryShow =~ s!:! / !g;
  print p("Searching for transporters for", HTML::Entities::encode($queryShow)),
    "\n" unless $format;
  my $chits = TransporterMatch($dbhC, $compoundSpec);
  if (@$chits == 0) {
    print p("Sorry, no transporters were found. Please try a different query.")
      unless $format;
  } else {
    print p("Found", scalar(@$chits), "transporters")
      unless $format;
  }
  @hits = @$chits;
} elsif ($query ne "") { # find similar proteins
  my $wordStatement = $wordMode ? " as complete word(s)" : "";
  print p("Searching for", b(HTML::Entities::encode($query)), $wordStatement),
    p("Or try", a({-href => "curatedClusters.cgi?set=${set}"}, "another search")),
    "\n"
      unless $format;
  my $chits = MatchRows($dbhC, $query, $wordMode);
  if (@$chits > $maxHits) {
    print p("Found over $maxHits proteins with matching descriptions. Please try a different query.");
    print end_html;
    exit(0);
  } elsif (@$chits == 0) {
    print p("Sorry, no hits were found")
      unless $format;
  } else {
    my $nHetero = scalar(grep defined GetHetero($dbhC, $_), map { $_->{curatedIds} } @$chits);
    print p("Found " . scalar(@$chits) . " characterized proteins with matching descriptions. $nHetero of these are heteromeric.")
      unless $format;
  }
  @hits = @$chits;
} elsif ($pathSpec eq "all") { # show list of pathways
  die if $step ne "";
  print
    h3("Cluster curated proteins for", $pathInfo{all}{desc}),
    start_form(-method => 'get', -action => 'curatedClusters.cgi'),
    hidden(-name => 'set', -value => $set, -override => 1),
    p("Search for step:", textfield(-name => 'stepSpec', -value => '', -size => 25, -maxLength => 200),
      small("&nbsp;Use % as a wild card character")),
    p("Cluster at",
      textfield(-name => "identity", -value => $minIdentity, -size => 3, -maxlength => 3),
      "%identity and",
      textfield(-name => "coverage", -value => $minCoverage, -size => 3, -maxlemgth => 3),
      "%coverage"),
    p(submit(-name => 'Search')),
    end_form,
    p("Or browse pathways:"),
    start_ul();
  foreach my $pathObj (@$pathInfo) {
    if ($pathObj->{pathwayId} ne "all") {
      print li(a({ -href => "curatedClusters.cgi?set=$set&path=".$pathObj->{pathwayId} },
                 $pathObj->{desc}));
    }
  }
  my $otherSet = $set eq "aa" ? "carbon" : "aa";
  my $otherSetDesc = $otherSet eq "aa" ? "amino acid biosynthesis" : "carbon catabolism";
  print end_ul(),
    p("Or", a({-href => "curatedClusters.cgi?set=$otherSet&path=all"},
              "cluster proteins for $otherSetDesc")),
    p("Or", a({-href => "curatedClusters.cgi?set=$set"}, "search for curated proteins by keyword"));
} elsif ($pathSpec ne "" && $step ne "") { # find proteins for this step
  my ($stepDesc) = $dbhS->selectrow_array("SELECT desc FROM Step WHERE pathwayId = ? AND stepId = ?",
                                          {}, $pathSpec, $step);
  die "Unknown step $step in pathway $pathSpec" unless defined $stepDesc;
  my $stepParts = $dbhS->selectall_arrayref("SELECT * from StepPart WHERE pathwayId = ? AND stepId = ?",
                                            { Slice => {} }, $pathSpec, $step);
  my $stepPartData = DataForStepParts($dbhS, $pathSpec, $step);
  my $queries = $dbhS->selectall_arrayref("SELECT * from StepQuery WHERE pathwayId = ? AND stepId = ?",
                                          { Slice => {} }, $pathSpec, $step);
  print p($closeMode ? "Finding" : "Clustering",
          "the characterized proteins for", i($step), "($stepDesc) in $pathInfo{$pathSpec}{desc}")
    unless $format;
  if ($closeMode) {
    print p("Or see",
            a({-href => "curatedClusters.cgi?set=$set&path=$pathSpec&step=$step"},
                        "clustering"),
            "for step $step");
  } else {
    print p("Or see",
            a({-href => "curatedClusters.cgi?set=$set&path=$pathSpec&step=$step&close=1"},
              "other characterized proteins similar to $step")) unless $format;
  }
  print p("Or see all steps for",
          a({-href => "curatedClusters.cgi?set=$set&path=$pathSpec"}, $pathInfo{$pathSpec}{desc})),
        p("Or",
          a({-href => "curatedClusters.cgi?set=$set"}, "cluster curated proteins matching a keyword"))
    unless $format;

  # Show step description
  if ($format eq "") {
    print p(b("Definition of", i($step))) if $format eq "";
    print start_ul();
    foreach my $row (@$stepParts) {
      print li(FormatStepPart($stepPartData, $row, $set, "", undef));
    }
    my $stepsObj = Steps::ReadSteps("../gaps/$set/$pathSpec.steps"); # for the comment
    my $stepObj = $stepsObj->{steps}{$step} || die;
    print li("Comment:", LinkifyComment($stepObj->{comment}))
      if $stepObj->{comment} ne "";
    print end_ul(), "\n";
  }

  foreach my $query (@$queries) {
    if ($query->{queryType} eq "curated") {
      my $curatedIds = $query->{curatedIds};
      my $row = $dbhC->selectrow_hashref("SELECT * from CuratedInfo WHERE curatedIds = ?",
                                         {}, $curatedIds);
      die "Unknown curatedIds $curatedIds" unless defined $row;
      push @hits, $row;
      $hitPart{$curatedIds} = 1;
    } elsif ($query->{queryType} eq "hmm") {
      my $hmm = $query->{hmmId};
      my $hmmFile = $queryPath. "/" . $query->{hmmFileName};
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
        my ($curatedIds, undef, $hitlen, $hmmName, undef, $hmmLen, $seqeval, $seqscore, undef, undef, undef, undef, $domeval, $domscore, undef, $qbeg, $qend, $hitbeg, $hitend) = @F;
        die "Cannot parse hmmsearch line $line"
          unless defined $hitend && $hitend > 0 && $hmmLen > 0 && $curatedIds ne "";
        my $row = $dbhC->selectrow_hashref("SELECT * from CuratedInfo WHERE curatedIds = ?",
                                          {}, $curatedIds);
        die "Unknown id $curatedIds reported by hmmsearch" unless defined $row;
        my $highCoverage = ($qend-$qbeg+1)/$hmmLen >= $minHmmCoverageGood;
        $hitHmmGood{$curatedIds} = 1 if $highCoverage;
        # if in close mode, ignore low-coverage HMM-only hits
        if ($highCoverage || ! $closeMode) {
          $hitHmm{$curatedIds} = 1;
          push @hits, $row;
        }
      }
      close($fhHits) || die "Error reading $hitsFile";
    } elsif ($query->{queryType} eq "uniprot" || $query->{queryType} eq "predicted") {
      my $id = "UniProt::" . $query->{uniprotId};
      my $desc = $query->{desc};
      $desc =~ s/^(rec|sub)name: full=//i;
      $desc =~ s/ *[{][A-Z0-9:|.-]+[}] *//g;
      $desc =~ s/altname: full=//i;
      $desc =~ s/;? *$//;
      push @hits, { 'curatedIds' => $id, 'descs' => $desc,
                    'seqLength' => length($query->{seq}) };
      $uniprotSeq{$id} = $query->{seq};
    } elsif ($query->{queryType} eq "curated2") {
      ; # curated2 means not actually characterized, so skip that
    } elsif ($query->{queryType} eq "ignore") {
      $ignore{ $query->{curatedIds} } = 1;
    } else {
      die "Unknown query type $query->{queryType} for step $step in $queryPath/$pathSpec.query";
    }
  }
  print p("Sorry, no curated sequences were found matching any of", scalar(@$queries),
          "queries for step $step") if @hits == 0;

  # Remove hits that match ignore (this should only happen for HMM hits)
  @hits = grep { !exists $ignore{ $_->{curatedIds} } } @hits;
} elsif ($pathSpec ne "") {
  die "Unknown pathway $pathSpec" unless exists $pathInfo{$pathSpec};
  my $steps = $dbhS->selectall_arrayref("SELECT * from Step WHERE pathwayId = ?",
                                        { Slice => {} }, $pathSpec);
  print h3("Steps in $pathInfo{$pathSpec}{desc}"),
    p("Cluster the charaterized proteins for a step:"), start_ul();
  foreach my $stepObj (@$steps) {
    print li(a({ -href => "curatedClusters.cgi?set=$set&path=$pathSpec&step=$stepObj->{stepId}" },
            $stepObj->{stepId}) . ":", $stepObj->{desc});
  }
  print end_ul(),
    p("Or see an",
      a({ -href => "gapView.cgi?set=$set&orgs=orgsFit&path=$pathSpec&showdef=1" },
        "overview of $pathInfo{$pathSpec}{desc}")),
    p("Or see",
      a({-href => "curatedClusters.cgi?set=$set&path=all"}, "all pathways"));
} else {
  die "Invalid mode for curatedClusters.cgi";
}

if (@hits == 0) {
  print end_html;
  exit(0);
}

# Fetch the hits' sequences
my %seqs = ();
foreach my $hit (@hits) {
  my $id = $hit->{curatedIds};
  if (exists $uniprotSeq{$id}) {
    $seqs{$id} = $uniprotSeq{$id};
  } else {
    my ($seq) = $dbhC->selectrow_array("SELECT seq FROM CuratedSeq WHERE curatedIds = ?",
                                      {}, $id);
    die "curatedIds $id has no sequence" unless defined $seq;
    $seqs{$id} = $seq;
  }
}
my %hitInfo = map { $_->{curatedIds} => $_ } @hits;
die "Seqs " . join(" ", keys %seqs) . "\nhits " . join(" ", keys %hitInfo)
  unless scalar(keys %seqs) == scalar(keys %hitInfo);
print p("Fetched " . scalar(keys %seqs) . " sequences"), "\n"
  unless $format;

open(my $fhFaa, ">", "$tmpPre.faa") || die "Error writing to $tmpPre.faa";
foreach my $id (sort keys %seqs) {
  print $fhFaa ">$id\n$seqs{$id}\n";
}
close($fhFaa) || die "Error writing to $tmpPre.faa";

if ($closeMode && $pathSpec && $step) {
  my $minIdClose = 0.4;
  my $minCovClose = 0.7;
  print p("Running ublast to find other characterized proteins with",
          ($minIdClose*100)."%",
          "identity and",
          ($minCovClose*100)."%",
          "coverage"), "\n";
  my $cmd = "$usearch -ublast $tmpPre.faa -db $curatedUdb -evalue 0.001 -blast6out $tmpPre.close -threads $nCPU -quiet -id $minIdClose -query_cov $minCovClose";
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
  }
  close($fhClose) || die "Error reading $tmpPre.close";
  unlink("$tmpPre.faa");
  unlink("$tmpPre.close");

  # Add queryInfo, score to the close objects
  foreach my $close (@close) {
    my $query = $close->{query};
    my $queryInfo = $hitInfo{$query} || die;
    $close->{queryInfo} = $queryInfo;
    my $cov = ($close->{qend} - $close->{qbeg} + 1) / $queryInfo->{seqLength};
    die "Invalid coverage $cov for $query subject $close->{subject}" unless $cov > 0 && $cov <= 1;
    $close->{score} = $cov * $close->{identity};
  }

  @close = sort { $b->{score} <=> $a->{score} } @close;
  my %subjectSeen = ();
  @close = grep { my $subject = $_->{subject};
                  my $seen = exists $subjectSeen{$subject};
                  $subjectSeen{$subject} = 1;
                  ! $seen; } @close;
  my $nHitsAll = scalar(@close);
  @close = grep !exists $hitInfo{ $_->{subject} }, @close;
  my $nPreIgnore = scalar(@close);
  my @ignore = grep exists $ignore{ $_->{subject} }, @close;
  @close = grep !exists $ignore{ $_->{subject} }, @close;
  print p("Found hits to $nPreIgnore  other characterized sequences.",
          "(Found $nHitsAll hits including self hits.)");
  print p("Removing ignored items reduces this to", scalar(@close), "close sequences.")
    if scalar(@close) < $nPreIgnore;
  # And also record hmm-only good hits
  my $hasSeqRule = keys(%uniprotSeq) > 0 || keys(%hitPart) > 0;
  my @hitHmmGoodOnly = ();
  if ($hasSeqRule) {
    @hitHmmGoodOnly = grep !exists $ignore{$_} && !exists $hitPart{$_}, keys %hitHmmGood;
    print p("There are also", scalar(@hitHmmGoodOnly), "HMM-only high-coverage hits.")
      if @hitHmmGoodOnly > 0;
  }
  if (@close > 0) {
    print h3("Close sequences");
    print p("(Sequences that are similar to these will not be high-confidence candidates for", i($step).".)");
  }
  foreach my $close (@close) {
    print FormatClose($dbhC, $close, \%seqs, "black");

  }
  if (@hitHmmGoodOnly > 0) {
    print h3("HMM-only sequences");
    print p("(Since these sequences' annotations are outside the definition for", i($step).",",
            "HMM hits that are over 40% similar to these sequences will be scored as moderate confidence.)");
    foreach my $hit (@hitHmmGoodOnly) {
      print p(CuratedToHtml($hitInfo{$hit}, $seqs{$hit}))."\n";
    }
  }

  if (@ignore > 0) {
    print h3("Close but ignored sequences");
    print p("(Sequences that are similar to these will still be high-confidence candidates for", i($step).".)");
    foreach my $close (@ignore) {
      print FormatClose($dbhC, $close, \%seqs, "darkgrey");
    }
  }

  print end_html;
  exit(0);
}

# else
# Clustering mode
print p("Running BLASTp"), "\n" unless $format;
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
        scalar(keys %sim), "of these sequences"), "\n"
  unless $format;

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
my $hasSeqRule = keys(%uniprotSeq) > 0 || keys(%hitPart) > 0;
my $hmmColor = $hasSeqRule ? "red" : "black";
# mark hits from Hmm only
foreach my $ids (keys %hitHmm) {
  if (!exists $hitPart{$ids}) {
    $hitInfo{$ids}{descs} .= exists $hitHmmGood{$ids} ?
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
unless ($format) {
  my $clustReport = "Found " . (scalar(@clusters) - scalar(@singletons)) . " clusters of similar sequences.";
  $clustReport .= " Another " . scalar(@singletons) . " sequences are not clustered." if @singletons > 0;
  my $baseURL = "curatedClusters.cgi?set=$set&path=$pathSpec&step=$step&query=$query&word=$wordMode&identity=$minIdentity&coverage=$minCoverage";
  my $viewBy;
  if ($byorg) {
    $viewBy = a({-href => $baseURL}, "by cluster");
  } else {
    $viewBy = a({-href => $baseURL."&byorg=1"}, "by organism");
  }
  print p($clustReport,
          "Download as", a({-href => "$baseURL&format=tsv"}, "table"),
          "or as",
	  a({-href => "$baseURL&format=rules"}, "draft rules"),
          "or view", $viewBy);
}

my @clustBySize = sort { scalar(keys %$b) <=> scalar(keys %$a) } @clusters;
my @singletonIds = sort map { (keys %{ $_ })[0] } @singletons;

my $nCluster = 0;
my $nSingle = 0;
my %idsToCluster = (); # ids to cluster name or singleton name
my %clusterNameToIds = (); # for clusters, cluster name to list of members
foreach my $cluster (@clustBySize) {
  my @ids = sort keys %$cluster;
  if (@ids > 1) {
    $nCluster++;
    my $name = "Cluster $nCluster";
    foreach my $ids (@ids) {
      $idsToCluster{$ids} = $name;
    }
    $clusterNameToIds{$name} = \@ids;
  }
}
foreach my $ids (@singletonIds) {
  $nSingle++;
  $idsToCluster{ $ids } = "Singleton $nSingle";
}

if ($byorg) { # show by organism
  # simplify organisms (so that they will hopefully be consistent across databases)
  my %orgIds = (); # organism => ids => 1
  foreach my $ids (keys %idsToCluster) {
    my $info = $hitInfo{$ids} || die;
    my @orgs = split /;; /, $info->{orgs};
    my %orgU = ();
    push @orgs, "Organism unknown" if @orgs == 0;
    foreach my $org (@orgs) {
      $org =~ s/^candidatus //i;
      $org =~ s/^([a-zA-Z]+) sp./$1/;
      my @words = split / /, $org;
      @words = splice(@words, 0, 2) if @words > 2;
      $org = join(" ", @words);
      $orgIds{$org}{$ids} = 1;
    }
  }

  foreach my $org (sort keys %orgIds) {
    my @ids = sort keys %{ $orgIds{$org} };
    print h3($org),"\n";
    foreach my $ids (@ids) {
      my $clusterName = $idsToCluster{$ids};
      my $clusterShow = $clusterName;
      if (exists $clusterNameToIds{$clusterName}) {
        my @idsLeft = grep { $_ ne $ids } @{ $clusterNameToIds{$clusterName} };
        my $idsLeftSpec = join("&", map "ids=$_", @idsLeft);
        my $URL = "curatedSim.cgi?set=$set&path=$pathSpec&ids=$ids&$idsLeftSpec";
        $clusterShow = a({ -href => $URL }, $clusterShow);
      }
      print p(CuratedToHtml($hitInfo{$ids}, $seqs{$ids}),
              "($clusterShow)"), "\n";
    }
  }
} else { # !$byorg, or, show by cluster
  foreach my $cluster (@clustBySize) {
    my @ids = sort keys %$cluster;
    if (@ids > 1) {
      my ($seed) = grep $cluster->{$_}, @ids;
      my $clusterId = $idsToCluster{$seed} || die;
      die unless defined $seed;
      my $nHetero = scalar(grep defined GetHetero($dbhC, $_), @ids);
      my $sz = scalar(@ids);
      my @len = map $hitInfo{$_}{seqLength}, @ids;
      my @clusterHeader = ($clusterId, min(@len) . "-" .max(@len), "amino acids");
      if ($nHetero ==  $sz) {
        push @clusterHeader, "(heteromeric)";
      } elsif ($nHetero > 0) {
        push @clusterHeader, "($nHetero/$sz heteromeric)";
      } else {
        push @clusterHeader, "(not heteromeric)";
      }
      print h3(@clusterHeader) unless $format;
      print "\n# " . join(" ", @clusterHeader) . "\n" if $format eq "rules";
      print small("The first sequence in each cluster is the seed.")
        if $nCluster == 1 && $format eq ""; 
      my @other = grep ! $cluster->{$_}, @ids;
      foreach my $id ($seed, @other) {
        if ($format eq "") {
          my @idsLeft = grep { $_ ne $id } @ids;
          my $idsLeftSpec = join("&", map "ids=$_", @idsLeft);
          print p(CuratedToHtml($hitInfo{$id}, $seqs{$id})
                  . ", "
                  . small(a({ -href => "curatedSim.cgi?set=$set&path=$pathSpec&ids=$id&$idsLeftSpec" },
                    "Compare to cluster")));
        } elsif ($format eq "rules") {
          my $het = HeteroToText(GetHetero($dbhC, $id));
          print "# $id $hitInfo{$id}{id2s} $hitInfo{$id}{descs} $hitInfo{$id}{orgs}$het\n";
        } elsif ($format eq "tsv") {
          $clusterId =~ s/ //;
          TSVPrint($clusterId, $id, $hitInfo{$id});
        }
      }
      if ($format eq "rules") {
        my $desc = $hitInfo{$seed}{descs} || "No description";
        $desc =~ s/ ;;.*//;
        $clusterId =~ s/ //;
        print join("\t", $clusterId,
                   $desc,
                   map { my $short = $_; $short =~ s/,.*//; "curated:$short" } ($seed, @other))."\n";
      }
    } # end if @ids > 1
  }
  if (@singletons > 0) {
    my $nHetero = scalar(grep defined GetHetero($dbhC, $_), @singletonIds);
    my $nSingle = scalar(@singletonIds);
    print h3("Singletons ($nHetero/$nSingle heteromeric)") unless $format;
    foreach my $id (@singletonIds) {
      my $singleId = $idsToCluster{$id} || die;
      if ($format eq "") {
        my $show = CuratedToHtml($hitInfo{$id}, $seqs{$id});
        if (exists $sim{$id}) {
          my @simto = sort keys %{ $sim{$id} };
          my $sim = $simto[0] || die;
          my $simCluster = $idsToCluster{$sim};
          $sim =~ s/,.*//; # keep 1st id only
          $show .= br() . small("(similar to $sim from ${simCluster},",
                                "but similarity to seed sequence is below thresholds)");
        }
        print p($show), "\n";
      } elsif ($format eq "rules") {
        my $het = HeteroToText(GetHetero($dbhC, $id));
        my $id2s = $hitInfo{$id}{id2s} || "";
        my $descs = $hitInfo{$id}{desc} || "";
        my $orgs = $hitInfo{$id}{orgs} || "";
        print "\n# $singleId $id $id2s $descs $orgs$het\n";
        my $short = $id; $short =~ s/,.*//;
        my $desc = $hitInfo{$id}{descs};
        $desc =~ s/;; .*//;
        $singleId =~ s/ //;
        print join("\t", $singleId, $desc, "curated:$short")."\n";
      } elsif ($format eq "tsv") {
        $singleId =~ s/ //;
        TSVPrint($singleId, $id, $hitInfo{$id});
      }
    }
  }
} # end else on $byorg
print end_html unless $format;
exit(0);

sub MatchRows($$$) {
  my ($dbhC, $query, $word) = @_;
  die "Searching for empty term"
    unless defined $query && $query ne "";
  if ($word && $query =~ m/^[0-9][.][0-9-]+[.][0-9-]+[.][A-Za-z]?[0-9-]*$/) {
    # use the EC table
    return $dbhC->selectall_arrayref(qq{ SELECT * from ECToCurated
                                         JOIN CuratedInfo USING (curatedIds)
                                         WHERE ec = ? },
                                     { Slice => {} }, $query);
  }
  #else
  my $chits = $dbhC->selectall_arrayref("SELECT * from CuratedInfo WHERE descs LIKE ?",
                                       { Slice => {} }, "%" . $query . "%");
  if ($word) {
    my @out = ();
    foreach my $chit (@$chits) {
      my @descs = split /;; /, $chit->{descs};
      my @list = map { { 'desc' => $_ } } @descs;
      my $matches = CuratedWordMatch(\@list, $query, 'desc');
      push @out, $chit if @$matches > 0;
    }
    return \@out;
  }
  #else
  return $chits;
}

sub CuratedToHtml($$) {
  my ($info, $seq) = @_;
  die "Undefined info" unless defined $info;
  die "No curatedIds for info" unless $info->{curatedIds};
  my $curatedIds = $info->{curatedIds};
  die "no seq for $curatedIds" unless $seq;
  my @ids = split /,/, $curatedIds;
  die unless @ids > 0;
  my @descs = split /;; /, $info->{descs};
  my @orgs = split /;; /, $info->{orgs} if defined $info->{orgs} && $info->{orgs} ne "";
  my @id2s = split /;; /, $info->{id2s} if defined $info->{id2s} && $info->{id2s} ne "";

  die "Mismatched length of ids and descs" unless scalar(@ids) == scalar(@descs);
  my $len = $info->{seqLength};
  die unless $len;
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
  my $pfamHits = $dbhC->selectall_arrayref("SELECT * FROM CuratedPFam WHERE curatedIds == ?",
                                           { Slice => {} }, $curatedIds);
  if (@$pfamHits > 0) {
    my @rows = sort { $a->{seqFrom} <=> $b->{seqFrom} } @$pfamHits;
    my @pfams = ();
    foreach my $row (@rows) {
      my $pfam = $row->{hmmAcc}; $pfam =~ s/[.]\d+$//;
      die $pfam unless $pfam =~ m/^PF\d+$/;
      my $cov = ($row->{hmmTo} - $row->{hmmFrom} + 1) / $row->{hmmLen};
      my $covPercent = int(0.5 + 100 * $cov);
      push @pfams, a({ -href => HMMToURL($pfam),
                       -title => "amino acids $row->{seqFrom}:$row->{seqTo} (${covPercent}% coverage of $row->{hmmAcc}, $row->{bits} bits)" },
                     $row->{hmmName});
    }
    push @pieces, small("PFams:", join(", ", @pfams));
  }
  my $hetComment = "";
  $hetComment = a({ -title => GetHetero($dbhC, $curatedIds) }, b("Heteromeric")) . ", "
    if defined GetHetero($dbhC, $curatedIds);
  return join("<BR>", @pieces,
                small($hetComment . $len, "amino acids: ",
                      join(", ", @links)));
}

my %heteroCache;

# Returns undef if not heteromeric; otherwise returns the comment (but it is usually empty)
sub GetHetero($$) {
  my ($dbhC, $curatedIds) = @_;
  unless (exists $heteroCache{$curatedIds}) {
    my $row = $dbhC->selectrow_arrayref("SELECT comment FROM Hetero WHERE curatedIds = ?",
                                        {}, $curatedIds);
    $heteroCache{$curatedIds} = defined $row ? $row->[0] : undef;
  }
  return $heteroCache{$curatedIds};
}

sub TransporterMatch($$) {
  my ($dbhC, $compoundSpec) = @_;
  my @compoundList = split /:/, $compoundSpec;
  my %compoundList = map { $_ => 1, lc($_) => 1 } @compoundList;

  # Find MetaCyc reactions that match the compound (as id or as compound name)
  # and have the compound in two different compartments.
  # (Punt on PTS.)
  my %metacycSaved = ();
  foreach my $compound (@compoundList) {
    my $cmpRxns = $dbhC->selectall_arrayref("SELECT * FROM CompoundInReaction WHERE cmpId = ? OR cmpDesc LIKE ?",
                                         { Slice => {} }, $compound, $compound);
    my %rxnAt = (); # rxnId => compartment => 1 for this compound
    foreach my $row (@$cmpRxns) {
      $rxnAt{ $row->{rxnId} }{ $row->{compartment} } = 1;
    }
    while (my ($rxnId, $compartments) = each %rxnAt) {
      next unless keys(%$compartments) > 1;
      my $curatedIdsList = $dbhC->selectcol_arrayref("SELECT * from EnzymeForReaction WHERE rxnId = ?",
                                                    {}, $rxnId);
      foreach my $curatedIds (@$curatedIdsList) {
        $metacycSaved{$curatedIds} = 1;
      }
    }
  }

  # Find entries in TCDB
  my %tcdbSaved = ();
  my $tcdbs = $dbhC->selectall_arrayref("SELECT * FROM TransporterSubstrate",
                                       { Slice => {} });
  foreach my $tcdb (@$tcdbs) {
    my @substrates = split /, /, $tcdb->{substrate};
    foreach my $substrate (@substrates) {
      $tcdbSaved{$tcdb->{curatedIds}} = 1
        if exists $compoundList{lc($substrate)};
    }
  }

  # Search the description for words transport, transporter, exporter, porter, permease, import
  my @out = ();
  my $chits = $dbhC->selectall_arrayref(
    qq{SELECT * FROM CuratedInfo WHERE
            descs LIKE "%transport%"
            OR descs LIKE "%porter%"
            OR descs LIKE "%import%"
            OR descs LIKE "%permease%"
            OR descs like "%PTS system%";},
    { Slice => {} });
  foreach my $chit (@$chits) {
    my @curatedIds = split /,/, $chit->{curatedIds};
    my @descs = split /;; /, $chit->{descs};
    my $keep = 0;
    foreach my $i (0..(scalar(@curatedIds) - 1)) {
      my $curatedId = $curatedIds[$i];
      my $db = $curatedId; $db =~ s/:.*//;
      # TCDB and MetaCyc were handled above
      next unless $db eq "SwissProt" || $db eq "CharProtDB" || $db eq "BRENDA"
        || $db eq "reanno" || $db eq "ecocyc";
      my $desc = $descs[$i];
      next unless $desc =~ m/transport/i
        || $desc =~ m/porter/i
          || $desc =~ m/import/i
            || $desc =~ m/permease/i
              || $desc =~ m/PTS system/i;
      foreach my $cmp (@compoundList) {
        my $pattern = quotemeta($cmp);
        # doing a word-based check this is tricky because substrates may be separated by "/", without spaces
        # or may appear as "glucose-binding"
        # or because another term like "6-phosphate" could be present as the next word.
        # That last issue is not handled.
        # Allow : because of cases like gluconate:H+ symporter
        my $b = "[,/ :]"; # pattern for word boundary
        if ($desc =~ m!^$pattern$b!i # at beginning
            || $desc =~ m!$b$pattern$b!i # in middle
            || $desc =~ m!$b$pattern$!i # at end
            || $desc =~ m!^$pattern-(binding|specific)!i # at beginning
            || $desc =~ m!$b$pattern-(binding|specific)!i) { # in middle
          $keep = 1;
        }
      }
    }
    push @out, $chit if $keep;
  }

  # Add items not already found via curated descriptions
  my %curatedIds = map { $_->{curatedIds} => 1 } @out;
  my @add = keys %metacycSaved;
  push @add, keys %tcdbSaved;
  foreach my $curatedIds (sort @add) {
    if (!exists $curatedIds{$curatedIds}) {
      my $row = $dbhC->selectrow_hashref("SELECT * from CuratedInfo WHERE curatedIds = ?",
                                         { Slice => {} }, $curatedIds);
      push @out, $row if defined $row;
      $curatedIds{$curatedIds} = 1;
    }
  }
  return \@out;
}

my $nSeqShow = 0;
sub TSVPrint($$$) {
  my ($cluster, $ids, $info) = @_;
  $nSeqShow++;
  my @ids = split /,/, $ids;
  die unless @ids > 0;
  my @descs = split /;; /, $info->{descs};
  my @orgs = split /;; /, $info->{orgs} if defined $info->{orgs};
  my @id2s = split /;; /, $info->{id2s} if defined $info->{id2s};
  # fields Sequence Cluster id id2 desc organism

  my $het = GetHetero($dbhC, $ids);
  $het = "heteromer" if defined $het && $het eq "";
  foreach my $i (0..(scalar(@ids)-1)) {
    print join("\t", $nSeqShow, $cluster,
               $ids[$i], $id2s[$i] || "", $descs[$i], $orgs[$i] || "", $het || "")."\n";
  }
}

sub FormatClose($$$$) {
  my ($dbhC, $close, $seqs, $color) = @_;
  my $query = $close->{query} || die;
  my $queryInfo = $close->{queryInfo} || die;
  my $subject = $close->{subject} || die;
  my $subjectInfo = $dbhC->selectrow_hashref("SELECT * FROM CuratedInfo WHERE curatedIds = ?",
                                             {}, $subject);
  die "Unknown close hit $subject" unless defined $subjectInfo;
  my ($subjectSeq) = $dbhC->selectrow_array("SELECT seq FROM CuratedSeq WHERE curatedIds = ?",
                                           {}, $subject);
  die "Unknown close hit $subject" unless defined $subjectSeq;

  my $idString = int($close->{identity})."%";
  my $covString =   "Amino acids $close->{sbeg}:$close->{send}/$subjectInfo->{seqLength}"
    . " are $idString identical to $close->{qbeg}:$close->{qend}/$queryInfo->{seqLength}"
    . " of the $step protein";

  my $subjectShort = $subject; $subjectShort =~ s/,.*//;
  return p({ -style => "margin-bottom: 0em; color: $color" },
          CuratedToHtml($subjectInfo, $subjectSeq)
          . ", " . small($subjectShort),
          br(),
          a({-title => $covString}, "$idString identical to")),
            p({-style => "margin-left: 5em; margin-top: 0em; font-size: 90%; color: $color"},
              CuratedToHtml($queryInfo, $seqs->{$query})),
           "\n";
}

sub HeteroToText($) {
  my ($het) = @_;
  return "" unless defined $het;
  $het = "heteromeric" if $het eq "";
  return " " . $het;
}

#!/usr/bin/perl -w
# Identify the best candidates for each step, score them, and score the rules
use strict;
use Getopt::Long;
use List::Util qw{min max};
use FindBin qw{$RealBin};
use DBI;
use lib "$RealBin/../lib";
use pbutils qw{ReadTable};
use Steps qw{ReadOrgTable ReadOrgProtein ParseOrgLocus};

my $maxCand = 5;
my @weightsDef = (-2,-0.1,1); # default weights for rules with low/medium/high candidates
my $fOverlap = 0.5;
my $usage = <<END
Usage: gapsummary.pl -set set -orgs orgprefix
  -hits hitsfile -revhits revhitsfile -out summary

The hits and revhits files are from gapsearch.pl and gaprevsearch.pl
Writes 1 line per organism x step to summary.steps
Writes 1 line per step x gene candidate to summary.cand
Writes 1 line per rule to summary.rules

Currently, the score of a candidate for a step is defined as

2: blast to a characterized protein at above 40% identity and 80%
   coverage and bits >= otherBits+10,
   or, hmm match and 80% coverage and not(otherIdentity >= 40 &
   otherCoverage >= 0.75)
1: blast to a characterized or curated protein at above 30% id. and
   80% cov. and bits >= otherBits, or blast above 40% id. and 70%
   coverage (ignoring otherBits)
   or hmm match (regardless of coverage or otherBits)
0: other blast hit with 50% coverage

This script also searches for split loci by trying to join two BLAST
hits to the same characterized or curated protein together (unless
there are already two candidates with score=2). It considers loci for
joining if they have at least 30% identity and bits >= otherBits+10
but do not meet the coverage requirements (score < 2). Split loci are
indicated by a non-empty locusId2 in the summary.cand file and by
comma-delimited locusIds or sysNames in the summary.steps file.

The "best" path for a rule is the one that has the no low-confidence
or no medium-confidence steps. If there are no such paths, or there is
a tie, it chooses the path with the highest total score (see -weights
below). If there is still a tie, it chooses the longer path.

Optional arguments:
-pathway pathwayId -- score only the queries for this pathway
-dbDir dbDir -- which directory the input databases are in
  (defaults to $RealBin/../tmp/path.set/)
-maxCand $maxCand -- number of candidates for each step to keep
-weights @weightsDef - the weight for
  each type of step.
-overlap $fOverlap -- ignore hits to other sequences that do not
   overlap at least this fraction of the original hit's alignment.
-noSplit -- do not search for split loci
END
;

# Given a list of hits and a list of reverse hits, all relating to the same orgId/locusId,
# return a hash that describes the best hit, with the fields
#   orgId, locusId, score,
#   (blast-based) blastBits, curatedIds, identity, blastCoverage, blastScore
#   (hmm-based) hmmBits, hmmId, hmmCoverage, hmmScore
#   (other-based) otherBits, otherIdentity, otherCoverage
# or else returns undef if there is no useful hit.
# Assumes that reverse hits have highest bit score first
# (but does not assume that the hits are sorted)
sub ScoreCandidate($$);

# locusBegin, locusEnd, and list of reverse hits => first one meeting overlap requirement
# If none are relevant, returns undef
sub RelevantRevhit($$$);

# Given a hash of locusId => hits, and a hash of locusId => revhits,
# try to find a split blast hit amongst the main hits.
# (Ignores all HMM hits in the input)
# The third argument is the existing candidate scores -- ignores merges
#   unless the score improves over both components.
# Returns a candidate object (but with no HMM fields filled out), or undef
sub FindSplit($$$);

# Given two hits of different loci to the same curated item,
# along with the best score of each locus and reverse
# hit information, returns a combined hit or undef.
# The order of arguments is:
# hit1, bestscore1, reverse hit list for hit1,
# hit2, bestscore2, most relevant reverse hit for hit2
sub MergeHits($$$$$$);

{

  my ($hitsFile, $revhitsFile, $orgprefix, $outpre);
  my ($queryDir, $set, $pathSpec, $dbDir, $noSplit);
  my @weights;
  die $usage
    unless GetOptions('set=s' => \$set,
                      'orgs=s' => \$orgprefix,
                      'hits=s' => \$hitsFile,
                      'revhits=s' => \$revhitsFile,
                      'out=s' => \$outpre,
                      'pathway=s' => \$pathSpec,
                      'dbDir=s' => \$dbDir,
                      'maxCand=i' => \$maxCand,
                      'weights=f{3,3}' => \@weights,
                      'overlap=f' => \$fOverlap,
                      'noSplit' => \$noSplit)
      && defined $set && defined $orgprefix
      && defined $hitsFile && defined $revhitsFile && defined $outpre;
  @weights = @weightsDef unless @weights;
  die "Must have 3 weights\n" unless @weights == 3;
  die "Weight(low) is above Weight(medium)\n" if $weights[0] > $weights[1];
  die "Weight(medium) is above Weight(high)\n" if $weights[1] > $weights[2];
  $dbDir = "$RealBin/../tmp/path.$set" unless defined $dbDir;
  die "No such directory: $dbDir\n" unless -d $dbDir;
  foreach my $file ($hitsFile, $revhitsFile) {
    die "No such file: $file\n" unless -s $file;
  }
  die "-maxCand must be at least 1\n" unless $maxCand >= 1;

  # Load all the inputs
  my @orgs = ReadOrgTable("$orgprefix.org");
  die "The organism table $orgprefix.org has no rows\n" unless @orgs > 0;
  my %orgs = map { $_->{orgId} => $_ } @orgs;
  my $aaIn = "$orgprefix.faa";
  die "No such file: $aaIn\n" unless -e $aaIn;
  my %locusInfo = (); # orgId => locusId => [length,sysName,desc]
  open(my $fhIn, "<", $aaIn) || die "Cannot read $aaIn\n";
  my $faastate = {};
  while (my $prot = ReadOrgProtein($fhIn,$faastate)) {
    my $orgId = $prot->{orgId};
    die "Unknown organism $orgId in $aaIn\n" unless exists $orgs{$orgId};
    $locusInfo{$orgId}{ $prot->{locusId} } = [ length($prot->{aaseq}), $prot->{sysName}, $prot->{desc} ];
  }
  close($fhIn) || die "Error reading $aaIn\n";

  my $dbhC = DBI->connect("dbi:SQLite:dbname=${dbDir}/curated.db","","",{ RaiseError => 1 }) || die $DBI::errstr;
  my $dbhS = DBI->connect("dbi:SQLite:dbname=${dbDir}/steps.db","","",{ RaiseError => 1 }) || die $DBI::errstr;

  my %pathwayDesc = (); # pathwayId to desc
  foreach my $row (@{ $dbhS->selectall_arrayref("SELECT * from Pathway", { Slice => {} }) }) {
    $pathwayDesc{$row->{pathwayId}} = $row->{desc};
  }
  if (defined $pathSpec) {
    die "Unknown pathway $pathSpec" unless exists $pathwayDesc{$pathSpec};
  }
  my %pathwaySteps = (); # pathwayId to list of stepids
  foreach my $pathwayId (keys %pathwayDesc) {
    $pathwaySteps{$pathwayId} = $dbhS->selectcol_arrayref("SELECT stepId FROM Step WHERE pathwayId = ?",
                                                          {}, $pathwayId);
  }

  # Load queries
  my %queryDesc = (); # blast queryId => desc
  my %queryLength = (); # blast queryId => protein length
  # a queryId is as in the output of gapsearch.pl -- curatedIds, uniprot:uniprotId, curated2:protId
  my %blastQuery = (); # pathwayId => stepId => queryId => 1
  my %blastToStep = (); # blast queryId => pathwayId => stepId => 1
  my %hmmQuery = (); # pathwayId => stepId => hmmId => 1
  my %hmmToStep = (); # hmmId => pathwayId => stepId => 1
  my %hmmDesc = ();
  # This includes both proteins assigned to the step (which implies that an other hit is ok)
  # and explicitly ignored items
  my %ignore = (); # pathwayId => stepId => queryId => 1

  foreach my $row (@{ $dbhS->selectall_arrayref("SELECT * from StepQuery", { Slice => {} }) }) {
    next if defined $pathSpec && $row->{pathwayId} ne $pathSpec;
    my $queryType = $row->{queryType} || die;
    my $pathwayId = $row->{pathwayId} || die;
    my $stepId = $row->{stepId} || die;

    if ($queryType eq "ignore") {
      my $id = $row->{curatedIds} || die;
      $ignore{$pathwayId}{$stepId}{$id} = 1;
    } elsif ($queryType eq "hmm") {
      my $hmmId = $row->{hmmId} || die;
      $hmmQuery{$pathwayId}{$stepId}{$hmmId} = 1;
      $hmmToStep{$hmmId}{$pathwayId}{$stepId} = 1;
      $hmmDesc{$hmmId} = $row->{desc};
    } else {
      my $id;
      if ($row->{queryType} eq "curated") {
        $id = $row->{curatedIds} || die;
      } elsif ($row->{queryType} eq "uniprot") {
        $id = "uniprot:" . ($row->{uniprotId} || die);
      } elsif ($row->{queryType} eq "curated2") {
        $id = "curated2:" . ($row->{protId} || die);
      } else {
        die "Unknown query type $row->{queryType}";
      }
      $blastQuery{$pathwayId}{$stepId}{$id} = 1;
      $blastToStep{$id}{$pathwayId}{$stepId} = 1;
      $ignore{$pathwayId}{$stepId}{$id} = 1;
      # And save the description
      if (exists $queryDesc{$id}) {
        die "Mismatching descriptions for $id"
          unless $queryDesc{$id} eq $row->{desc};
      } else {
        $queryDesc{$id} = $row->{desc} || die;
        $queryLength{$id} = length($row->{seq})
          || die "length of seq for $id is 0";
      }
    }
  }

  # Both hits and revhits initially have locusId of the form orgId:locusId
  # This is converted to an orgId and a locusId below
  my @hitFields = qw{locusId type queryId bits locusBegin locusEnd qBegin qEnd qLength identity};
  my @hits = ReadTable($hitsFile, \@hitFields);
  # add the fromCurated2 field to all blast hits -- these are lower priority because
  # the reference protein is not actually characterized
  foreach my $hit (@hits) {
    if ($hit->{type} eq "blast") {
      $hit->{fromCurated2} = $hit->{queryId} =~ m/^curated2:/ ? 1 : 0;
    }
  }
  my @revFields = qw{locusId otherId bits locusBegin locusEnd otherBegin otherEnd otherIdentity};
  my @revhits = ReadTable($revhitsFile, \@revFields);

  my %hits; # orgId => pathwayId => stepId => locusId => list of relevant hits
  foreach my $hit (@hits) {
    my $queryId = $hit->{queryId};
    my ($orgId, $locusId) = ParseOrgLocus($hit->{locusId});
    $hit->{orgId} = $orgId;
    $hit->{locusId} = $locusId;
    my $hash; # pathway => step => 1
    if ($hit->{type} eq "hmm") { # queryId is an hmm id
      # If looking at one pathway, not all queries are loaded,
      # so do not check for unexpected queries
      die "Unknown hmm query $queryId found in hits file $hitsFile\n"
        unless defined $pathSpec || exists $hmmToStep{$queryId};
      $hash = $hmmToStep{$queryId};
    } else {
      die "Unknown hit type $hit->{type} in hits file $hitsFile\n"
        unless defined $pathSpec || $hit->{type} eq "blast";
      die "Unknown blast query $queryId found in hits file $hitsFile\n"
        unless defined $pathSpec || exists $blastToStep{$queryId};
      $hash = $blastToStep{$queryId};
    }
    while (my ($pathwayId, $stephash) = each %$hash) {
      foreach my $stepId (keys %$stephash) {
        push @{ $hits{$orgId}{$pathwayId}{$stepId}{$locusId} }, $hit;
      }
    }
  }

  my %revhits = (); # orgId => locusId => list of relevant reverse hits
  foreach my $revhit (@revhits) {
    my ($orgId, $locusId) = ParseOrgLocus($revhit->{locusId});
    $revhit->{orgId} = $orgId;
    $revhit->{locusid} = $locusId;
    push @{ $revhits{$orgId}{$locusId} }, $revhit;
  }

  # add locusLength to each hit (it already has qLength)
  # also add coverage (of the curated item)
  foreach my $hit (@hits) {
    my $orgId = $hit->{orgId};
    my $locusId = $hit->{locusId};
    my $queryId = $hit->{queryId} || die;
    die "Unknown locus $orgId $locusId in $hitsFile but not in $aaIn\n"
      unless exists $locusInfo{$orgId}{$locusId};
    my $len = $locusInfo{$orgId}{$locusId}[0];
    $hit->{locusLength} = $len;
    $hit->{coverage} = ($hit->{qEnd} - $hit->{qBegin} + 1) / $hit->{qLength};
    unless (exists $hmmToStep{$queryId}) {
      # this is a blast query, should already know its length
      die "Unexpected query $queryId of unknown length"
        unless defined $pathSpec || exists $queryLength{$queryId};
      die "Mismatched lengths for $queryId"
        unless defined $pathSpec || $hit->{qLength} == $queryLength{$queryId};
    }
  }

  # Fetch all sequence lengths, for computing other coverage
  my %curatedLength = (); # curatedIds => sequence length
  my $curatedLengths = $dbhC->selectall_arrayref("SELECT curatedIds, seqLength FROM CuratedInfo");
  foreach my $row (@$curatedLengths) {
    $curatedLength{ $row->[0] } = $row->[1];
  }

  # add otherLength and otherCoverage to each revhit
  # (Note -- most of this work is unnecessary -- otherCoverage
  #  is only important for HMM hits. This can be fixed later.)
  foreach my $revhit (@revhits) {
    my $otherId = $revhit->{otherId};
    if (!exists $queryLength{$otherId}) {
      my $length = $curatedLength{$otherId};
      die "Unknown curated id $otherId in rev hits $revhitsFile\n"
        unless $length;
      $queryLength{$otherId} = $length;
    }
    $revhit->{otherLength} = $queryLength{$otherId};
    die "Illegal reverse hit length for $otherId: $revhit->{otherEnd} vs $queryLength{$otherId}"
      if $revhit->{otherEnd} > $queryLength{$otherId};
    $revhit->{otherCoverage} = ($revhit->{otherEnd} - $revhit->{otherBegin} + 1) / $queryLength{$otherId};
  }

  # Score each candidate based on its various hits, the reverse hits,
  # and the ignore information
  my %cand = (); # orgId => pathwayId => stepId => sorted list of candidates
  while (my ($orgId, $pathhash) = each %hits) {
    while (my ($pathwayId, $stephash) = each %$pathhash) {
      next if defined $pathSpec && $pathwayId ne $pathSpec;
      my %locusRev = ();
      while (my ($step, $locushash) = each %$stephash) {
        # score the candidates for this step
        my @cand = ();
        while (my ($locusId, $hits) = each %$locushash) {
          my @rev = ();
          @rev = @{ $revhits{$orgId}{$locusId} } if exists $revhits{$orgId}{$locusId};
          @rev = grep !exists $ignore{$pathwayId}{$step}{ $_->{otherId} }, @rev;
          @rev = sort { $b->{bits} <=> $a->{bits}
                        || $a->{otherId} cmp $b->{otherId} } @rev;
          $locusRev{$locusId} = \@rev; # saved for use in splitting below
          my $cand = ScoreCandidate($hits, \@rev);
          push @cand, $cand if defined $cand;
        }
        # select the best candidates for this step
        # push curated-not-characterized items after other items with the same score
        foreach my $cand (@cand) {
          $cand->{maxBits} = max($cand->{blastBits} || 0, $cand->{hmmBits} || 0);
        }
        @cand = sort { $b->{score} <=> $a->{score}
                         || $a->{fromCurated2} <=> $b->{fromCurated2}
                         || $b->{maxBits} <=> $a->{maxBits}
                         || $b->{locusId} cmp $a->{locusId} } @cand;
        @cand = splice(@cand, 0, $maxCand) if @cand > $maxCand;

        # look for split ORFs unless there are two high-confidence candidates
        unless (defined $noSplit || (@cand >= 2 && $cand[1]{score} == 2)) {
          my %locusCand = map { $_->{locusId} => $_ } @cand;
          my $merge = FindSplit($locushash, \%locusRev, \%locusCand);
          if (defined $merge) {
            # Put the HMM and fromCurated2 info from 1st (higher scoring) locus into $merge
            my $cand1 = $locusCand{ $merge->{locusId} };
            die unless defined $cand1;
            foreach my $key (qw{hmmBits hmmId hmmCoverage hmmScore fromCurated2}) {
              $merge->{$key} = $cand1->{$key};
            }
            $merge->{maxBits} = max($merge->{blastBits}, $merge->{hmmBits} || 0);
            # replace the previous entries for the merged genes with the merge
            @cand = grep { $_->{locusId} ne $merge->{locusId}
                             && $_->{locusId} ne $merge->{locusId2} } @cand;
            push @cand, $merge;
            # Re-sort and re-truncate the list
            @cand = sort { ($b->{score} || 0) <=> ($a->{score} || 0)
                             || ($a->{fromCurated2} || 0) <=> ($b->{fromCurated2}||0)
                             || $b->{maxBits} <=> $a->{maxBits}
                             || $b->{locusId} cmp $a->{locusId} } @cand;
            @cand = splice(@cand, 0, $maxCand) if @cand > $maxCand;
          }
        }
        # Make sure all metadata fields are present
        foreach my $cand (@cand) {
          $cand->{orgId} = $orgId;
          $cand->{pathway} = $pathwayId;
          $cand->{step} = $step;
        }
        $cand{$orgId}{$pathwayId}{$step} = \@cand;
      } # end loop over steps
    }
  }

  # Load the rule instances and instance components
  my $instances = $dbhS->selectall_arrayref("SELECT * from RuleInstance ORDER BY instanceId",
                                            { Slice => {} });
  my $components = $dbhS->selectall_arrayref("SELECT * from InstanceComponent ORDER BY instanceId,componentId",
                                            { Slice => {} });
  my %ruleOrder = (); # pathwayId to rules, in order
  my %ruleSeen = (); # pathwayId => ruleId => 1

  my %instanceComponents = (); # instance to list of component objects
  foreach my $component (@$components) {
    my $pathwayId = $component->{pathwayId};
    my $ruleId = $component->{ruleId};
    my $instanceId = $component->{instanceId};
    push @{ $instanceComponents{$instanceId} }, $component;
    unless (exists $ruleSeen{$pathwayId}{$ruleId}) {
      $ruleSeen{$pathwayId}{$ruleId} = 1;
      push @{ $ruleOrder{$pathwayId} }, $ruleId;
    }
  }
  my %rulePaths = (); # pathway to rule to list of paths, each of which is a
  # list of component objects
  foreach my $instance (@$instances) {
    my $pathwayId = $instance->{pathwayId};
    my $ruleId = $instance->{ruleId};
    my $instanceId = $instance->{instanceId};
    my $components = $instanceComponents{$instanceId}
      || die "No components of instance $instance for $pathwayId $ruleId";
    push @{ $rulePaths{$pathwayId}{$ruleId} }, $components;
  }

  # And score the rules
  # To do this, need to put the rules in dependency order.
  # I.e., first the terminals (no instance has any dependencies),
  # then the ones with a "height" of 1 (no instance has any non-terminal children), etc.
  my %dependencyOrder = (); # pathwayId to rules, in order to score (most-specific first, all last)
  foreach my $pathwayId (sort keys %pathwayDesc) {
    next if $pathwayId eq "all";
    next if defined $pathSpec && $pathSpec ne $pathwayId;
    my $rulePaths = $rulePaths{$pathwayId};

    # In each round, compute heights 1 level farther up
    my %ruleToHeight = ();
    my $maxRounds = 10000;
    my $nRound;
    for($nRound = 0; $nRound < $maxRounds; $nRound++) {
      my $nHeightsOld = scalar(keys %ruleToHeight);
      foreach my $ruleId (sort keys $rulePaths) {
        my $maxHeight = 0;
        my $paths = $rulePaths->{$ruleId} || die "Unknown $pathwayId $ruleId";
        foreach my $path (@$paths) {
          foreach my $component (@$path) {
            if ($component->{subRuleId} ne "") {
              if (exists $ruleToHeight{ $component->{subRuleId} }) {
                $maxHeight = max($maxHeight, 1 + $ruleToHeight{ $component->{subRuleId} });
              } else {
                $maxHeight = undef;
                last;
              }
            }
          }
          last if !defined $maxHeight;
        }
        $ruleToHeight{$ruleId} = $maxHeight if defined $maxHeight;
      }
      my $nHeights = scalar(keys %ruleToHeight);
      die $nRound if $nHeights == $nHeightsOld;
      last if exists $ruleToHeight{"all"};
    }
    die "Too many rounds, cyclic dependencies?" if $nRound == $maxRounds;
    die "Wrong #heights" unless scalar(keys %ruleToHeight) == scalar(keys %$rulePaths);
    my @rulesInOrder = sort { $ruleToHeight{$a} <=> $ruleToHeight{$b}
                                || $a cmp $b } keys %ruleToHeight;
    $dependencyOrder{$pathwayId} = \@rulesInOrder;
  }

  my %ruleScores = (); # orgId => pathway => rulename => hash of n, n012, score, path, path2, stepsUsed, pathExpanded
  # where n is #steps, n012 is a vector with #of steps at each score,
  # score is the total weighted score,
  # path and path2 are lists of paths (lists of component objects)
  #	path2 is used only if there is a tie
  # stepsUsed is a hash that includes all steps used, including steps for dependencies
  # pathExpanded is the list of all stepIds including those in dependencies (with duplicates removed)
  # Note that n and score need to be recomputed each time because of potential overlap
  foreach my $orgId (sort keys %orgs) {
    foreach my $pathwayId (sort keys %pathwayDesc) {
      next if defined $pathSpec && $pathSpec ne $pathwayId;
      foreach my $ruleId (@{ $dependencyOrder{$pathwayId} }) {
        my @scoredPaths = (); # list of hashes with n, score, n012, stepsUsed, and path
        die "No instances for rule $ruleId in pathway $pathwayId"
          unless exists $rulePaths{$pathwayId}{$ruleId};
        foreach my $path (@{ $rulePaths{$pathwayId}{$ruleId} }) {
          my %stepsUsed = (); # steps used in this best approach for this path
          # (organism dependent due to scoring of dependencies)
          foreach my $component (@$path) {
            my $subRuleId = $component->{subRuleId};
            if ($component->{stepId} ne "") {
              $stepsUsed{ $component->{stepId} } = 1;
            } elsif (exists $ruleScores{$orgId}{$pathwayId}{$subRuleId}) {
              foreach my $stepId (keys %{ $ruleScores{$orgId}{$pathwayId}{$subRuleId}{stepsUsed} }) {
                $stepsUsed{$stepId} = 1;
              }
            } else {
              die "No score yet for dependency $subRuleId while scoring $ruleId for pathway $pathwayId\n";
            }
          }
          # score this path
          my $n = 0;
          my @n012 = (0,0,0); # number of steps with score of 0, 1, or 2
          my $scoreTot = 0;
          foreach my $stepId (keys %stepsUsed) {
            my $cand = $cand{$orgId}{$pathwayId}{$stepId}[0]
              if exists $cand{$orgId}{$pathwayId}{$stepId} && @{ $cand{$orgId}{$pathwayId}{$stepId} } > 0;
            my $scoreThis = defined $cand ? $cand->{score} : 0;
            $n++;
            die unless $scoreThis >= 0 && $scoreThis <= 2;
            $n012[$scoreThis]++;
            $scoreTot += $weights[$scoreThis];
          }

          my $minScore = 0;
          $minScore = 1 if $n012[0] == 0; # all are medium or above
          $minScore = 2 if $n012[0] == 0 && $n012[1] == 0;

          push @scoredPaths, { 'n' => $n,
                               'score' => $scoreTot, 'minScore' => $minScore,
                               'n012' => \@n012,
                               'stepsUsed' => \%stepsUsed,
                               'path' => $path };
        }
        # select the best 1 or 2 scored paths, based on highest minScore,
        # or highest weighted (total) score, or highest #steps
        @scoredPaths = sort { $b->{minScore} <=> $a->{minScore}
                              || $b->{score} <=> $a->{score}
                              || $b->{n} <=> $a->{n} } @scoredPaths;
        die unless @scoredPaths > 0;
        $ruleScores{$orgId}{$pathwayId}{$ruleId} = $scoredPaths[0];
        $ruleScores{$orgId}{$pathwayId}{$ruleId}{path2} = $scoredPaths[1]{path} if @scoredPaths >= 2;
        # and compute pathExpanded for the main path
        my $path = $ruleScores{$orgId}{$pathwayId}{$ruleId}{path};
        my @pathExpanded = ();
        foreach my $component (@$path) {
          if ($component->{stepId} ne "") {
            push @pathExpanded, $component->{stepId};
          } else {
            my $pathPart = $ruleScores{$orgId}{$pathwayId}{ $component->{subRuleId} }{pathExpanded};
            die unless defined $pathPart;
            push @pathExpanded, @$pathPart;
          }
        }
        my %pathSoFar = ();
        @pathExpanded = grep { my $keep = !exists $pathSoFar{$_};
                               $pathSoFar{$_} = 1;
                               $keep; } @pathExpanded;
        $ruleScores{$orgId}{$pathwayId}{$ruleId}{pathExpanded} = \@pathExpanded;
      }
    }
  }

  # Record which steps are on the best path for the 'all' rule
  my %onBestPath = (); # orgId => pathId => step => 1
  foreach my $orgId (sort keys %orgs) {
    foreach my $pathwayId (sort keys %pathwayDesc) {
      next if $pathwayId eq "all";
      next if defined $pathSpec && $pathwayId ne $pathSpec;
      my $path = $ruleScores{$orgId}{$pathwayId}{'all'}{pathExpanded};
      die "No path for all for $orgId $pathwayId" unless defined $path && @$path > 0;
      foreach my $stepId (@$path) {
        $onBestPath{$orgId}{$pathwayId}{$stepId} = 1;
      }
    }
  }

  # Write out the candidates and steps
  open(my $fhCand, ">", "$outpre.cand") || die "Cannot write to $outpre.cand\n";
  open(my $fhStep, ">", "$outpre.steps") || die "Cannot write to $outpre.steps";
  my @candfields = qw{orgId gdb gid pathway step score
                      locusId sysName desc
                      locusId2 sysName2 desc2
                      blastBits curatedIds identity blastCoverage blastScore curatedDesc
                      hmmBits hmmId hmmCoverage hmmScore hmmDesc
                      otherIds otherBits otherIdentity otherCoverage};
  print $fhCand join("\t", @candfields)."\n";
  my @stepfields = qw{orgId gdb gid pathway step onBestPath score locusId sysName score2 locusId2 sysName2};
  print $fhStep join("\t", @stepfields)."\n";
  foreach my $orgId (sort keys %cand) {
    foreach my $pathwayId (sort keys %pathwayDesc) {
      next if defined $pathSpec && $pathwayId ne $pathSpec;
      my $steps = $pathwaySteps{$pathwayId};
      my $stephash = $cand{$orgId}{$pathwayId} || {};
      foreach my $step (@$steps) {
        my $cands = $stephash->{$step} || [];
        # per-candidate output
        foreach my $cand (@$cands) {
          my (undef, $sysName, $desc) = @{ $locusInfo{$orgId}{$cand->{locusId}} };
          $cand->{sysName} = $sysName;
          $cand->{desc} = $desc;
          if ($cand->{locusId2}) {
            my (undef, $sysName2, $desc2) = @{ $locusInfo{$orgId}{$cand->{locusId2}} };
            $cand->{sysName2} = $sysName2;
            $cand->{desc2} = $desc2;
          }
          $cand->{curatedDesc} = $queryDesc{ $cand->{curatedIds} } if $cand->{curatedIds};
          $cand->{hmmDesc} = $hmmDesc{ $cand->{hmmId} } if $cand->{hmmId};
          $cand->{gdb} = $orgs{$orgId}{gdb};
          $cand->{gid} = $orgs{$orgId}{gid};
          my @out = map { defined $cand->{$_} ? $cand->{$_} : "" } @candfields;
          print $fhCand join("\t", @out)."\n";
        }
        # per-step line
        my @stepout = ($orgId, $orgs{$orgId}{gdb}, $orgs{$orgId}{gid},
                       $pathwayId, $step, exists $onBestPath{$orgId}{$pathwayId}{$step} ? 1 : 0);
        foreach my $i (0,1) {
          my $c = $cands->[$i] if $i < @$cands;
          if (!defined $c) {
            push @stepout, ("", "", "");
          } else {
            my $locusId = $c->{locusId};
            my $sysName = $locusInfo{$orgId}{$locusId}[1];
            if ($c->{locusId2}) {
              $locusId .= "," . $c->{locusId2};
              $sysName .= "," . $locusInfo{$orgId}{ $c->{locusId2} }[1];
            }
            push @stepout, $c->{score}, $locusId, $sysName;
          }
        }
        print $fhStep join("\t", @stepout)."\n";
      }
    }
  }
  close($fhCand) || die "Error writing to $outpre.cand\n";
  close($fhStep) || die "Error writing to $outpre.steps\n";
  print STDERR "Wrote $outpre.cand and $outpre.steps\n";

  # Write out the best paths
  open (my $fhR, ">", "$outpre.rules") || die "Cannot write to $outpre.rules\n";
  print $fhR join("\t", qw{orgId gdb gid pathway rule nHi nMed nLo score expandedPath path path2})."\n";
  foreach my $orgId (sort keys %orgs) {
    foreach my $pathwayId (sort keys %pathwayDesc) {
      next if $pathwayId eq "all";
      next if defined $pathSpec && $pathwayId ne $pathSpec;
      foreach my $ruleId (@{ $ruleOrder{$pathwayId} }) {
        my $sc = $ruleScores{$orgId}{$pathwayId}{$ruleId};
        die unless exists $sc->{n};
        die unless exists $sc->{path};
        my @pathShow = map { $_->{subRuleId} || $_->{stepId} } @{ $sc->{path} };
        my @path2Show = map { $_->{subRuleId} || $_->{stepId} } @{ $sc->{path2} }
          if defined $sc->{path2};
        print $fhR join("\t", $orgId, $orgs{$orgId}{gdb}, $orgs{$orgId}{gid},
                        $pathwayId, $ruleId,
                        $sc->{n012}[2], $sc->{n012}[1], $sc->{n012}[0], $sc->{score},
                        join(" ", @{ $sc->{pathExpanded} }),
                        join(" ", @pathShow),
                        join(" ", @path2Show))."\n";
      }
    }
  }
  close($fhR) || die "Error writing to $outpre.rules\n";
  print STDERR "Wrote $outpre.rules\n";
}

sub ScoreCandidate($$) {
  my ($hits, $revhits) = @_;
  die "No hits in ScoreCandidate\n" unless @$hits > 0;
  my $orgId = $hits->[0]{orgId};
  my $locusId = $hits->[0]{locusId};

  my @blastHits = sort { $b->{bits} <=> $a->{bits} || $a->{queryId} cmp $b->{queryId} }
    grep { $_->{type} eq "blast" } @$hits;
  my $bestBlastHit;
  foreach my $hit (@blastHits) {
    next unless $hit->{coverage} >= 0.5;
    my $rh = RelevantRevhit($hit->{locusBegin}, $hit->{locusEnd}, $revhits);
    my $otherBits = defined $rh ? $rh->{bits} : -100;
    my $score = 0;
    if ($hit->{identity} >= 40 && $hit->{coverage} >= 0.8 && $hit->{bits} >= $otherBits + 10) {
      $score = 2;
    } elsif ($hit->{identity} >= 30 && $hit->{coverage} >= 0.8 && $hit->{bits} >= $otherBits) {
      $score = 1;
    } elsif ($hit->{identity} >= 40 && $hit->{coverage} >= 0.7) {
      $score = 1;
    }
    $score = 1 if $score == 2 && $hit->{fromCurated2};
    unless (defined $bestBlastHit && $score <= $bestBlastHit->{blastScore}) {
      $bestBlastHit = { 'blastBits' => $hit->{bits},
                        'curatedIds' => $hit->{queryId},
                        'fromCurated2' => $hit->{fromCurated2},
                        'identity' => $hit->{identity},
                        'blastCoverage' => $hit->{coverage},
                        'blastScore' => $score,
                        'otherIds' => $rh ? $rh->{otherId} : "",
                        'otherBits' => $rh ? $rh->{bits} : "",
                        'otherIdentity' => $rh ? $rh->{otherIdentity} : "" ,
                        'otherCoverage' => $rh ? $rh->{otherCoverage} : ""};
      last if $score == 2;
    }
  }

  my @hmmHits = sort { $b->{bits} <=> $a->{bits} } grep { $_->{type} eq "hmm" } @$hits;
  my $bestHMMHit;
  foreach my $hit (@hmmHits) {
    my $rh = RelevantRevhit($hit->{locusBegin}, $hit->{locusEnd}, $revhits);
    my $otherIdentity = $rh ? $rh->{otherIdentity} : 0;
    my $otherCoverage = $rh ? $rh->{otherCoverage} : 0;
    my $score = 1;
    $score = 2 if $hit->{coverage} >= 0.8 && !($otherIdentity >= 40 && $otherCoverage >= 0.75);
    unless (defined $bestHMMHit && $score <= $bestHMMHit->{hmmScore}) {
      $bestHMMHit = { 'hmmBits' => $hit->{bits},
                      'hmmId' => $hit->{queryId},
                      'hmmCoverage' => $hit->{coverage},
                      'hmmScore' => $score,
                      'otherIds' => $rh ? $rh->{otherId} : "",
                      'otherBits' => $rh ? $rh->{bits} : "",
                      'otherIdentity' => $rh ? $rh->{otherIdentity} : "" ,
                      'otherCoverage' => $rh ? $rh->{otherCoverage} : ""};
    }
    last if $score == 2;
  }

  # Possible that none of the BLAST hits meet coverage criteria
  return undef if !defined $bestBlastHit && !defined $bestHMMHit;

  my $out = { 'orgId' => $orgId, 'locusId' => $locusId };
  if ($bestBlastHit) {
    while (my ($key, $value) = each %$bestBlastHit) {
      die "Undefined value for $key" unless defined $value;
      $out->{$key} = $value;
    }
  }
  if ($bestHMMHit) {
    while (my ($key, $value) = each %$bestHMMHit) {
      die "Undefined value for $key" unless defined $value;
      die "Clashing key $key" if exists $out->{$key} && $key !~ m/^other/;
      $out->{$key} = $value;
    }
  }
  $out->{score} = max($bestBlastHit ? $bestBlastHit->{blastScore} : 0,
                      $bestHMMHit ? $bestHMMHit->{hmmScore} : 0);
  $out->{fromCurated2} = 0 unless defined $out->{fromCurated2};
  return $out;
}

sub RelevantRevhit($$$) {
  my ($locusBegin, $locusEnd, $revhits) = @_;
  foreach my $rh (@$revhits) {
    # compute overlap region, if any
    my $begin = max($locusBegin, $rh->{locusBegin});
    my $end = min($locusEnd, $rh->{locusEnd});
    return $rh
      if $end >= $begin
        && $end-$begin+1 >= ($locusEnd-$locusBegin+1) * $fOverlap;
  }
  return undef
}

# Try combining all pairs of blast hits at above 30% and with no
# other-hit that scores 10 higher -- meaning that these could be high-confidence
# hits once merged.
# (Assumes that there are no high-confidence hits)
sub FindSplit($$$) {
  my ($hithash, $revhash, $candhash) = @_;

  my %locusBest = (); # locusId => best hits
  my %locusTarget = (); # locusId => queryId => best hit object
  while (my ($locusId, $hits) = each %$hithash) {
    my @hits = grep { $_->{type} eq "blast" && $_->{identity} >= 30 } @$hits;
    @hits = sort { $a->{fromCurated2} <=> $b->{fromCurated2}
                     || $b->{bits} <=> $a->{bits} } @hits;
    @hits = splice(@hits, 0, $maxCand) if @hits > $maxCand;
    $locusBest{$locusId} = \@hits if @hits > 0;
    foreach my $hit (@hits) {
      $locusTarget{$locusId}{ $hit->{queryId} } = $hit
        unless exists $locusTarget{$locusId}{ $hit->{locusId} };
    }
  }

  my @loci = sort { $locusBest{$b}[0]{bits} <=> $locusBest{$a}[0]{bits} } keys %locusBest;
  @loci = splice(@loci, 0, $maxCand) if @loci > $maxCand;

  # And test all remaining pairs
  my @comb = ();
  for (my $j = 1; $j < scalar(@loci); $j++) {
    my $locus2 = $loci[$j];
    foreach my $hit2 (@{ $locusBest{$locus2} }) {
      # Skip if there is an other hit that is close
      my $rh2 = RelevantRevhit($hit2->{locusBegin}, $hit2->{locusEnd}, $revhash->{$locus2});
      next if defined $rh2 && $hit2->{bits} <= $rh2->{bits} + 10;

      my $queryId = $hit2->{queryId} || die;
      for (my $i = 0; $i < $j; $i++) {
        my $locus1 = $loci[$i]; # the higher-scoring one
        die if $locus1 eq $locus2;
        next unless exists $locusTarget{$locus1}{$queryId};
        my $hit1 = $locusTarget{$locus1}{$queryId};
        die unless $hit1 && $hit2 && $hit1->{queryId} eq $hit2->{queryId};
        my $score1 = exists $candhash->{$locus1}{score} ? $candhash->{$locus1}{score} : -1;
        my $score2 = exists $candhash->{$locus2}{score} ? $candhash->{$locus2}{score} : -1;
        my $combHit = MergeHits($hit1, $score1, $revhash->{$locus1},
                                $hit2, $score2, $rh2);
        push @comb, $combHit if defined $combHit;
      }
    }
  }
  @comb = sort { $b->{score} <=> $a->{score}
                   || $b->{blastBits} <=> $a->{blastBits} } @comb;
  # Return just one hit as finding these is uncommon and we don't know if hits are
  # redundant (i.e., same pair of loci but different queryId)
  return $comb[0]; # may be undef
}

sub MergeHits($$$$$$){
  my ($hit1, $bestscore1, $revhits1,
      $hit2, $bestscore2, $revhit2) = @_;
  die unless $hit1->{queryId} eq $hit2->{queryId};
  my $queryId = $hit1->{queryId};

  my $locus1 = $hit1->{locusId};
  my $locus2 = $hit2->{locusId};
  die unless defined $locus1 && defined $locus2 && $locus1 ne $locus2;

  my $alnLen1 = $hit1->{qEnd} - $hit1->{qBegin}+1;
  my $alnLen2 = $hit2->{qEnd} - $hit2->{qBegin}+1;

  # Overlap on the curated sequence should be at most 20% of either alignment
  # Combined alignment should cover at least 70% of the curated protein
  my ($overlap, $combLen);
  if ($hit1->{qEnd} < $hit2->{qBegin} || $hit1->{qBegin} > $hit2->{qEnd}) {
    # No overlap
    $overlap = 0;
    $combLen = $alnLen1 + $alnLen2;
  } else {
    my $combBegin = min($hit1->{qBegin}, $hit2->{qBegin});
    my $combEnd = max($hit1->{qEnd}, $hit2->{qEnd});
    $combLen = $combEnd - $combBegin + 1;
    $overlap = $alnLen1 + $alnLen2 - $combLen;
  }
  return undef unless $overlap < 0.2 * $alnLen1
    && $overlap < 0.2 * $alnLen2
    && $combLen >= 0.7 * $hit1->{qLength};


  # Check that there is no better other-hit for locus1 and locus2
  return undef if defined $revhit2 && $hit2->{bits} <= $revhit2->{bits} + 10;

  my $revhit1 = RelevantRevhit($hit1->{locusBegin}, $hit1->{locusEnd}, $revhits1);
  return undef if defined $revhit1 && $hit1->{bits} <= $revhit1->{bits} + 10;

  # Score the hit and check that it is better than either component
  my $combIdentity = ($alnLen1 * $hit1->{identity} + $alnLen2 * $hit2->{identity}) / ($alnLen1 + $alnLen2);
  my $combCoverage = $combLen / $hit1->{qLength};
  my $combScore = $combIdentity >= 40 && $combCoverage >= 0.8 ? 2 : 1;
  $combScore = 1 if $combScore == 2 && $hit1->{fromCurated2};

  return undef unless $combScore > $bestscore1 && $combScore > $bestscore2;

  # Ensure that the hash references below succeed
  $revhit1 = {} if !defined $revhit1;
  $revhit2 = {} if !defined $revhit2;

  # Return a merged hit
  return { 'locusId' => $locus1,
           'locusId2' => $locus2,
           'curatedIds' => $queryId,
           'blastBits' => $hit1->{bits} + $hit2->{bits},
           'identity' => $combIdentity,
           'blastCoverage' => $combCoverage,
           'blastScore' => $combScore,
           'score' => $combScore,
           'otherIds' => $revhit1->{otherId} || $revhit2->{otherId} || "",
           'otherBits' => ($revhit1->{bits} || 0) + ($revhit2->{bits} || 0),
           'otherIdentity' => max($revhit1->{otherIdentity} || 0, $revhit2->{otherIdentity} || 0),
           'otherCoverage' => max($revhit1->{otherCoverage} || 0, $revhit2->{otherCoverage} || 0)
         };
}

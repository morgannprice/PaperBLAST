#!/usr/bin/perl -w
# Identify the best candidates for each step, score them, and score the rules
use strict;
use Getopt::Long;
use List::Util qw{min max};
use FindBin qw{$Bin};
use lib "$Bin/../lib";
use pbutils qw{ReadTable};
use Steps qw{ReadSteps ReadOrgTable ReadOrgProtein ParseOrgLocus};

my $maxCand = 5;
my @weightsDef = (-2,-0.1,1); # default weights for rules with low/medium/high candidates
my $infoFile = "curated.faa.info";
my $fOverlap = 0.5;
my $stepDir = ".";
my $queryDir = ".";
my $usage = <<END
Usage: gapsummary.pl -pathways his met ... pro -orgs orgprefix
  -hits hitsfile -revhits revhitsfile -out summary

The hits and revhites files are from gapsearch.pl and gaprevsearch.pl
Writes 1 line per organism x step to summary.steps
Writes 1 line per step x gene candidate to summary.cand
Writes 1 line per rule to summary.rules

Currently, the score of a candidate for a step is defined as

2: blast to a characterized protein at above 40% identity and 80%
   coverage and bits >= otherBits+10,
   or hmm match and 80% coverage and !(otherIdentity >= 40 &
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
-info $infoFile -- made by curatedFaa.pl with the -curatedids
  option
-stepDir $stepDir -- which directory the *.steps files are in
-queryDir $queryDir -- which directory the *.query files are in
-maxCand $maxCand -- number of candidates for each step to keep
-weights weightLow weightMedium weightHigh -- the weight for
  each type of step. The default is @weightsDef.
-overlap $fOverlap -- ignore hits to other sequences that do not
   overlap at least this fraction of the original hit's alignment.
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
  my @pathways;
  my @weights;
  my ($hitsFile, $revhitsFile, $orgprefix, $outpre);
  die $usage
    unless GetOptions('pathways=s{1,}' => \@pathways,
                      'orgs=s' => \$orgprefix,
                      'hits=s' => \$hitsFile,
                      'revhits=s' => \$revhitsFile,
                      'out=s' => \$outpre,
                      'info=s' => \$infoFile,
                      'stepDir=s' => \$stepDir,
                      'queryDir=s' => \$queryDir,
                      'maxCand=i' => \$maxCand,
                      'weights=f{3,3}' => \@weights,
                      'overlap=f' => \$fOverlap)
      && @pathways > 0 && defined $orgprefix
      && defined $hitsFile && defined $revhitsFile && defined $outpre;
  @weights = @weightsDef unless @weights;
  die "Must have 3 weights\n" unless @weights == 3;
  die "Weight(low) is above Weight(medium)\n" if $weights[0] > $weights[1];
  die "Weight(medium) is above Weight(high)\n" if $weights[1] > $weights[2];
  foreach my $dir ($stepDir, $queryDir) {
    die "No such directory: $dir\n" unless -d $dir;
  }
  foreach my $file ($hitsFile, $revhitsFile, $infoFile) {
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

  my %pathways = (); # steps and rules items for each pathway
  my %queryDesc = (); # id => desc
  my %ignore = (); # pathwayId => step => id => 1
  my %blastQuery = (); # pathwayId => step => id => 1
  my %hmmQuery = (); # pathwayId => step => hmm
  my %hmmToStep = (); # hmm => pathwayId => step => 1
  my %queryToStep = (); # curated id => pathwayId => step => 1
  foreach my $pathId (@pathways) {
    die "Duplicate pathway $pathId\n" if exists $pathways{$pathId};
    my $stepFile = "$stepDir/$pathId.steps";
    my $queryFile = "$queryDir/$pathId.query";
    foreach my $file ($stepFile, $queryFile) {
      die "No such file: $file\n" unless -e $file;
    }
    $pathways{$pathId} = ReadSteps($stepFile);
    my @header = qw{step type query desc file sequence};
    my @queries = ReadTable($queryFile, \@header);
    foreach my $query (@queries) {
      my $id = $query->{query};
      my $type = $query->{type};
      my $step = $query->{step};
      $id = "uniprot:" . $id if $type eq "uniprot";
      $id = "curated2:" . $id if $type eq "curated2";
      if ($type eq "curated" || $type eq "curated2" || $type eq "uniprot") {
        $blastQuery{$pathId}{$step}{$id} = 1;
        $queryToStep{$id}{$pathId}{$step} = 1;
        $ignore{$pathId}{$step}{$id} = 1;
      } elsif ($type eq "hmm") {
        $hmmQuery{$pathId}{ $step }{ $id } = 1;
        $hmmToStep{$id}{$pathId}{$step} = 1;
      } elsif ($type eq "ignore") {
        $ignore{$pathId}{$step}{$id} = 1;
      } else {
        die "Unknown query type $type in $queryFile\n";
      }
      # And save descriptions
      if ($type eq "curated" || $type eq "curated2" || $type eq "uniprot" || $type eq "hmm") {
        if (exists $queryDesc{$id}) {
          die "Mismatch for description of $id in $queryFile -- $queryDesc{$id} vs. $query->{desc}\n"
            unless $queryDesc{$id} eq $query->{desc};
        } else {
          $queryDesc{$id} = $query->{desc};
        }
      }
    }
  }

  # Both hits and revhits initially have locusId of the form orgId:locusId
  # This is converted to an orgId and a locusId below
  my @hitFields = qw{locusId type curatedId bits locusBegin locusEnd cBegin cEnd cLength identity};
  my @hits = ReadTable($hitsFile, \@hitFields);
  # add the fromCurated field to all blast hits -- these are lower priority because
  # the reference protein is not actually characterized
  foreach my $hit (@hits) {
    if ($hit->{type} eq "blast") {
      $hit->{fromCurated} = $hit->{curatedId} =~ m/^curated2:/ ? 1 : 0;
    }
  }
  my @revFields = qw{locusId otherId bits locusBegin locusEnd otherBegin otherEnd otherIdentity};
  my @revhits = ReadTable($revhitsFile, \@revFields);

  my %curatedInfo = (); # id => length, desc
  open(my $fhInfo, "<", $infoFile) || die "Cannot read $infoFile\n";
  while (my $line = <$fhInfo>) {
    chomp $line;
    my ($ids, $length, $descs) = split /\t/, $line;
    die "Not enough columns in $infoFile: $line\n" unless defined $descs;
    $curatedInfo{$ids} = [$length, $descs] unless $ids eq "ids";
  }
  close($fhInfo) || die "Error reading $infoFile\n";

  my %hits; # orgId => pathwayId => step => locusId => list of relevant hits
  foreach my $hit (@hits) {
    my $curatedId = $hit->{curatedId};
    my ($orgId, $locusId) = ParseOrgLocus($hit->{locusId});
    $hit->{orgId} = $orgId;
    $hit->{locusId} = $locusId;
    my $hash; # pathway => step => 1
    if ($hit->{type} eq "hmm") { # curatedId is an hmm id
      die "Unknown hmm query $curatedId found in hits file $hitsFile\n"
        unless exists $hmmToStep{$curatedId};
      $hash = $hmmToStep{$curatedId};
    } else {
      die "Unknown hit type $hit->{type} in hits file $hitsFile\n"
        unless $hit->{type} eq "blast";
      die "Unknown blast query $curatedId found in hits file $hitsFile\n"
        unless exists $queryToStep{$curatedId};
      $hash = $queryToStep{$curatedId};
    }
    while (my ($pathwayId, $stephash) = each %$hash) {
      foreach my $step (keys %$stephash) {
        push @{ $hits{$orgId}{$pathwayId}{$step}{$locusId} }, $hit;
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

  # add locusLength to each hit (it already has cLength)
  # also add coverage (of the curated item)
  foreach my $hit (@hits) {
    my $orgId = $hit->{orgId};
    my $locusId = $hit->{locusId};
    die "Unknown locus $orgId $locusId in $hitsFile but not in $aaIn\n"
      unless exists $locusInfo{$orgId}{$locusId};
    my $len = $locusInfo{$orgId}{$locusId}[0];
    $hit->{locusLength} = $len;
    $hit->{coverage} = ($hit->{cEnd} - $hit->{cBegin} + 1) / $hit->{cLength};
  }

  # add otherLength and otherCoverage to each revhit
  foreach my $revhit (@revhits) {
    my $otherId = $revhit->{otherId};
    die "Unknown curated id $otherId in rev hits $revhitsFile\n"
      unless exists $curatedInfo{$otherId};
    my $len = $curatedInfo{$otherId}[0];
    $revhit->{otherLength} = $len;
    $revhit->{otherCoverage} = ($revhit->{otherEnd} - $revhit->{otherBegin} + 1) / $len;
  }

  # Score each candidate based on its various hits, the reverse hits,
  # and the ignore information
  my %cand = (); # orgId => pathwayId => step => sorted list of candidates
  while (my ($orgId, $pathhash) = each %hits) {
    while (my ($pathwayId, $stephash) = each %$pathhash) {
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
                         || $a->{fromCurated} <=> $b->{fromCurated}
                         || $b->{maxBits} <=> $a->{maxBits}
                         || $b->{locusId} cmp $a->{locusId} } @cand;
        @cand = splice(@cand, 0, $maxCand) if @cand > $maxCand;

        # look for split ORFs unless there are two high-confidence candidates
        unless (@cand >= 2 && $cand[1]{score} == 2) {
          my %locusCand = map { $_->{locusId} => $_ } @cand;
          my $merge = FindSplit($locushash, \%locusRev, \%locusCand);
          if (defined $merge) {
            # Put the HMM and fromCurated info from 1st (higher scoring) locus into $merge
            my $cand1 = $locusCand{ $merge->{locusId} };
            die unless defined $cand1;
            foreach my $key (qw{hmmBits hmmId hmmCoverage hmmScore fromCurated}) {
              $merge->{$key} = $cand1->{$key};
            }
            $merge->{maxBits} = max($merge->{blastBits}, $merge->{hmmBits} || 0);
            # replace the previous entries for the merged genes with the merge
            @cand = grep { $_->{locusId} ne $merge->{locusId}
                             && $_->{locusId} ne $merge->{locusId2} } @cand;
            push @cand, $merge;
            # Re-sort and re-truncate the list
            @cand = sort { $b->{score} <=> $a->{score}
                             || $a->{fromCurated} <=> $b->{fromCurated}
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

  # And score the rules
  my %ruleScores = (); # orgId => pathway => rulename => hash of n, n012, score, path, path2, stepsUsed, expandedPath
  # where n is #steps, n012 is a vector with #of steps at each score,
  # score is the total weighted score,
  # path and path2 are lists of steps (path2 is used only if there is a tie),
  # and stepsUsed is a hash that includes all steps used, including steps for dependencies
  # expandedPath is the list of all steps including those in dependencies (with duplicates removed)
  # Note that n and score need to be recomputed each time because of potential overlap
  foreach my $orgId (sort keys %orgs) {
    foreach my $pathId (@pathways) {
      my $st = $pathways{$pathId};
      my $steps = $st->{steps};
      my $rules = $st->{rules};
      foreach my $rule (@{ $st->{ruleOrder} }) {
        my @scoredPaths = ();
        foreach my $path (@{ $rules->{$rule} }) {
          my %stepsUsed = (); # organism dependent due to scoring of dependencies
          foreach my $piece (@$path) {
            if (exists $steps->{$piece}) {
              $stepsUsed{$piece} = 1;
            } elsif (exists $ruleScores{$orgId}{$pathId}{$piece}) {
              foreach my $step (keys %{ $ruleScores{$orgId}{$pathId}{$piece}{stepsUsed} }) {
                $stepsUsed{$step} = 1;
              }
            } else {
              die "No score yet for dependency $piece while scoring $rule for pathway $pathId\n";
            }
          }
          # score this path
          my $n = 0;
          my @n012 = (0,0,0); # number of steps with score of 0, 1, or 2
          my $scoreTot = 0;
          foreach my $step (keys %stepsUsed) {
            my $cand = $cand{$orgId}{$pathId}{$step}[0]
              if exists $cand{$orgId}{$pathId}{$step} && @{ $cand{$orgId}{$pathId}{$step} } > 0;
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
        $ruleScores{$orgId}{$pathId}{$rule} = $scoredPaths[0];
        $ruleScores{$orgId}{$pathId}{$rule}{path2} = $scoredPaths[1]{path} if @scoredPaths >= 2;
        # and compute pathExpanded for the main path
        my $path = $ruleScores{$orgId}{$pathId}{$rule}{path};
        my @pathExpanded = ();
        foreach my $piece (@$path) {
          if (exists $steps->{$piece}) {
            push @pathExpanded, $piece;
          } else {
            my $pathPart = $ruleScores{$orgId}{$pathId}{$piece}{pathExpanded};
            die unless defined $pathPart;
            push @pathExpanded, @$pathPart;
          }
        }
        my %pathSoFar = ();
        @pathExpanded = grep { my $keep = !exists $pathSoFar{$_};
                               $pathSoFar{$_} = 1;
                               $keep; } @pathExpanded;
        $ruleScores{$orgId}{$pathId}{$rule}{pathExpanded} = \@pathExpanded;
      }
    }
  }

  # Record which steps are on the best path for the 'all' rule
  my %onBestPath = (); # orgId => pathId => step => 1
  foreach my $orgId (sort keys %orgs) {
    foreach my $pathId (@pathways) {
      my $path = $ruleScores{$orgId}{$pathId}{'all'}{pathExpanded};
      die "No path for all for $orgId $pathId" unless @$path > 0;
      foreach my $step (@$path) {
        $onBestPath{$orgId}{$pathId}{$step} = 1;
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
                      otherIds otherBits otherIdentity otherCoverage otherDesc};
  print $fhCand join("\t", @candfields)."\n";
  my @stepfields = qw{orgId gdb gid pathway step onBestPath score locusId sysName score2 locusId2 sysName2};
  print $fhStep join("\t", @stepfields)."\n";
  foreach my $orgId (sort keys %cand) {
    foreach my $pathId (@pathways) {
      my $steps = $pathways{$pathId}{steps};
      my $stephash = $cand{$orgId}{$pathId} || {};
      foreach my $step (sort keys %$steps) {
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
          $cand->{hmmDesc} = $queryDesc{ $cand->{hmmId} } if $cand->{hmmId};
          $cand->{otherDesc} = $curatedInfo{ $cand->{otherIds } }[1] if $cand->{otherIds};
          $cand->{gdb} = $orgs{$orgId}{gdb};
          $cand->{gid} = $orgs{$orgId}{gid};
          my @out = map { defined $cand->{$_} ? $cand->{$_} : "" } @candfields;
          print $fhCand join("\t", @out)."\n";
        }
        # per-step line
        my @stepout = ($orgId, $orgs{$orgId}{gdb}, $orgs{$orgId}{gid},
                       $pathId, $step, exists $onBestPath{$orgId}{$pathId}{$step} ? 1 : 0);
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

  # And write out the best paths
  open (my $fhR, ">", "$outpre.rules") || die "Cannot write to $outpre.rules\n";
  print $fhR join("\t", qw{orgId gdb gid pathway rule nHi nMed nLo score expandedPath path path2})."\n";
  foreach my $orgId (sort keys %orgs) {
    foreach my $pathId (@pathways) {
      my $st = $pathways{$pathId};
      foreach my $rule (@{ $st->{ruleOrder} }) {
        my $sc = $ruleScores{$orgId}{$pathId}{$rule};
        die unless exists $sc->{n};
        print $fhR join("\t", $orgId, $orgs{$orgId}{gdb}, $orgs{$orgId}{gid},
                        $pathId, $rule,
                        $sc->{n012}[2], $sc->{n012}[1], $sc->{n012}[0], $sc->{score},
                        join(" ", @{$sc->{pathExpanded}}),
                        join(" ", @{$sc->{path}}),
                        exists $sc->{path2} ? join(" ", @{$sc->{path2}}) : "") . "\n";
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

  my @blastHits = sort { $b->{bits} <=> $a->{bits} || $a->{curatedId} cmp $b->{curatedId} }
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
    $score = 1 if $score == 2 && $hit->{curatedId} =~ m/^curated2:/;
    unless (defined $bestBlastHit && $score <= $bestBlastHit->{blastScore}) {
      $bestBlastHit = { 'blastBits' => $hit->{bits},
                        'curatedIds' => $hit->{curatedId},
                        'fromCurated' => $hit->{fromCurated},
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
                      'hmmId' => $hit->{curatedId},
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
  $out->{fromCurated} = 0 unless defined $out->{fromCurated};
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
  my %locusTarget = (); # locusId => curatedId => best hit object
  while (my ($locusId, $hits) = each %$hithash) {
    my @hits = grep { $_->{type} eq "blast" && $_->{identity} >= 30 } @$hits;
    @hits = sort { $a->{fromCurated} <=> $b->{fromCurated}
                     || $b->{bits} <=> $a->{bits} } @hits;
    @hits = splice(@hits, 0, $maxCand) if @hits > $maxCand;
    $locusBest{$locusId} = \@hits if @hits > 0;
    foreach my $hit (@hits) {
      $locusTarget{$locusId}{ $hit->{curatedId} } = $hit
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

      my $curatedId = $hit2->{curatedId} || die;
      for (my $i = 0; $i < $j; $i++) {
        my $locus1 = $loci[$i]; # the higher-scoring one
        die if $locus1 eq $locus2;
        next unless exists $locusTarget{$locus1}{$curatedId};
        my $hit1 = $locusTarget{$locus1}{$curatedId};
        die unless $hit1 && $hit2 && $hit1->{curatedId} eq $hit2->{curatedId};
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
  # redundant (i.e., same pair of loci but different curatedId)
  return $comb[0]; # may be undef
}

sub MergeHits($$$$$$){
  my ($hit1, $bestscore1, $revhits1,
      $hit2, $bestscore2, $revhit2) = @_;
  die unless $hit1->{curatedId} eq $hit2->{curatedId};
  my $curatedId = $hit1->{curatedId};

  my $locus1 = $hit1->{locusId};
  my $locus2 = $hit2->{locusId};
  die unless defined $locus1 && defined $locus2 && $locus1 ne $locus2;

  my $alnLen1 = $hit1->{cEnd} - $hit1->{cBegin}+1;
  my $alnLen2 = $hit2->{cEnd} - $hit2->{cBegin}+1;

  # Overlap on the curated sequence should be at most 20% of either alignment
  # Combined alignment should cover at least 70% of the curated protein
  my ($overlap, $combLen);
  if ($hit1->{cEnd} < $hit2->{cBegin} || $hit1->{cBegin} > $hit2->{cEnd}) {
    # No overlap
    $overlap = 0;
    $combLen = $alnLen1 + $alnLen2;
  } else {
    my $combBegin = min($hit1->{cBegin}, $hit2->{cBegin});
    my $combEnd = max($hit1->{cEnd}, $hit2->{cEnd});
    $combLen = $combEnd - $combBegin + 1;
    $overlap = $alnLen1 + $alnLen2 - $combLen;
  }
  return undef unless $overlap < 0.2 * $alnLen1
    && $overlap < 0.2 * $alnLen2
    && $combLen >= 0.7 * $hit1->{cLength};


  # Check that there is no better other-hit for locus1 and locus2
  return undef if defined $revhit2 && $hit2->{bits} <= $revhit2->{bits} + 10;

  my $revhit1 = RelevantRevhit($hit1->{locusBegin}, $hit1->{locusEnd}, $revhits1);
  return undef if defined $revhit1 && $hit1->{bits} <= $revhit1->{bits} + 10;

  # Score the hit and check that it is better than either component
  my $combIdentity = ($alnLen1 * $hit1->{identity} + $alnLen2 * $hit2->{identity}) / ($alnLen1 + $alnLen2);
  my $combCoverage = $combLen / $hit1->{cLength};
  my $combScore = $combIdentity >= 40 && $combCoverage >= 0.8 ? 2 : 1;
  $combScore = 1 if $combScore == 2 && $curatedId =~ m/^curated2:/;

  return undef unless $combScore > $bestscore1 && $combScore > $bestscore2;

  # Ensure that the hash references below succeed
  $revhit1 = {} if !defined $revhit1;
  $revhit2 = {} if !defined $revhit2;

  # Return a merged hit
  return { 'locusId' => $locus1,
           'locusId2' => $locus2,
           'curatedIds' => $curatedId,
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

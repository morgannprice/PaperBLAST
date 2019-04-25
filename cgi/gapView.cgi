#!/usr/bin/perl -w

# Parameters:
# Required: base=subdir/orgs -- where the input files are, relative to tmp
#
# Optional:
# orgId -- which organism
# path -- which pathway
# step -- which step in that pathway

use strict;
use CGI qw(:standard Vars start_ul);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use HTML::Entities;
use lib "../lib";
use Steps;
use pbutils;
use pbweb qw{start_page};

sub ScoreToStyle($);
sub ScoreToLabel($);
sub ShowScoreShort($);
sub HMMToURL($);

my $base = param("base") || die "Must specify base";
$base =~ m!^[a-zA-Z0-9._-]+/?[a-zA-Z0-9._-]+$! || die "Invalid base $base";
my $pre = "../tmp/$base";

my @orgs = ReadOrgTable("$pre.org");
die "No organisms for $pre.org" unless @orgs > 0;
my %orgs = map { $_->{orgId} => $_ } @orgs;
my $faafile = "$pre.faa";
die "No such file: $faafile" unless -e $faafile;

my $pathSpec = param("path");
$pathSpec = "" if !defined $pathSpec;
die "Invalid path parameter $pathSpec"
  unless $pathSpec eq "" || $pathSpec =~ m/^[a-zA-Z0-9._-]+$/;
my @path = (); # all pathways, or the one specified

if ($pathSpec eq "") {
  @path = sort qw{his met ser thr gly cys ile leu val chorismate phe tyr trp asn gln arg lys pro};
} else {
  @path = ($pathSpec);
}

foreach my $path (@path) {
  my $sumfile = "$pre.$path.sum.steps";
  die "No such file: $sumfile" unless -e $sumfile;
}

my $orgId = param("orgId"); # the organism specified, or, ""
$orgId = "" if !defined $orgId;
if ($orgId ne "") {
  die "Unknown orgId $orgId" unless exists $orgs{$orgId};
}

my $step = param("step");
if (!defined $step) {
  $step = "";
} else {
  $step =~ m/^[a-zA-Z0-9._'"-]+$/ || die "Invalid step $step";
}

my ($steps, $rules);
if (@path == 1) {
  my $st = ReadSteps("../tmp/gaps/$path[0].steps");
  $steps = $st->{steps};
  $rules = $st->{rules};
  die "Non-existent step $step" if $step ne "" && !exists $steps->{$step};
}

my $title = "Gaps";
$title .= " for $pathSpec" if $pathSpec ne "";
$title = "Finding step $step for $pathSpec"
  if $step ne "" && $orgId ne "" && $pathSpec ne "";
$title .= " in $orgs{$orgId}{genomeName}"
  if $orgId ne "";
start_page('title' => $title, 'banner' => 'Gap viewer', 'bannerURL' => "gapView.cgi?base=$base");

if ($orgId eq "") {
  my @orgsSorted = sort { $a->{genomeName} cmp $b->{genomeName} } @orgs;
  print p(scalar(@orgsSorted), "organisms");
  print start_ul;
  foreach my $org (@orgsSorted) {
    my $orgId = $org->{orgId};
    my $URL = "gapView.cgi?base=$base&orgId=$orgId";
    $URL .= "&path=$pathSpec" if @path == 1;
    print li(a({ -href => $URL }, $org->{genomeName} ));
  }
  print end_ul;
} elsif (@path > 1) {
  # overview of pathways for this organism
  print start_ul;
  foreach my $path (@path) {
    print li(a({-href => "gapView.cgi?base=$base&orgId=$orgId&path=$path"}, $path));
  }
  print end_ul;
  print p("Or see",
          a({ -href => "gapView.cgi?base=$base"},
            scalar(@orgs), "organisms"))
    if @orgs > 1;
} elsif ($step eq "") {
  # overview of this pathway, first the rules, then all of the steps
  my @sumRules = ReadTable("$pre.$pathSpec.sum.rules", qw{orgId gdb gid rule score expandedPath});
  @sumRules = grep { $_->{orgId} eq $orgId } @sumRules;
  my @sumSteps = ReadTable("$pre.$pathSpec.sum.steps", qw{orgId gdb gid step score locusId sysName});
  @sumSteps = grep { $_->{orgId} eq $orgId } @sumSteps;
  my %sumSteps = map { $_->{step} => $_ } @sumSteps;
  print h3(scalar(@sumRules), "rules");
  print start_ul;
  foreach my $rule (@sumRules) {
    my @stepList = split / /, $rule->{expandedPath};
    my @parts = ();
    foreach my $step (@stepList) {
      my $stepDef = $steps->{$step} || die "Invalid step $step";
      my $stepS = exists $sumSteps{$step} ? $sumSteps{$step} : {};
      my $score = $stepS->{score} || 0;
      my $label = ScoreToLabel($score);
      my $id = $stepS->{sysName} || $stepS->{locusId} || "";
      push @parts, a({ -href => "gapView.cgi?base=$base&orgId=$orgId&path=$pathSpec&step=$step",
                       -style => ScoreToStyle($score),
                       -title => "$stepDef->{desc} -- $id ($label)" },
                     $step);
    }
    print li($rule->{rule}.":", @parts);
  }
  print end_ul;
  my @stepsSorted = sort { $a->{i} <=> $b->{i} } (values %$steps);
  print h3(scalar(@stepsSorted) . " steps (" . scalar(@sumSteps) . " with candidates)");
  my %sumSteps = map { $_->{step} => $_ } @sumSteps;
  my @tr = ();
  my @header = qw{Step Description Best 2nd-best};
  push @tr, Tr(th(\@header));
  # For each step, show the step name and description, the best candidate (if any), and the 2nd best candidate(s) if any
  # Use all of the steps that are defined, not just the ones in the rules file, as ones with no candidates are missing
  foreach my $stepS (@stepsSorted) {
    my $step = $stepS->{name};
    die "invalid step $step" unless exists $steps->{$step};
    my @cand = ();
    if (exists $sumSteps{$step}) {
      push @cand, [ $sumSteps{$step}{locusId}, $sumSteps{$step}{sysName}, $sumSteps{$step}{score} ]
        if $sumSteps{$step}{locusId} ne "";
      push @cand, [ $sumSteps{$step}{locusId2}, $sumSteps{$step}{sysName2}, $sumSteps{$step}{score2} ]
        if $sumSteps{$step}{locusId2} ne "";
    }
    my @show = ();
    foreach my $cand (@cand) {
      my ($locusId,$sysName,$score) = @$cand;
      my $id = $sysName || $locusId;
      push @show, span({-style => ScoreToStyle($score), -title => ScoreToLabel($score) }, $id);
    }
    while(@show < 2) {
      push @show, "";
    }
    push @tr, Tr(td(a({-href => "gapView.cgi?base=$base&orgId=$orgId&path=$pathSpec&step=$step"}, $step),
                    td($stepS->{desc}),
                    td($show[0]), td($show[1])));
  }
  print table({-cellpadding=>2, -cellspacing=>0, -border=>1}, @tr);
} else {
  # overview of this step in this organism
  my @req = qw{orgId gdb gid step score locusId sysName desc locusId2 sysName2 desc2
               blastBits curatedIds identity blastCoverage blastScore curatedDesc
               hmmBits hmmId hmmCoverage hmmScore hmmDesc
               otherIds otherBits otherIdentity otherCoverage otherDesc};
  my @cand = ReadTable("$pre.$pathSpec.sum.cand", \@req);
  @cand = grep { $_->{orgId} eq $orgId && $_->{step} eq $step } @cand;
  if (@cand == 0) {
    print h3("No candidates for $step: $steps->{$step}{desc}");
  } else {
    print h3(scalar(@cand), "candidates for $step:", $steps->{$step}{desc});
    my @header = qw{Score Gene Description Similar-to Id. Cov. Bits  Other-hit Other-id. Other-bits};
    $header[-1] = span({-title => "A characterized protein that is similar to the gene but is not associated with step $step"},
                       $header[-1]);
    foreach (@header) { s/-/ /; }
    my @tr = Tr(th({-valign => "bottom"}, \@header));
    foreach my $cand (@cand) {
      # potentially make two rows, one for BLAST and one for HMM
      my $id = $cand->{sysName} || $cand->{locusId};
      my $desc = $cand->{desc};
      if ($cand->{locusId2}) { # (this should only happen for BLAST hits)
        $id .= "; " . ($cand->{sysName2} || $cand->{locusId2});
        $desc .= "; " . $cand->{desc2};
      }
      if ($cand->{blastScore} ne "") {
        my $otherIdentity = "";
        $otherIdentity = span({ -title => "coverage: " . int(0.5 + 100 *$cand->{otherCoverage})."%"},
                              int(0.5 + $cand->{otherIdentity})."%")
          if $cand->{otherBits};
        my @hitIds = split /,/, $cand->{curatedIds};
        my $URL = "http://papers.genomics.lbl.gov/cgi-bin/litSearch.cgi?query=" . $hitIds[0];
        my $idShowHit = $hitIds[0];
        $idShowHit =~ s/^.*://;
        my @otherIds = split /,/, $cand->{otherIds};
        my $URLother = "";
        my $idShowOther = "";
        if ($cand->{otherBits}) {
          $URLother = "http://papers.genomics.lbl.gov/cgi-bin/litSearch.cgi?query=" . $otherIds[0];
          $idShowOther = $otherIds[0];
          $idShowOther =~ s/^.*://;
        }
        my $descShowCurated = $cand->{curatedDesc}; $descShowCurated =~ s/;;.*//;
        my $descShowOther = $cand->{otherDesc}; $descShowOther =~ s/;;.*//;
        push @tr, Tr(td({-valign => "top"},
                         [ ShowScoreShort($cand->{blastScore}),
                           $id, $desc,
                          a({-href => $URL, -title => "View $idShowHit in PaperBLAST"}, $descShowCurated),
                          int(0.5 + $cand->{identity})."%",
                          int(0.5 + 100 * $cand->{blastCoverage})."%",
                          $cand->{blastBits},
                          small($cand->{otherBits} ? a({-href => $URLother,
                                                  -title => "view $idShowOther in PaperBLAST"},
                                                 $descShowOther) : ""),
                          small($otherIdentity),
                          small($cand->{otherBits} > 0 ? $cand->{otherBits} : "")
                         ]));
      }
      if ($cand->{hmmScore} ne "") {
        my $hmmURL = HMMToURL($cand->{hmmId});
        push @tr, Tr(td({ -valign => "top" },
                        [ ShowScoreShort($cand->{hmmScore}), $id, $desc,
                          a({-href => $hmmURL, -title => $cand->{hmmDesc} }, $cand->{hmmId}),
                          "", int(0.5 + 100 * $cand->{hmmCoverage})."%",
                          $cand->{hmmBits}, "", "", ""]));
      }
    }
    print table({-cellpadding=>2, -cellspacing=>0, -border=>1}, @tr);
  }
  # Arguably should be showing the queries as well
  print h3("Definition of $step");
  print start_ul();
  foreach my $search (@{ $steps->{$step}{search} }) {
    my ($type,$value) = @$search;
    my $show;
    if ($type eq "EC") {
      my $URL = "http://papers.genomics.lbl.gov/cgi-bin/genomeSearch.cgi?gdb=" . $orgs{$orgId}{gdb}
        . "&gid=" . $orgs{$orgId}{gid}
        . "&query=$value&word=1";
      $show = "Curated proteins or TIGRFams with EC " . a({-href => $URL, -title => "Run Curated BLAST"}, $value);
    } elsif ($type eq "hmm") {
      $show = "HMM " . a({-href => HMMToURL($value) }, $value);
    } elsif ($type eq "term") {
      my $URL = "http://papers.genomics.lbl.gov/cgi-bin/genomeSearch.cgi?gdb=" . $orgs{$orgId}{gdb}
        . "&gid=" . $orgs{$orgId}{gid}
        . "&word=1"
        . "&query=" . encode_entities($value);
      $show = "Curated proteins matching " . a({-href => $URL, -title => "Run Curated BLAST"}, $value);
    } elsif ($type eq "curated") {
      my $URL = "http://papers.genomics.lbl.gov/cgi-bin/litSearch.cgi?query=".$value;
      my $show_id = $value; $show_id =~ s/^.*://;
      $show = "Curated sequence " . a({-href => $URL, -title => "View in PaperBLAST"}, $show_id);
    } elsif ($type eq "uniprot") {
      my $URL = "https://www.uniprot.org/uniprot/".$value;
      $show = "UniProt sequence" . a({-href => $URL, -title => "View in UniProt"}, $value);
    } elsif ($type eq "ignore_other") {
      my $URL = "http://papers.genomics.lbl.gov/cgi-bin/curatedSearch.cgi?word=1"
        . "&query=" . encode_entities($value);
      $show = "Ignore hits to items matching "
        . a({-href => $URL}, $value)
          . " when looking for 'other' hits";
    } elsif ($type eq "ignore") {
      my $URL = "http://papers.genomics.lbl.gov/cgi-bin/litSearch.cgi?query=$value";
      my $showId = $value; $showId =~ s/^.*://;
      $show = "Ignore hits to "
        . a({-href => $URL, -title => "View in PaperBLAST"}, $showId)
        . " when looking for 'other' hits";
    }
    print li($show);
  }
  print end_ul();
}

print p("Or see all pathways for",
        a({ -href => "gapView.cgi?base=$base&orgId=$orgId" }, $orgs{$orgId}{genomeName}))
  if $orgId ne "" && $pathSpec ne "" && $step eq "";
print p("Or see all steps for",
        a({ -href => "gapView.cgi?base=$base&orgId=$orgId&path=$pathSpec" }, $pathSpec),
        "in",
        a({ -href => "gapView.cgi?base=$base&orgId=$orgId" }, $orgs{$orgId}{genomeName}))
  if $orgId ne "" && $pathSpec ne "" && $step ne "";

print end_html;

sub ScoreToStyle($) {
  my ($score) = @_;
  my $color = $score > 1 ? "#007000" : ($score < 1 ? "#CC4444" : "#000000");
  return "color: $color; font-weight: bold;" if $score > 1;
  return "color: $color;";
}

sub ScoreToLabel($) {
  my ($score) = @_;
  return $score > 1 ? "high confidence" : ($score < 1 ? "low confidence" : "medium confidence");
}

sub ShowScoreShort($) {
  my ($score) = @_;
  return span({ -style => ScoreToStyle($score), -title => ScoreToLabel($score) },
              $score > 1 ? "hi" : ($score < 1 ? "lo" : "med"));
}

sub HMMToURL($) {
  my ($hmmId) = @_;
  if ($hmmId =~ m/^TIGR/) {
    return "http://tigrfams.jcvi.org/cgi-bin/HmmReportPage.cgi?acc=".$hmmId;
  } elsif ($hmmId =~ m/^PF/) {
    my $hmmIdShort = $hmmId; $hmmIdShort =~ s/[.]\d+$//;
    return "http://pfam.xfam.org/family/".$hmmIdShort;
  }
  return "";
}

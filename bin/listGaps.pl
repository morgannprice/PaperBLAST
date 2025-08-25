#!/usr/bin/perl -w
use strict;

use strict;
use Getopt::Long;
use FindBin qw{$RealBin};
use lib "$RealBin/../lib";
use pbutils qw{ReadTable};
use DBI;

my $stepBase = "$RealBin/../gaps";
my $resultsDir = "$RealBin/../tmp";
my ($stepDir, $queryDir);

my $usage = <<END
checkCuratedGaps.pl -set aa -orgs orgsDef.208 orgsDef.208.6f -out prefix

Given a pathway set, an org directory with computed results, and a
6-frame directory with computed results, build a table of known gaps
and spurious ("curated") gaps, written to prefix.curated.gaps.tsv and
prefix.known.gaps.tsv
END
  ;

my $set;
my @orgSpecs;
my $outPre;

die $usage
  unless GetOptions('set=s' => \$set, 'orgs=s{2,2}' => \@orgSpecs,
                    'out=s' => \$outPre)
  && @ARGV == 0
  && defined $set && defined $outPre && @orgSpecs == 2;

my @orgDirs = map $resultsDir . "/" . $_, @orgSpecs;
foreach my $dir (@orgDirs) {
  die "No such directory: $dir\n" unless -d $dir;
  die "No such file: $dir/$set.sum.steps\n" unless -e "$dir/$set.sum.steps";
}
my ($dir, $dir6) = @orgDirs;

my @orgs = ReadTable("$dir/orgs.org", ["orgId", "gdb","gid","genomeName"]);
my @orgs6 = ReadTable("$dir6/orgs.org", ["orgId", "gdb","gid","genomeName"]);
die "Different numbers of organisms for $dir and $dir6\n"
  unless scalar(@orgs) == scalar(@orgs6);
my %orgs = map { $_->{orgId} => $_ } @orgs;

my @stepHeader = qw{orgId gdb gid pathway step onBestPath score};
my @steps = ReadTable("$dir/$set.sum.steps", \@stepHeader);
my @steps6 = ReadTable("$dir6/$set.sum.steps", \@stepHeader);

my %steps6 = (); # orgId => pathway => step => score
my %steps6best = (); # orgId => pathway => step => score but only for best path
foreach my $row (@steps6) {
  $steps6{ $row->{orgId} }{ $row->{pathway} }{ $row->{step} } = $row->{score};
  $steps6best{ $row->{orgId} }{ $row->{pathway} }{ $row->{step} } = $row->{score}
    if $row->{onBestPath};
}

my %nogap6 = (); # orgId => pathway if the best path has no low-confidence steps
foreach my $orgId (sort keys %steps6best) {
  my $hash = $steps6best{$orgId};
  foreach my $pathway (sort keys %$hash) {
    my @scores = values %{ $hash->{$pathway} };
    my @scores0 = grep { $_ ne "1" && $_ ne "2" } @scores;
    $nogap6{$orgId}{$pathway} = 1 if scalar(@scores0) == 0;
  }
}

my $outCurated = "$outPre.curated.gaps.tsv"; # spurious gaps
my $outKnown = "$outPre.known.gaps.tsv"; # actual gaps
open(my $fhCurated, ">", $outCurated) || die "Cannot write to $outCurated\n";
open(my $fhKnown, ">", $outKnown) || die "Cannot write to $outKnown\n";

print $fhCurated join("\t", qw{gdb gid genomeName pathway step class comment}) . "\n";
print $fhKnown join("\t", qw{gdb gid genomeName pathway step}) . "\n";

my $nGenuineGap = 0;
my $nSpuriousSame = 0;
my $nSpuriousOther = 0;
foreach my $row (@steps) {
  my $orgId = $row->{orgId};
  die "Unknown orgId $orgId" unless exists $orgs{$orgId};
  my $org = $orgs{$orgId};
  next unless $row->{onBestPath};
  next if $row->{score} eq "1" || $row->{score} eq "2"; # not a gap
  die "No sixframe analysis for $orgId" unless exists $steps6{ $orgId };
  my $score6 = $steps6{ $row->{orgId} }{ $row->{pathway} }{ $row->{step} };
  if ($score6 eq "1" || $score6 eq "2") {
    # spurious
    $nSpuriousSame++;
    my $confString = $score6 eq "2" ? "it was found" : "a medium-confidence candidate was found";
    print $fhCurated join("\t", $org->{gdb}, $org->{gid}, $org->{genomeName},
                          $row->{pathway}, $row->{step}, "spurious",
                          "Although this step was not found in the annotated proteins, $confString when analyzing the six-frame translation"
                          . " of the genome. This may indicate a missing gene call"
                          . " or a frameshift error in the genome sequence.")."\n";
  } elsif ($nogap6{ $row->{orgId} }{ $row->{pathway} }) {
    # spurious but different step
    $nSpuriousOther++;
    print $fhCurated join("\t", $org->{gdb}, $org->{gid}, $org->{genomeName},
                          $row->{pathway}, $row->{step}, "spurious",
                          "Although this step was not found in the annotated proteins, when analyzing the six-frame translation of the genome,"
                          . " another step was found to complete the pathway. This may indicate a missing gene call"
                          . " or a frameshift error in the genome sequence.")."\n";
  } else {
    $nGenuineGap++;
    print $fhKnown join("\t", $org->{gdb}, $org->{gid}, $org->{genomeName},
                        $row->{pathway}, $row->{step})."\n";
  }
}
close($fhCurated) || die "Error writing to $outCurated\n";
close($fhKnown) || die "Error writing to $outKnown\n";
print STDERR "Wrote $outCurated\nWrote $outKnown\n";
print STDERR "Total low-confidence steps: genuine gaps $nGenuineGap spurious (step found in 6-frame) $nSpuriousSame other-spurious (other step found in 6-frame) $nSpuriousOther\n";

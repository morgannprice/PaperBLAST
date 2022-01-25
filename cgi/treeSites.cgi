#!/usr/bin/perl -w
# treeSites.cgi -- show key sites along a tree
use strict;
use CGI qw(:standard Vars start_ul);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use URI::Escape;
use HTML::Entities;
use List::Util qw{sum max};
use lib "../lib";
use pbutils qw{ReadFastaEntry};
use MOTree;

# In rendering mode, these options are required:
# aln or alnFile -- alignment in multi-fasta format
#   In header lines, anything after the initial space is assumed to be a description.
#   Either "." or "-" are gap characters.
# tree or treeFile -- tree in newick format.
#   Labels on internal nodes are ignored.
#   Multi-furcations are allowed. It is treated as rooted.
# Optional arguments:
# anchor -- sequence to choose positiosn from
# pos -- comma-delimited list of positions (in anchor if set, or else, in alignment)
#
# If alnFile or treeFile is missing, it shows an input form.

sub fail($);
sub warning($);

# maximum size of posted data, in bytes
my $maxMB = 100;
$CGI::POST_MAX = $maxMB*1024*1024;

print
  header(-charset => 'utf-8'),
  start_html(-title => "Sites on a Tree"),
  h2("Sites on a Tree"),
  "\n";

my $alnSet = (defined param('aln') && param('aln') ne "")
  || (defined param('alnFile') && ref param('alnFile'));
my $treeSet = (defined param('tree') && param('tree') ne "")
  || (defined param('treeFile') && ref param('treeFile'));

if (!$alnSet || !$treeSet) {
  # show the form
  print
    start_form(-name => 'input', -method => 'POST', -action => 'treeSites.cgi'),
    p("Alignment in multi-fasta format:",
      br(),
      textarea(-name => 'aln', -value => '', -cols => 70, -rows => 10),
      br(),
      "Or upload:", filefield(-name => 'alnFile', -size => 50)),
    p("Rooted tree in newick format:",
      br(),
      textarea(-name => 'tree', -value => '', -cols => 70, -rows => 4),
      br(),
      "Or upload:", filefield(-name => 'treeFile', -size => 50)),
    p("Select positions to show:",
      textfield(-name => 'pos', -size => 50),
      br(), br(),
      "from sequence:", textfield(-name => 'anchor', -size => 25)),
    p(submit()),
    end_form;
  print end_html;
  exit(0);
}

# else rendering mode
my $fhAln;
if (defined param('aln') && param('aln') ne "") {
  my $aln = param('aln');
  open $fhAln, "<", \$aln;
} elsif (defined param('alnFile') && ref param('alnFile')) {
  $fhAln = param('alnFile')->handle;
  die unless $fhAln;
} else {
  fail("No alignment specified");
}

my %alnSeq = ();
my %alnDesc = ();
my $state = {};
while (my ($id, $seq) = ReadFastaEntry($fhAln, $state)) {
  if ($id =~ m/^(\S+)\s(.*)/) {
    $id = $1;
    $alnDesc{$id} = $2;
  }
  fail("Duplicate sequence for " . encode_entities($id))
    if exists $alnSeq{$id};
  $seq =~ s/[.]/-/g;
  $seq = uc($seq);
  $seq =~ m/^[A-Z-]+$/ || fail("Invalid characters in sequence for " . encode_entities($id));
  $alnSeq{$id} = $seq;
}
fail("No sequences in the alignment")
  if (scalar(keys %alnSeq) == 0);

my $alnLen;
while (my ($id, $seq) = each %alnSeq) {
  $alnLen = length($seq) if !defined $alnLen;
  fail("Inconsistent sequence length for " . encode_entities($id))
    unless length($seq) == $alnLen;
}

my $moTree;
if (defined param('tree') && param('tree') ne "") {
  eval { $moTree = MOTree::new('newick' => param('tree')) };
} elsif (defined param('treeFile') && ref param('treeFile')) {
  my $fh = param('treeFile')->handle;
  die unless $fh;
  eval { $moTree = MOTree::new('fh' => $fh) };
} else {
  fail("No tree specified");
}
fail("Could not parse tree") unless $moTree;

# Fail if there are any leaves in the tree without sequences
my $nodes = $moTree->depthfirst(); # depth-first ordering of all nodes
my $root = $moTree->get_root_node;
my @leaves = grep $moTree->is_Leaf($_), @$nodes;
my %idToLeaf = ();
foreach my $leaf (@leaves) {
  my $id = $moTree->id($leaf);
  $idToLeaf{$id} = $leaf;
  fail("Leaf named " . encode_entities($id) . " is not in the alignment")
    unless exists $alnSeq{$id};
}

# Issue a warning error for any sequences not in the tree
foreach my $id (keys %alnSeq) {
  warning("Sequence " . encode_entities($id) . " is not in the tree")
    unless exists $idToLeaf{$id};
}

my %branchLen = ();             # node to branch length
my $missingLen = 0; # set if any leaf or internal node had no branch length
# Convert any negative branch lengths to 0
# Convert any missing branch lengths to 1
foreach my $node (@$nodes) {
  next if $node == $root;
  my $len = $moTree->branch_length($node);
  if ($len eq "") {
    print p($node);
    $missingLen = 1;
    $len = 1;
  }
  $len = 0 if $len < 0;
  $branchLen{$node} = $len;
}
warning("Missing branch lengths were set to 1") if $missingLen;

my $anchorId = param('anchor');
$anchorId = "" if !defined $anchorId;
fail("Unknown anchor id " . $anchorId)
  if $anchorId ne "" && !exists $alnSeq{$anchorId};
my ($anchorAln, $anchorSeq, $anchorLen);
if ($anchorId ne "") {
  $anchorAln = $alnSeq{$anchorId};
  $anchorSeq = $anchorAln; $anchorSeq =~ s/[-]//g;
  $anchorLen = length($anchorSeq);
}

my @anchorPos;               # 1-based, and in the anchor if it is set
if (defined param('pos') && param('pos') ne "") {
  my $pos = param('pos');
  $pos =~ s/\s//g;
  @anchorPos = split /,/, $pos;
  foreach my $i (@anchorPos) {
    fail("Invalid position $i")
      unless $i =~ m/^\d+$/ && $i >= 1 && $i <= $alnLen;
    fail("position $i is past end of anchor " . encode_entities($anchorId))
      if $anchorId ne "" && $i > $anchorLen;
  }
}

my @alnPos = ();                # 0-based, and in the alignment
if ($anchorId eq "") {
  @alnPos = map { $_ - 1 } @anchorPos;
} else {
  my %anchorToAln = ();    # 1-based in anchor to 0-based in alignment
  my $at = 0;
  for (my $i = 0; $i < $alnLen; $i++) {
    my $c = substr($anchorAln, $i, 1);
    if ($c ne "-") {
      $at++;
      $anchorToAln{$at} = $i;
    }
  }
  @alnPos = map $anchorToAln{$_}, @anchorPos;
}

my @lines = ("Drawing the tree for", scalar(keys %idToLeaf), "proteins");
if (@alnPos > 0) {
  push @lines, ("with", scalar(@alnPos), "aligned positions");
  push @lines, ("numbered as in", encode_entities($anchorId))
    if $anchorId ne "";
}
print p(@lines);

# Build an svg
# Layout:
# y axis (0 at top):
# 1 blank row at top, of height $padTop
# 1 row per leaf, of height $rowHeight
# space for the scale bar at bottom, of height padBottom
# x axis: $padLeft to $padLeft + $treeWidth has the tree
# spacer of width $pdMiddle
# then 1 column for each position
my $padTop = 45;
my $rowHeight = 8;
my $padBottom = 25;
my $padLeft = 10;
my $treeWidth = 150;
my $padMiddle = 10;
my $posWidth = 30;
my $svgHeight = $padTop + scalar(@leaves) * $rowHeight + $padBottom;
my $svgWidth = $padLeft + $treeWidth + $padMiddle + scalar(@alnPos) * $posWidth;

my %rawY; # Unscaled y for each node (0 to nLeaves-1)
for (my $i = 0; $i < scalar(@leaves); $i++) {
  $rawY{ $leaves[$i] } = $i;
}
foreach my $node (reverse @$nodes) {
  if (!exists $rawY{$node}) {
    my @values = map $rawY{$_}, $moTree->children($node);
    die unless @values > 0;
    foreach my $value (@values) {
      die "rawY not set yet for child of $node" if !defined $value;
    }
    $rawY{$node} = sum(@values) / scalar(@values);
  }
}
my $maxY = max(values %rawY);
$maxY = 1 if $maxY == 0;
my %nodeY;
while (my ($node, $rawY) = each %rawY) {
  $nodeY{$node} = $padTop + $rowHeight * (0.5 + (scalar(@leaves)-1) * $rawY / $maxY);
}

my %rawX = ($root => 0); # Unscaled y, with root at 0
foreach my $node (@$nodes) {
  next if $node eq $root;
  my $parentX = $rawX{ $moTree->ancestor($node) };
  die $node unless defined $parentX;
  die $node unless defined $branchLen{$node};
  $rawX{$node} = $parentX + $branchLen{$node};
}
my %nodeX;
my $maxX = max(values %rawX);
while (my ($node, $rawX) = each %rawX) {
  $nodeX{$node} = $padLeft + $treeWidth * $rawX / $maxX;
}

my @svg = ();
foreach my $node (@$nodes) {
  if ($node != $root) {
    # draw line to ancestor
    my $parent = $moTree->ancestor($node);
    die unless defined $parent;
    push @svg, qq{<line x1="$nodeX{$node}" y1="$nodeY{$node}" x2="$nodeX{$parent}" y2="$nodeY{$parent}" stroke="black" />};
  }
  # draw node with popup info, if any
  my $id = $moTree->id($node);
  my $radius = 1.5;
  $radius = 3 if $moTree->is_Leaf($node) && $id eq $anchorId;
  push @svg, qq{<circle cx="$nodeX{$node}" cy="$nodeY{$node}" r="$radius">};
  my $title = $id;
  if ($moTree->is_Leaf($node) && @alnPos > 0) {
    my $longDesc = $id;
    $longDesc = "${id}: $alnDesc{$id}" if exists $alnDesc{$id};
    my $seq = $alnSeq{$id} || die;
    my @val = map substr($seq, $_, 1), @alnPos;
    $title = "$longDesc (has " . join("", @val) . ")";
  }
  $title = "(no label)" if $title eq "";
  push @svg, "<title>" . encode_entities($title) . "</title>"
    if defined $title  && $title ne "";
  push @svg, "</circle>";
}

# Scale bar
if (! $missingLen) {
  my @scales = reverse qw{0.001 0.002 0.005 0.01 0.02 0.05 0.1 0.2 0.5 1};
  while ($scales[0] > 0.8 * $maxX && @scales > 1) {
    shift @scales;
  }
  my $scaleSize = $scales[0];
  my $scaleLeft = $padLeft;
  my $scaleRight = $padLeft + $treeWidth * $scaleSize/$maxX;
  my $scaleY = $padTop + scalar(@leaves) * $rowHeight + $padBottom * 0.8;
  push @svg, qq{<line x1="$scaleLeft" y1="$scaleY" x2="$scaleRight" y2="$scaleY" stroke="black" />};
  my $scaleMid = ($scaleLeft+$scaleRight)/2;
  my $scaleY2 = $scaleY - 4;
  push @svg, qq{<text text-anchor="middle" x="$scaleMid" y="$scaleY2">$scaleSize /site</text>};
}

# The "taylor" color scheme is based on
# "Residual colours: a proposal for aminochromography" (Taylor 1997)
# and https://github.com/omarwagih/ggseqlogo/blob/master/R/col_schemes.r
# I added a dark-ish grey for gaps.
my %taylor = split /\s+/,
  qq{D #FF0000 S #FF3300 T #FF6600 G #FF9900 P #FFCC00 C #FFFF00 A #CCFF00 V #99FF00
     I #66FF00 L #33FF00 M #00FF00 F #00FF66 Y #00FFCC W #00CCFF H #0066FF R #0000FF
     K #6600FF N #CC00FF Q #FF00CC E #FF0066 - #555555};

# Show selected alignment positions (if any)
my $pos0X = $padLeft + $treeWidth + $padMiddle;
for (my $i = 0; $i < @alnPos; $i++) {
  my $left = $pos0X + $posWidth * $i;
  my $x = $left + $posWidth/2;
  my $labelY = $padTop - 8;
  my $labelChar = "#";
  $labelChar = substr($anchorAln, $alnPos[$i], 1) if $anchorId ne "";
  my $colLabel = $labelChar . $anchorPos[$i];
  push @svg, qq{<text text-anchor="left" transform="translate($x,$labelY) rotate(-45)">$colLabel</text>"};

  foreach my $leaf (@leaves) {
    my $id = $moTree->id($leaf);
    my $seq = $alnSeq{$id} || die;
    my $char = substr($seq, $alnPos[$i], 1);
    my $top = $nodeY{$leaf} - $rowHeight/2 - 0.5;
    my $heightUse = $rowHeight + 1; # slight overlap to avoid subtle white lines
    my $color = exists $taylor{$char} ? $taylor{$char} : "black";
    push @svg, qq{<rect x="$left" y="$top" width="$posWidth" height="$heightUse" style="fill:$color; stroke-width: 0;" />};
    # text-anchor centers horizontally, dominant-baseline centers vertically
    push @svg, qq{<text text-anchor="middle" dominant-baseline="middle" x="$x" y="$nodeY{$leaf}">$char</text>};
  }
}

print join("\n",
  "<DIV>",
  qq{<SVG width="$svgWidth" height="$svgHeight" style="position: relative; left: 1em;">},
  qq{<g transform scale(1)>},
  @svg,
  "</g>",
  "</SVG>",
  "</DIV>");

print end_html;
exit(0);

sub fail($) {
  my ($notice) = @_;
  print
    p(b($notice)),
      p(a({-href => "treeSites.cgi"}, "Try again")),
        end_html;
  exit(0);
}

sub warning($) {
  my ($notice) = @_;
  print p({-style => "color: red;"}, "Warning:", $notice), "\n";
}

#!/usr/bin/perl -w
# treeSites.cgi -- show key sites along a tree
use strict;
use CGI qw(:standard Vars start_ul);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use URI::Escape;
use HTML::Entities;
use Digest::MD5 qw{md5_hex};
use List::Util qw{sum min max};
use lib "../lib";
use pbutils qw{ReadFastaEntry ParseClustal};
use MOTree;

# In rendering mode, these options are required:
# aln or alnFile or alnMD5 -- alignment in multi-fasta format, or the file name (usually, md5 hash)
#   In header lines, anything after the initial space is assumed to be a description.
#   Either "." or "-" are gap characters.
# tree or treeFile or treeMD5 -- tree in newick format, or the file name (usually, md5 hash)
#   Labels on internal nodes are ignored.
#   Multi-furcations are allowed. It is treated as rooted.
# Optional arguments:
# anchor -- sequence to choose positiosn from
# pos -- comma-delimited list of positions (in anchor if set, or else, in alignment)
# tsvFile or tsvId -- descriptions for the ids (as uploaded file or file name)
# (This tool automatically saves uploads using alnId, treeId, or tsvId, and creates links using
#  the MD5 arguments.)
#
# If alignment or tree is missing, it shows an input form.

sub fail($); # print error to HTML and exit
sub warning($); # print warning to HTML
sub HandleTsvLines; # handle tab-delimited description lines

# maximum size of posted data, in bytes
my $maxMB = 25;
$CGI::POST_MAX = $maxMB*1024*1024;
my $maxN = 2000; # maximum number of sequences in alignment or leaves in tree

print
  header(-charset => 'utf-8'),
  start_html(-title => "Sites on a Tree",
             -script => [{ -type => "text/javascript", -src => "../static/treeSites.js"}]),
  h2("Sites on a Tree"),
  "\n";

my $alnSet = param('aln') || param('alnFile') || param('alnId');
my $treeSet = param('tree') || param('treeFile') || param('treeId');

my $tmpDir = "../tmp/aln";
if (!$alnSet || !$treeSet) {
  # show the form
  warning("Please give an alignment and a tree")
    if $alnSet xor $treeSet;
  print
    p("This tool will show a phylogenetic tree along with selected sites from a protein alignment."),
    start_form(-name => 'input', -method => 'POST', -action => 'treeSites.cgi'),
    p("Alignment in multi-fasta or clustal format (up to $maxN sequences or $maxMB megabytes):",
      br(),
      textarea(-name => 'aln', -value => '', -cols => 70, -rows => 10),
      br(),
      "Or upload:", filefield(-name => 'alnFile', -size => 50)),
    p("Rooted tree in newick format:",
      br(),
      textarea(-name => 'tree', -value => '', -cols => 70, -rows => 4),
      br(),
      "Or upload:", filefield(-name => 'treeFile', -size => 50)),
    p("Positions to show (optional):",
      textfield(-name => 'pos', -size => 50),
      br(), br(),
      "Positions are relative to sequence (optional):", textfield(-name => 'anchor', -size => 25)),
    p(submit()),
    end_form,
    p("Or try",
      a({-href => "sites.cgi" }, "SitesBLAST").":",
      "find homologs with known functional residues and see if they are conserved");
  print end_html;
  exit(0);
}

# else rendering mode
my @alnLines = ();
my $alnId;
if (param('aln')) {
  @alnLines = split '\n', param('aln');
} elsif (param('alnFile')) {
  my $fhAln = param('alnFile')->handle;
  fail("alnFile not a file") unless $fhAln;
  @alnLines = <$fhAln>;
} elsif (param('alnId')) {
  $alnId = param('alnId');
  fail("Invalid alnId") unless $alnId =~ m/^[a-zA-Z0-9]+$/;
  my $alnFile = "$tmpDir/$alnId.aln";
  open(my $fhAln, "<", $alnFile) || fail("Unknown alnId $alnId");
  @alnLines = <$fhAln>;
  close($fhAln) || die "Error reading $alnFile";
} else {
  fail("No alignment specified");
}

my %alnSeq;
my %alnDesc;

if (my $hash = ParseClustal(@alnLines)) {
  %alnSeq = %$hash;
} else {
  my $alnString = join("\n", @alnLines);
  open (my $fh, "<", \$alnString);
  my $state = {};
  while (my ($id, $seq) = ReadFastaEntry($fh, $state, 1)) {
    if ($id =~ m/^(\S+)\s(.*)/) {
      $id = $1;
      $alnDesc{$id} = $2;
    }
    fail("Duplicate sequence for " . encode_entities($id))
      if exists $alnSeq{$id};
    $seq =~ s/[.]/-/g;
    $seq = uc($seq);
    $seq =~ m/^[A-Z-]+$/ || fail("Invalid sequence for " . encode_entities($id));
    $alnSeq{$id} = $seq;
  }
  fail($state->{error}) if defined $state->{error};
}
fail("No sequences in the alignment")
  if (scalar(keys %alnSeq) == 0);

fail("Too many sequences in the alignment") if scalar(keys %alnSeq) > $maxN;

my $alnLen;
while (my ($id, $seq) = each %alnSeq) {
  $alnLen = length($seq) if !defined $alnLen;
  fail("Inconsistent sequence length for " . encode_entities($id))
    unless length($seq) == $alnLen;
}

my ($moTree, $treeId);

if (param('tree')) {
  eval { $moTree = MOTree::new('newick' => param('tree')) };
} elsif (param('treeFile')) {
  my $fh = param('treeFile')->handle;
  fail("treeFile not a file") unless $fh;
  eval { $moTree = MOTree::new('fh' => $fh) };
} elsif (param('treeId')) {
  $treeId = param('treeId');
  fail("Invalid treeId") unless $treeId =~ m/^[a-zA-Z0-9]+$/;
  my $file = "$tmpDir/$treeId.tree";
  open(my $fh, "<", $file) || fail ("Unknown treeId $treeId");
  eval { $moTree = MOTree::new('fh' => $fh) };
  close($fh) || die "Error reading $file";
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

# Issue a warning error for any sequences not in the tree, if this
# is the first time they were used together
if (!defined $alnId || !defined $treeId) {
  foreach my $id (keys %alnSeq) {
    warning("Sequence " . encode_entities($id) . " is not in the tree")
      unless exists $idToLeaf{$id};
  }
}

# Save the tree and alignment, if necessary
if (!defined $alnId) {
  $alnId = md5_hex(@alnLines);
  my $file = "$tmpDir/$alnId.aln";
  if (! -e $file) {
    open(my $fh, ">", $file) || die "Cannot write to $file";
    foreach my $line (@alnLines) {
      $line =~ s/[\r\n]+$//;
      print $fh $line."\n";
    }
    close($fh) || die "Error writing to $file";
  }
}

if (!defined $treeId) {
  my $newick = $moTree->toNewick();
  $treeId = md5_hex($newick);
  my $file = "$tmpDir/$treeId.tree";
  if (! -e $file) {
    open (my $fh, ">", $file) || die "Cannot write to $file";
    print $fh $newick."\n";
    close($fh) || die "Error writing to $file";
  }
}

# Try to parse the table to get descriptions
my $tsvId;
if (param('tsvFile')) {
  my $fh = param('tsvFile')->handle;
  fail("tsvFile is not a file") unless $fh;
  my @lines = <$fh>;
  my $n = HandleTsvLines(@lines);
  if ($n == 0) {
    warn("No descriptions found for matching ids in the uploaded table");
  } else {
    print p("Found $n descriptions in the uploaded table"),"\n";
    $tsvId = md5_hex(@lines);
    my $file = "$tmpDir/$tsvId.tsv";
    if (! -e $file) {
      open (my $fh, ">", $file) || die "Cannot write to $file";
      foreach my $line (@lines) {
        $line =~ s/[\r\n]+$//;
        print $fh $line."\n";
      }
      close($fh) || die "Error writing to $file";
    }
  }
} elsif (param('tsvId')) {
  $tsvId = param('tsvId');
  fail("Invalid tsvId") unless $tsvId =~ m/^[a-zA-Z0-9]+$/;
  my $file = "$tmpDir/$tsvId.tsv";
  open(my $fh, "<", $file) || fail("Unknown tsvId $tsvId");
  my @lines = <$fh>;
  close($fh) || die "Error reading $file";
  my $n = HandleTsvLines(@lines);
  warn("No descriptions found for matching ids in the table") if $n == 0;
}

# Finished loading input

my %branchLen = (); # node to branch length
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
fail("Unknown anchor id " . encode_entities($anchorId))
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

my @drawing = ("Drawing a tree for " . scalar(keys %idToLeaf) . " proteins.");
if (@alnPos > 0 && $anchorId ne "") {
  push @drawing, "Position numbering is from " . encode_entities($anchorId) . ".";
}
print p(@drawing);
print p("Download",
        a({ -href => "$tmpDir/$treeId.tree" }, "tree"),
        "or",
        a({ -href => "$tmpDir/$alnId.aln" }, "alignment").",",
        "or see",
        a({ -href => join("",
                          "treeSites.cgi?anchor=", uri_escape($anchorId),
                          "&pos=", join(",",@anchorPos),
                          "&treeId=", $treeId,
                          "&alnId=", $alnId,
                          "&tsvId=", $tsvId || "") },
          "permanent link"),
        "to this page, or",
        a({ -href => "treeSites.cgi" }, "upload new data").".");

print p(start_form(-method => 'GET', -action => 'treeSites.cgi'),
        hidden( -name => 'alnId', -default => $alnId, -override => 1),
        hidden( -name => 'treeId', -default => $treeId, -override => 1),
        hidden( -name => 'tsvId', -default => $tsvId, -override => 1),
        "Select positions",
        textfield(-name => "pos", -default => join(",",@anchorPos), -size => 30, -maxlength => 200),
        "in",
        textfield(-name => "anchor", -default => $anchorId, -size => 20, -maxlength => 200),
        submit(),
        end_form);

print p(start_form(-method => 'POST', -action => 'treeSites.cgi'),
        hidden( -name => 'alnId', -default => $alnId, -override => 1),
        hidden( -name => 'treeId', -default => $treeId, -override => 1),
        hidden( -name => 'anchor', -default => $anchorId, -override => 1),
        hidden( -name => 'pos', -default => join(",",@anchorPos), -override => 1),
        "Upload descriptions:", filefield(-name => 'tsvFile', -size => 50),
        submit(),
        br(),
        small("Descriptions should be tab-delimited with the sequence id in the 1st column and the description in the 2nd column"),
        end_form);

print p(start_form( -onsubmit => "return leafSearch();" ),
        "Search for sequences to highlight:",
        textfield(-name => 'query', -id => 'query', -size => 20),
        button(-name => 'Search', -onClick => "leafSearch()"),
        button(-name => 'Clear', -onClick => "leafClear()"),
        br(),
        span({-style => "font-size: 80%;", -id => "searchStatement"}, ""),
        end_form);

# Build an svg
# Layout:
# y axis (0 at top):
# 1 blank row at top, of height $padTop
# 1 row per leaf, of height $rowHeight
# space for the scale bar at bottom, of height padBottom
# x axis: $padLeft to $padLeft + $treeWidth has the tree
# spacer of width $pdMiddle
# then 1 column for each position
# and padRight
my $padTop = 45;
my $renderSmall = scalar(@leaves) > 100;
my $rowHeight = $renderSmall ? 3 : (scalar(@leaves) <= 20 ? 20 : 8);
my $minShowHeight = 20; # minimum height of a character to draw
my $padBottom = 45;
my $padLeft = 10;
my $treeWidth = 250;
my $padMiddle = 50;
my $padRight = 24;
my $posWidth = 30;
my $svgHeight = $padTop + scalar(@leaves) * $rowHeight + $padBottom;
my $svgWidth = $padLeft + $treeWidth + $padMiddle + scalar(@alnPos) * $posWidth + $padRight;

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

my %leafHas = (); # all the shown positions for that leaf
foreach my $leaf (@leaves) {
  my $id = $moTree->id($leaf);
  my $seq = $alnSeq{$id} || die;
  my @val = map substr($seq, $_, 1), @alnPos;
  $leafHas{$leaf} = join("", @val);
}

my %nodeTitle = ();
foreach my $node (@$nodes) {
  my $id = $moTree->id($node);
  my $title = encode_entities($id);
  if ($moTree->is_Leaf($node)) {
    $title .= ": " . encode_entities($alnDesc{$id}) if exists $alnDesc{$id};
    $title .= " (has $leafHas{$node})" if @alnPos > 0;
  }
  $nodeTitle{$node} = $title;
}

# The "taylor" color scheme is based on
# "Residual colours: a proposal for aminochromography" (Taylor 1997)
# and https://github.com/omarwagih/ggseqlogo/blob/master/R/col_schemes.r
# I added a dark-ish grey for gaps.
my %taylor = split /\s+/,
  qq{D #FF0000 S #FF3300 T #FF6600 G #FF9900 P #FFCC00 C #FFFF00 A #CCFF00 V #99FF00
     I #66FF00 L #33FF00 M #00FF00 F #00FF66 Y #00FFCC W #00CCFF H #0066FF R #0000FF
     K #6600FF N #CC00FF Q #FF00CC E #FF0066 - #555555};

my @svg = (); # lines in the svg

# For leaves, add an invisible horizontal bar with more opportunities for popup text
# These need to be output first to ensure they are behind everything else
for (my $i = 0; $i < @leaves; $i++) {
  my $leaf = $leaves[$i];
  my $x1 = 0; # $nodeX{$leaf} - 2;
  my $x2 = $svgWidth;
  my $width = $x2 - $x1;
  my $y1 = $nodeY{$leaf} - $rowHeight/2;
  push @svg, qq{<rect x="$x1" y="$y1" width="$x2" height="$rowHeight" fill="white" stroke="none" >};
  push @svg, "<TITLE>$nodeTitle{$leaf}</TITLE>\n</rect>";
}

# Show selected alignment positions (if any)
my $pos0X = $padLeft + $treeWidth + $padMiddle;
for (my $i = 0; $i < @alnPos; $i++) {
  my $pos = $alnPos[$i];
  my $left = $pos0X + $posWidth * $i;
  my $x = $left + $posWidth/2;
  my $labelY = $padTop - 3;
  my $labelChar = "#";
  $labelChar = substr($anchorAln, $pos, 1) if $anchorId ne "";
  my $colLabel = $labelChar . $anchorPos[$i];
  push @svg, qq{<text text-anchor="left" transform="translate($x,$labelY) rotate(-45)">$colLabel</text>};
  if (@leaves >= 20) {
    my $labelY2 = $svgHeight - $padBottom + 3;
    push @svg, qq{<text transform="translate($x,$labelY2) rotate(90)">$colLabel</text>"};
  }

  # draw boxes for every position
  foreach my $leaf (@leaves) {
    my $id = $moTree->id($leaf);
    my $seq = $alnSeq{$id} || die;
    my $char = uc(substr($seq, $pos, 1));
    my $top = $nodeY{$leaf} - $rowHeight/2;
    my $heightUse = $rowHeight + 0.2; # extra height to prevent thin white lines
    my $color = exists $taylor{$char} ? $taylor{$char} : "grey";
    my $encodedId = encode_entities($id);
    my $boxLeft = $left + $posWidth * 0.1;
    my $boxWidth = $posWidth * 0.8;
    push @svg, qq{<rect x="$boxLeft" y="$top" width="$boxWidth" height="$heightUse" style="fill:$color; stroke-width: 0;" >};
    push @svg, qq{<TITLE>$encodedId has $leafHas{$leaf}</TITLE>};
    push @svg, qq{</rect>};
  }

  # compute conservation of the position up the tree
  my %conservedAt = (); # node => value if it is conserved within this subtree, or ""
  foreach my $node (reverse @$nodes) {
    if ($moTree->is_Leaf($node)) {
      my $id = $moTree->id($node);
      die unless exists $alnSeq{$id};
      $conservedAt{$node} = substr($alnSeq{$id}, $pos, 1);
    } else {
      my @children = $moTree->children($node);
      my $char;
      foreach my $child (@children) {
        die unless exists $conservedAt{$child};
        if ($conservedAt{$child} eq "") {
          $char = "";
          last;
        } elsif (!defined $char) {
          $char = $conservedAt{$child};
        }  elsif ($char ne $conservedAt{$child}) {
          $char = "";
          last;
        }
      }
      $conservedAt{$node} = $char;
    }
  }

  # draw text for each clade at best location (if there is space)
  foreach my $node (@$nodes) {
    my $ancestor = $moTree->ancestor($node);
    if ($conservedAt{$node} ne "" && ($node == $root || $conservedAt{$ancestor} eq "")) {
      # Check if the height of this subtree is at least $minShowHeight
      my @leavesBelow;
      if ($moTree->is_Leaf($node)) {
        @leavesBelow = ($node);
      } else {
        @leavesBelow = @{ $moTree->all_leaves_below($node) };
      }
      # XXX show range of nodes, or number?? Or remove it?
      my @leavesBelowY = map $nodeY{$_}, @leavesBelow;
      my $height = $rowHeight + max(@leavesBelowY) - min(@leavesBelowY);
      next unless $height >= $minShowHeight;
      my $title = "";
      my $leafUse = $leavesBelow[0];
      my $id = encode_entities($moTree->id($leafUse));
      my $n1 = scalar(@leavesBelow) - 1;
      $title = ($n1 > 0 ? "$id and $n1 similar proteins have" : "$id has") . " " . $conservedAt{$leafUse}
        . " at " . $anchorPos[$i];
      push @svg, qq{<text text-anchor="middle" dominant-baseline="middle" x="$x" y="$nodeY{$node}"><TITLE>$title</TITLE>$conservedAt{$node}</text>};
    }
  }
} # End loop over positions

# Draw the tree after drawing the positions, so that text for leaf names (if displayed)
# goes on top of the color bars
foreach my $node (@$nodes) {
  my $radius = $renderSmall ? 2 : 3;
  my $style = "";
  if ($moTree->is_Leaf($node) && $moTree->id($node) eq $anchorId) {
    $radius = $renderSmall ? 2.5 : 4;
    $style = qq{fill="red"};
  }

  if ($node != $root) {
    # draw lines left and then up or down to ancestor
    my $parent = $moTree->ancestor($node);
    die unless defined $parent;
    push @svg, qq{<line x1="$nodeX{$node}" y1="$nodeY{$node}" x2="$nodeX{$parent}" y2="$nodeY{$node}" stroke="black" />};
    push @svg, qq{<line x1="$nodeX{$parent}" y1="$nodeY{$node}" x2="$nodeX{$parent}" y2="$nodeY{$parent}" stroke="black" />};
  }

  # draw node with popup info, if any
  # If it is a leaf, also make an (invisible) label; the group is to join these together
  push @svg, "<g>";
  push @svg, qq{<circle cx="$nodeX{$node}" cy="$nodeY{$node}" r="$radius" $style onclick="leafClick(this)">};
  push @svg, "<TITLE>$nodeTitle{$node}</TITLE>" if $nodeTitle{$node} ne "";
  push @svg, "</circle>";
  if ($moTree->is_Leaf($node)) {
    my $xLabel = $nodeX{$node} + $radius + 2;
    my $id = $moTree->id($node);
    push @svg, qq{<text dominant-baseline="middle" x="$xLabel" y="$nodeY{$node}" text-anchor="left" style="display:none; font-size:80%;">$id<TITLE>$nodeTitle{$node}</TITLE></text>};
  }
  push @svg, "</g>";
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

print join("\n",
  "<DIV>",
  qq{<SVG width="$svgWidth" height="$svgHeight" style="position: relative; left: 1em;">},
  qq{<g transform="scale(1)">},
  @svg,
  "</g>",
  "</SVG>",
  "</DIV>");

print end_html;
exit(0);

sub fail($) {
  my ($notice) = @_;
  my $URL = "treeSites.cgi";
  my $URL = "treeSites.cgi?alnId=$alnId&treeId=$treeId"
    if defined $alnId && defined $treeId;
  $URL .= "&tsvId=$tsvId" if defined $tsvId;
  print
    p(b($notice)),
      p(a({-href => $URL}, "Try again")),
        end_html;
  exit(0);
}

sub warning($) {
  my ($notice) = @_;
  print p({-style => "color: red;"}, "Warning:", $notice), "\n";
}

sub HandleTsvLines {
  my @lines = @_;
  my $n = 0;
  foreach my $line (@lines) {
    $line =~ s/[\r\n]+$//;
    my ($id, $desc) = split /\t/, $line;
    if (exists $alnSeq{$id} && defined $desc && $desc =~ m/\S/) {
      $alnDesc{$id} = $desc;
      $n++;
    }
  }
  return $n;
}


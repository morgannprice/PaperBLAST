#!/usr/bin/perl -w
# Given a UniProt identifier and a bl2seq alignment (with that identifier as the subject),
# make a tab-delimited table of features and whether they are conserved.

use strict;
use Bio::AlignIO;
use LWP::Simple;

my $debug = 0;
if (@ARGV > 0 && $ARGV[0] eq "-debug") {
  $debug = 1;
  shift @ARGV;
}
die "Usage: alnreport.pl [ -debug ] subject bl2seq.out\n" unless @ARGV==2;
my ($subjectId, $blfile) = @ARGV;

my $str = Bio::AlignIO->new(-file => $blfile, -format => 'bl2seq');
my $aln = $str->next_aln();
exit(0) if !defined $aln; # no hits in input file
$aln->num_sequences() == 2 || die "Do not have two esquences in $blfile";
my $numpos = $aln->length();
my $qseq = $aln->get_seq_by_pos(1); # 1st sequence in alignment (1-based!)
my $sseq = $aln->get_seq_by_pos(2);

my %sposToQpos = ();
my %sposToAlnpos = ();
foreach my $i (1..$numpos) {
  my $spos = $sseq->location_from_column($i);
  if (defined $spos && $spos->location_type() eq 'EXACT') {
    $sposToAlnpos{ $spos->start() } = $i;
    my $qpos = $qseq->location_from_column($i);
    if (defined $qpos && $qpos->location_type() eq 'EXACT') {
      $sposToQpos{ $spos->start() } = $qpos->start();
    }
  }
}

# Fetch the UniProt entry
my $url = "https://www.uniprot.org/uniprot/$subjectId.txt";
print STDERR "Fetching $url\n" if $debug;
my $uniprot = LWP::Simple::get($url);
chomp $uniprot;
die "Cannot fetch $url" unless $uniprot;
my @lines = split /\n/, $uniprot;
print STDERR "Received " . scalar(@lines) . " lines\n" if $debug;
die "Cannot fetch $url" unless $lines[-1] eq "//";
my @ft = grep m/^FT /, @lines;
my @ft2 = (); # combining the continuation lines
foreach my $ft (@ft) {
  if ($ft =~ m/^FT   [A-Z]/) { # 3 spaces
    push @ft2, $ft;
  } elsif ($ft =~ m/^FT    +(.*)/) { # at least 4 spaces
    die "Unexpected continuation FT line in the UniProt entry for $subjectId"
      if @ft2 == 0;
    $ft2[-1] .= " $1";
  } else {
    die "Unexpected FT line\n$ft\nin the UniProt entry for $subjectId";
  }
}

my %ft = (); # type => list of [ begin, end, comment ]
foreach my $ft (@ft2) {
  $ft =~ m/^FT +([A-Z_]+) +([0-9<>]+) +([0-9<>]+) *(.*)$/
    || die "Invalid FT\n$ft\nin UniProt entry for $subjectId";
  my ($type,$begin,$end,$comment) = ($1,$2,$3,$4);
  # remove evidence codes from comments. The top two are curated; the remainder are not;
  # but decided not to distinguish them.
  $comment .= " (by similarity)"
    if $comment =~ m/ECO:0000250/ || $comment =~ m/ECO:0000255/
      || $comment =~ m/ECO:0000256/ || $comment =~ m/ECO:0000259/
        || $comment =~ m/ECO:0000259/ || $comment =~ m/ECO:0000213/;
  $comment =~ s/[{]ECO:\d+[^}]+[}][.]?//g;
  $comment =~ s/ +/ /g;
  push @{ $ft{$type} }, [ $begin, $end, $comment ]
    if $begin =~ m/^\d+$/ && $end =~ m/^\d+$/;
}

# Which types of features to report on, in what order
my @types = qw/ACT_SITE BINDING CA_BIND ZN_FING METAL DNA_BIND NP_BIND MUTAGEN DOMAIN MOTIF REGION REPEAT PROPEP SIGNAL TRANSIT TRANSMEM INTRAMEM NON_STD MOD_RES CARBOHYD LIPID DISULFID CROSSLNK HELIX STRAND COILED TURN SITE/;

print join("\t", qw/type sbegin send qbegin qend nAligned nMatch subjectSeq querySeq comment/)."\n";

foreach my $type (@types) {
  next unless exists $ft{$type};
  foreach my $o (@{ $ft{$type} }) {
    my ($begin, $end, $comment) = @$o;
    # What fraction are in the alignment?
    my $nPos = $end-$begin+1;
    my $nAln = "";
    my $nMatch = "";
    my @qpos = ();
    my ($begAln, $endAln);
    my @pos = ($begin..$end);
    @pos = ($begin,$end) if ($type eq "DISULFID");
      


    foreach my $i (@pos) {
      if (exists $sposToQpos{$i}) {
        $nAln++;
        my $alnpos = $sposToAlnpos{$i};
        $begAln = $alnpos if !defined $begAln;
        $endAln = $alnpos;
        my $sval = substr($sseq->seq, $alnpos-1, 1);
        my $qval = substr($qseq->seq, $alnpos-1, 1);
        $nMatch++ if $sval eq $qval;
        push @qpos, $sposToQpos{$i};
      }
    }

    my $qbeg = @qpos > 0 ? $qpos[0] : "";
    my $qend = @qpos > 0 ? $qpos[-1] : "";
    my $partS = @qpos > 0 ? substr($sseq->seq, $begAln-1, $endAln-$begAln+1) : "";
    my $partQ = @qpos > 0 ? substr($qseq->seq, $begAln-1, $endAln-$begAln+1) : "";
    if ($type eq "DISULFID" && @pos > 1 && $partS ne "") {
      $partS = substr($partS, 0, 1) . "," . substr($partS, length($partS)-1, 1);
      $partQ = substr($partQ, 0, 1) . "," . substr($partQ, length($partQ)-1, 1)
        if $partQ ne "";
    }
    print join("\t", $type, $begin, $end, $qbeg, $qend, $nAln, $nMatch, $partS, $partQ, $comment)."\n";
  }
}


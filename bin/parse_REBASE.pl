#!/usr/bin/perl -w
use strict;

my $usage = <<END
Usage: parse_REBASE.pl rebase_directory > REBASE.curated_parsed

The rebase directory must include
bairoch.txt (in swissprot-like format)
protein_seqs.txt (in fasta-like format)
END
;
die $usage unless @ARGV == 1;
my ($dir) = $ARGV[0];
die $usage unless -d $dir;

die "No such file: $dir/bairoch.txt\n" unless -e "$dir/bairoch.txt";
die "No such file: $dir/protein_seqs.txt\n" unless -e "$dir/protein_seqs.txt";

my %entries = (); # ID => hash with various attributes
$/ = "\/\/\n";
open(BAIROCH, "<", "$dir/bairoch.txt") || die "Cannot read $dir/bairoch.txt";
my $nDupEntry = 0;
while(my $record = <BAIROCH>) {
  my @lines = split /\n/, $record;
  my $entry = {};
  foreach my $line (@lines) {
    chomp $line;
    if ($line =~ m/^([A-Z][A-Z])   (.*)$/) {
      my ($key,$value) = ($1,$2);
      next if $key eq "CC"; #ignore comments
      if ($key eq "RN" || $key eq "RA" || $key eq "RL") {
        # allow duplicates;
      } else {
        die "Duplicate value for $key" if exists $entry->{$key};
      }
      $entry->{$key} = $value;
    } else {
      die "Unexpected line $line" unless $line eq "" || $line eq "//";
    }
  }
  if (exists $entry->{ID}) {
    my $id = $entry->{ID};
    if (exists $entries{$id}) {
      $nDupEntry++;
    } else {
      $entries{$id} = $entry;
    }
  } else {
    # allow empty records
    die "No ID in entry" if keys %$entry > 0;
  }
}
close(BAIROCH) || die "Error reading $dir/bairoch.txt";

print STDERR "Read " . scalar(keys %entries) . " entries (ignoring $nDupEntry duplicates) from $dir/bairoch.txt\n";

my %printed = ();
open(SEQS, "<", "$dir/protein_seqs.txt") || die "Cannot read $dir/protein_seqs.txt";
$/ = "<>";
while (my $record = <SEQS>) {
  my @lines = split /\n/, $record;
  next if @lines == 0;
  chomp $lines[-1];
  $lines[-1] =~ s/<>$//;
  @lines = grep { $_ ne "" } @lines;
  shift @lines if $lines[0] eq "REBASE protein sequences";
  my $header = shift @lines;
  die "Invalid header $header" unless $header =~ m/^>(.*)$/;
  $header = $1;
  foreach my $line (@lines) {
    die "Invalid sequence line $line" unless $line =~ m/^[A-Z ]+$/;
  }
  my $seq = join("", @lines);
  $seq =~ s/ //g;

  # Parse the header into a list of key:value pairs
  my %attr = ();
  my @pieces = split /\t/, $header;
  foreach my $piece (@pieces) {
    $piece =~ m/^([a-zA-Z]+):(.*)$/ || die "Invalid header piece $piece";
    die "Strange header $header for $1" if exists $attr{$1};
    $attr{$1} = $2;
  }
  my $id = $attr{REBASE};
  if (exists $entries{$id} && ! $printed{$id}) {
    my $entry = $entries{$id};
    die "No SeqLength for $id" unless exists $attr{SeqLength};
    die "Wrong SeqLength for $id" unless $attr{SeqLength} == length($seq);
    # The curated_parsed fields are:
    # database, identifier, id2, short name, description, organism, sequence, comment, pmids
    die "No EnzType for $id" unless exists $attr{EnzType};
    my $comment = ""; # should include more information such as references
    die "No RS entry" unless $entry->{RS};
    next unless $entry->{RS} =~ m/[ACGT]/; # must have some sequence specificity
    my $desc = $attr{EnzType} . " recognizing " . $entry->{RS};
    $desc =~ s/, [?];$//;
    $desc .= " with methylation at " . $entry->{MS} if $entry->{MS};
    $desc =~ s/;$//;
    print join("\t", "REBASE", $id, $attr{GenBank} || "", "",
               $desc,
               $entry->{OS}, $seq, $comment, "") . "\n";
    $printed{$id} = 1;
  }
}
close(SEQS) || die "Error reading $dir/protein_seqs.txt";

# Utilities for PaperBLAST
package pbutils;
require Exporter;
use strict;
use File::stat;

our (@ISA,@EXPORT);
@ISA = qw(Exporter);
@EXPORT = qw(read_list wget ftp_html_to_files write_list mkdir_if_needed ReadFasta ParsePTools ReadFastaEntry ReadTable ReadColumnNames NewerThan);

sub read_list($) {
  my ($file) = @_;
  open(LIST, "<", $file) || die "Cannot read $file";
    my @lines = <LIST>;
    close(LIST) || die "Error reading $file";
    @lines = map { chomp; $_ } @lines;
    return(@lines);
}

# Retry wget calls to increase odds of success
my $nFailures = 0;
my $nFailMax = 200;

sub wget($$) {
    my ($url, $file) = @_;
    if ($ENV{PB_DOWNLOAD_PASS} && -s $file) {
      print STDERR "Using existing non-empty file for $file\n";
    } else {
      while(system("wget", "-nv", "-O", $file, $url) != 0) {
        $nFailures++;
        die "Failed to load $url\n" if $nFailures >= $nFailMax;
        print STDERR "Waiting and retrying\n";
        sleep(2);
      }
    }
}

sub ftp_html_to_files($) {
    my ($listfile) = @_;
    my @files = ();
    open(IN, "<", $listfile) || die;
    while(<IN>) {
        chomp;
        next unless m!>([A-Za-z0-9_.-]+)<!;
        my $file = $1;
        push @files, $file;
    }
    close(IN) || die "Error reading $listfile";
    return @files;
}

sub write_list($$) {
    my ($files, $outfile) = @_;
    open(LIST, ">", $outfile) || die "Cannot write to $outfile";
    foreach my $file (@$files) {
        print LIST "$file\n";
    }
    close(LIST) || die "Error writing to $outfile";
    print STDERR "Wrote $outfile\n";

}

sub mkdir_if_needed($) {
    my ($subdir) = @_;
    (-d $subdir) || mkdir($subdir) || die "Cannot mkdir $subdir";
}

# returns a reference to a hash of name => sequence
sub ReadFasta ($) {
    my ($filename) = @_;
    open(IN, "<", $filename) || die "Cannot read $filename";
    my %seqs = ();
    my $name = undef;
    while(<IN>) {
	chomp;
	if (m/>(\S+)/) {
	    die "Duplicate entry for $1 in filename" if exists $seqs{$1};
	    $name = $1;
	} else {
	    die "Unexpected line $_" unless defined $name;
	    $seqs{$name} .= $_;
	}
    }
    close(IN) || die "Error reading $filename";
    return(\%seqs);
}

# Read a record from a pathway tools attribute-value file
# and returns the result as a reference to a hash of
# attribute_name => list of { value => value, annotation_name => anno }
# where the key 'value' is the main attribute value and the annotation names vary.
#
# Multiple values for an attribute are supported, but if there are
# multiple annotations with the same name for an attribute/value pair,
# then only the first one is reported.
#
# Parsing is based on documentation from
# https://bioinformatics.ai.sri.com/ptools/flatfile-format.html
sub ParsePTools {
  my ($fh) = shift;
  die unless defined $fh;
  local($/) = "\n//\n";
  my $record = <$fh>;
  return undef if !defined $record;
  my @lines = split /\n/, $record;
  foreach (@lines) { chomp; }
  # Remove comments
  @lines = grep { ! m/^#/ } @lines;
  return undef if @lines == 0;
  $lines[-1] eq "//" || die "Last lines in record does not match //";
  pop @lines;
  # Merge trailing lines
  my @merged = ();
  foreach my $line (@lines) {
    if ($line =~ m!^/!) {
      die "Continuation line with no preceding line: $line"
        unless @merged > 0;
      $merged[-1] .= " " . substr($line, 1);
    } else {
      push @merged, $line;
    }
  }
  my $last_attr = undef;
  my $out = {};
  foreach my $line (@merged) {
    $line =~ m/^(\S+) - (.*)$/ || die "Cannot parse attribute or annotation from $line";
    my ($attr, $value) = ($1,$2);
    if ($attr =~ m/^\^/) {
      my $anno = substr($attr, 1);
      die "Annotation with no preceding attribute: $line" unless defined $last_attr;
      my $h = $out->{$last_attr}[-1];
      if (exists $h->{$anno}) {
        # print STDERR "Duplicate annotation for $last_attr $anno\n";
        ;
      } else {
        $h->{$anno} = $value;
      }
    } else {
      push @{ $out->{$attr} }, { 'value' => $value };
      $last_attr = $attr;
    }
  }
  return $out;
}

# Read one entry at a time from a fasta file
# The first argument is a hash to keep track of saved state, i.e.:
#   my $state = {};
#   while(my ($header,$sequence) = ReadFastaEntry($fh,$state)) { ... }
# (header will have the ">" removed)
sub ReadFastaEntry {
  my ($fh, $state) = @_;
  die unless ref $state;
  return () if exists $state->{DONE}; # end state
  # initialization
  if (!defined $state->{header}) {
    $state->{header} = "";
    $state->{sequence} = "";
  }
  while (my $line = <$fh>) {
    chomp $line;
    if ($line =~ m/^>(.*)/) {
      my $old_header = $state->{"header"};
      my $old_sequence = $state->{"sequence"};
      $state->{"header"} = $1;
      die "Empty header in $line" if $state->{header} eq "";
      $state->{"sequence"} = "";
      return ($old_header, $old_sequence) if $old_header ne "";
    } else {
      die "Unexpected sequence with no header" if $state->{"header"} eq "";
      $line =~ s/ //g;
      $line = uc($line);
      # allow - or . as used in alignments and * as used for stop codons
      die "Invalid sequence line $line" unless $line =~ m/^[A-Z*.-]*$/;
      $state->{sequence} .= $line;
    }
  }
  # reached end of file
  $state->{DONE} = 1;
  return () if $state->{header} eq ""; # empty file
  return ($state->{header}, $state->{sequence});
}

# filename and list of required fields => list of hashes, each with field->value
# The list can be a single name or a reference to a list
sub ReadTable {
    my ($filename,@required) = @_;
    if (scalar(@required) == 1 && ref $required[0]) {
        @required = @{ $required[0] };
    }
    open(IN, "<", $filename) || die "Cannot read $filename";
    my $headerLine = <IN>;
    $headerLine =~ s/[\r\n]+$//; # for DOS
    # Check for Mac style files -- these are not supported, but give a useful error
    die "Tab-delimited input file $filename is a Mac-style text file, which is not supported\n"
        . "Use\ndos2unix -c mac $filename\nto convert it to a Unix text file.\n"
        if $headerLine =~ m/\t/ && $headerLine =~ m/\r/;
    my @cols = split /\t/, $headerLine;
    my %cols = map { $cols[$_] => $_ } (0..(scalar(@cols)-1));
    foreach my $field (@required) {
	die "No field $field in $filename" unless exists $cols{$field};
    }
    my @rows = ();
    while(my $line = <IN>) {
	$line =~ s/[\r\n]+$//;
	my @F = split /\t/, $line, -1;
	die "Wrong number of columns in:\n$line\nin $filename"
	    unless scalar(@F) == scalar(@cols);
	my %row = map { $cols[$_] => $F[$_] } (0..(scalar(@cols)-1));
	push @rows, \%row;
    }
    close(IN) || die "Error reading $filename";
    return @rows;
}

# filename to list of column names
sub ReadColumnNames($) {
    my ($filename) = @_;
    open(IN, "<", $filename) || die "Cannot read $filename";
    my $line = <IN>;
    close(IN) || die "Error reading $filename";

    $line =~ s/[\r\n]+$//; # for DOS
    my @cols = split /\t/, $line;
    return @cols;
}

sub NewerThan($$) {
    my ($file1, $file2) = @_;
    die "Invalid arguments to NewerThan" unless defined $file1 && defined $file2;
    die "No such file: $file2" unless -e $file2;
    return 0 unless -e $file1;
    return stat($file1)->mtime >= stat($file2)->mtime ? 1 : 0;
}

1;

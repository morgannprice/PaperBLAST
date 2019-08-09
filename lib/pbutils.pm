# Utilities for PaperBLAST
package pbutils;
require Exporter;
use strict;
use File::stat;

our (@ISA,@EXPORT);
@ISA = qw(Exporter);
@EXPORT = qw!read_list wget ftp_html_to_files write_list mkdir_if_needed
             ReadFasta ReadFastaDesc ReadFastaEntry
             ParsePTools
             ReadTable ReadColumnNames
             NewerThan
             CuratedMatch CuratedWordMatch
             IdToUniqId FetchSeqs UniqIdToSeq
             FetchCuratedInfo!;

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

# Returns a hash containing either "error"
# or hashes of "seq", "desc", and "len"
sub ReadFastaDesc($) {
    my ($file) = @_;
    my %seq = ();
    my %desc = ();
    my $name = undef;
    open(FAA, "<", $file) || return('error' => "Cannot read $file" );
    while(<FAA>) {
        s/[\r\n]+$//;
        if (m/^>(.*)$/) {
            my $header = $1;
            if ($header =~ m/^(\S+)\s+(\S.*)$/) {
                $name = $1;
                $desc{$name} = $2;
            } else {
                return('error' => "bad header for sequence:\n$header\n") unless $header =~ m/^\S+$/;
                $name = $header;
                $desc{$name} = $header;
            }
            return('error' => "Duplicate sequence id:\n$name\n") if exists $seq{$name};
            $seq{$name} = "";
        } else {
            return('error' => "sequence before header:\n$_\n") unless defined $name;
            s/\s//g;
            $seq{$name} .= $_;
        }
    }
    close(FAA) || return('error' => "Error reading $file");
    my %len = ();
    while (my ($name,$seq) = each %seq) {
        $len{$name} = length($seq);
        return('error' => "No sequence for id:\n$name\n") if ($len{$name} == 0);
    }
    return("desc" => \%desc, "seq" => \%seq, "len" => \%len);
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
#
# If using $return_error, then on an error it sets $state->{error} and returns 0

sub ReadFastaEntry {
  my ($fh, $state, $return_error) = @_;
  die unless ref $state;
  return () if exists $state->{DONE}; # end state
  return () if exists $state->{error};
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
      if ($state->{header} eq "") {
        $state->{error} = "Empty header in $line";
        return () if $return_error;
        die $state->{error};
      }
      $state->{"sequence"} = "";
      return ($old_header, $old_sequence) if $old_header ne "";
    } else {
      if ($state->{"header"} eq "") {
        $state->{error} = "Unexpected sequence with no header" ;
        return () if $return_error;
        die $state->{error};
      }
      $line = uc($line);
      # allow - or . as used in alignments and * as used for stop codons
      unless ($line =~ m/^[A-Z*.-]*$/) {
        $state->{error} = "Invalid sequence line $line";
        return () if $return_error;
        die $state->{error};
      }
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

# Returns reference to a list of rows from CuratedGene
sub CuratedMatch($$$) {
  my ($dbh, $query, $limit) = @_;
  return $dbh->selectall_arrayref("SELECT * FROM CuratedGene WHERE desc LIKE ? LIMIT ?",
                                     { Slice => {} }, "%" . $query . "%", $limit);
}

# Given rows that include the desc field, return a reference to a
# list with the rows whose desc matches the query as words
# % is the wild card; all other characters are treated literally
# The optional third argument is which field to use (defaults to "desc")
sub CuratedWordMatch {
  my ($chits, $query, $field) = @_;
  die "Invalid chits argument to CuratedWordMatch"
    unless ref $chits eq "ARRAY";
  die "Invalid query argument to CuratedWordMatch"
    unless ref $query eq "" && defined $query;
  $field = "desc" if !defined $field;
  die "Invalid field argument to CuratedWordMatch"
    unless ref $field eq "";

  # bracket characters are trouble, so remove from both query and tested description
  my $sub = "[" . quotemeta("[]{}") . "]";
  $query =~ s/$sub//g;
  my $quoted = quotemeta($query); # this will quote % as well
  # turn % into a separation of words; note quoting of \\ so that \b it appears in the new string
  $quoted =~ s/\\%/\\b.*\\b/g;
  my @keep = grep { my $test = $_->{$field};
                    die "Missing desc field in CuratedWordMatch()" unless defined $test;
                    $test =~ s/$sub//g;
                    $test =~ m/\b$quoted\b/i } @$chits;
  return \@keep;
}

# Replace a sequence id with its uniq id, if it is a known duplicate in the SeqToDuplicate table
sub IdToUniqId($$) {
  my ($dbh, $seqid) = @_;
  my $uniqRef = $dbh->selectcol_arrayref("SELECT sequence_id FROM SeqToDuplicate WHERE duplicate_id = ? LIMIT 1",
                                         {}, $seqid);
  return $uniqRef->[0] if @$uniqRef;
  return $seqid;
}

sub FetchSeqs($$$$) {
  my ($blastdir, $blastdb, $uniqids, $outfile) = @_;
  if (@$uniqids==0) {
    # fastacmd fails if given an empty list
    open(my $fh, ">", $outfile) || die "Cannot write to $outfile";
    close($fh) || die "Error writing to $outfile";
    return;
  }
  my $fastacmd = "$blastdir/fastacmd";
  die "No such executable: $fastacmd" unless -x $fastacmd;
  my $tmpdir = $ENV{TMP} || "/tmp";
  die "Not a directory: $tmpdir" unless -d $tmpdir;
  my $listFile = "$tmpdir/list.$$.in";
  open(my $fh, ">", $listFile) || die "Cannot write to $listFile";
  foreach my $uniqid (@$uniqids) {
    print $fh "$uniqid\n";
  }
  close($fh) || die "Error writing to $listFile";
  my $cmd = "$fastacmd -i $listFile -d $blastdb -p T > $outfile";
  system($cmd) == 0
    || die "$cmd failed: $!";
  unlink($listFile);
}

sub UniqIdToSeq($$$) {
  my ($pbdir, $blastdir, $uniqId) = @_;
  my $tmpfile = "/tmp/$$.uniqtoseq.faa";
  FetchSeqs($blastdir, "$pbdir/uniq.faa", [$uniqId], $tmpfile);
  my $seqs = ReadFasta($tmpfile);
  unlink($tmpfile);
  my @v = values(%$seqs);
  die unless @v == 1;
  return $v[0];
}

# PaperBLAST database handle and uniqId => reference to a list of (db,protId,desc) pairs
sub FetchCuratedInfo($$) {
  my ($pb, $uniqId) = @_;
  my $dups = $pb->selectcol_arrayref(qq{ SELECT duplicate_id FROM SeqToDuplicate
                                         WHERE sequence_id = ? },
                                     {}, $uniqId);
  my @ids = ( $uniqId );
  push @ids, @$dups;
  my @out = ();
  foreach my $id (@ids) {
    if ($id =~ m/^([^:]+)::(.*)$/) {
      my ($db,$protId) = ($1,$2);
      my ($desc) = $pb->selectrow_array("SELECT desc FROM CuratedGene WHERE db = ? AND protId = ?",
                                        {}, $db, $protId);
      die "No CuratedGene entry for ${db}::${protId}" unless defined $desc;
      push @out, [$db,$protId, $desc];
    }
  }
  @out = sort { $a->[0] cmp $b->[0] || $a->[1] cmp $b->[1] } @out;
  return \@out;
}

1;

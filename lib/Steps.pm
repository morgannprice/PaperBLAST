package Steps;
require Exporter;
use strict;
our (@ISA,@EXPORT);
@ISA = qw(Exporter);
@EXPORT = qw(ReadSteps FetchUniProtSequence AssemblyToProt
             WriteAssemblyAsOrgProteins WriteSixFrameAsOrgProteins
             ReadOrgTable ReadOrgProtein ParseOrgLocus);
use LWP::Simple qw{get};
use LWP::UserAgent;

sub ReadSteps($); # stepfile => hash of steps, rules, ruleOrder
# steps is a hash of step name to a hash of i, name, desc, and search, which is a list of type/value pairs
#	(use i to sort them by the same order as in the stepfile)
# rules is a hash of rule name to a list of lists, each element being a step name or another rule name
#	The requirement is met if any of the sublists are met (OR at the top level)
#	The sublist is met if all of its components are met (AND at the 2nd level)
# ruleOrder is a list of rule names in order
sub ReadSteps($) {
  my ($stepsfile) = @_;
  my  @stepTypes = qw{term hmm EC uniprot ignore ignore_other curated};
  my %stepTypes = map { $_ => 1 } @stepTypes;
  my $steps = {}; # step name to hash of name, desc, and search, which is a list of pairs
  my $rules = {}; # rule name to list of lists, each element is a step name or another rule name
  my @ruleOrder = (); # list of rules in input order
  open(my $fhSteps, "<", $stepsfile) || die "Cannot read $stepsfile";
  my $nSteps = 0;
  while (my $line = <$fhSteps>) {
    $line =~ s/[\r\n]+$//;
    $line =~ s/#.*//;
    $line =~ s/\s+$//;
    next if $line eq ""; # empty or comment line
    if ($line =~ m/^(\S+):\s+(\S.*)$/) {
      my ($rulename, $pieces) = ($1, $2);
      die "Cannot use step name $rulename as rule name in $stepsfile\n"
        if exists $steps->{$rulename};
      my @pieces = split /\s+/, $pieces;
      foreach my $piece (@pieces) {
        die "Rule $rulename includes $piece which has not been defined yet in $stepsfile\n"
          unless exists $rules->{$piece} || exists $steps->{$piece};
        die "Rule $rulename cannot include itself in $stepsfile\n"
          if $piece eq $rulename;
      }
      push @ruleOrder, $rulename unless exists $rules->{$rulename};
      push @{ $rules->{$rulename} }, \@pieces;
    } elsif ($line =~ m/\t/) {
      my @F = split /\t/, $line;
      die "Invalid rule line: $line\n" unless @F >= 3;
      my $stepname = shift @F;
      die "Duplicate definition of step $stepname in $stepsfile\n"
        if exists $steps->{$stepname};
      die "Spaces or other whitespace within step name $stepname in $stepsfile are not allowed\n"
        if $stepname =~ m/\s/;
      my $desc = shift @F;
      my @search = ();
      foreach my $search (@F) {
        $search =~ m/^([a-zA-Z_]+):(.*)$/
          || die "Do not recognize search specifier $search for step $stepname in $stepsfile\n";
        my ($type, $value) = ($1,$2);
        die "Do not recognize the type in search specifier $search for step $stepname in $stepsfile\n"
          unless exists $stepTypes{$type};
        push @search, [$type, $value];
      }
      die "Invalid step has no search items: $line\n"
        unless @search > 0;
      $steps->{$stepname} = { 'name' => $stepname, 'desc' => $desc, 'search' => \@search, 'i' => $nSteps++ };
    } else {
      die "Do not recognize line $line in $stepsfile as either rule or step\n";
    }
  }
  close($fhSteps) || die "Error reading $stepsfile";

  # Check for cyclic dependencies within the rules
  # Specifically, when going through rules in order, should
  # only depend on rules that have already been defined.
  my %sofar = ();
  foreach my $rule (@ruleOrder) {
    foreach my $reqs (@{ $rules->{$rule} }) {
      foreach my $piece (@$reqs) {
        if (!exists $steps->{$piece} && !exists $sofar{$piece}) {
          die "Out-of-order dependency error in ${stepsfile}\n"
            . "A rule for $rule depends on $piece, but $piece is defined later\n";
        }
      }
    }
    $sofar{$rule} = 1;
  }
  return { 'steps' => $steps, 'rules' => $rules, 'ruleOrder' => \@ruleOrder };
}

# Returns sequence and description (just joining the DE lines)
sub FetchUniProtSequence($) {
  my ($id) = @_;
  # fetch the full entry because this allows either ids or accessions (i.e., Q72CU7_DESVH or Q72CU7) to be used
  my $URL = "https://www.uniprot.org/uniprot/${id}.txt";

#  my $ua = LWP::UserAgent->new( ssl_opts => { SSL_verify_mode => 'SSL_VERIFY_NONE' } );
#  my $response = $ua->get($URL);
#  die "redirect" if $response->is_redirect;
#  die "Cannot contact $URL: " . $response->as_string unless $response->is_success;
#  my $content = $response->decoded_content;
#  return undef unless $content;
  my $content = get($URL);
  die "Failed to fetch $URL\n" unless $content;
  my @lines = split /\n/, $content;

  my @de = ();
  my $seqlen = undef;
  while (@lines > 0) {
    my $line = shift @lines;
    if ($line =~ m/^DE +(.*)$/) {
      push @de, $1;
    }
    if ($line =~ m/^SQ +SEQUENCE +(\d+)/) {
      $seqlen = $1;
      last;
    }
  }
  die "Invalid uniprot entry $URL -- no SQ line found\n"
    unless defined $seqlen;
  my $seq = "";
  while (@lines > 0) {
    my $line = shift @lines;
    if ($line =~ m!^//!) {
      last;
    } else {
      die "Invalid sequence line $line\n"
        unless $line =~ m/^[A-Z ]+$/;
      $line =~ s/ //g;
      $seq .= $line;
    }
  }
  die "Error processing sequence from $URL -- sequence length mismatch\n"
    unless length($seq) == $seqlen;
  return ($seq, join(" ", @de));
}

# Given a fetched assembly, such as from FetchAssembly::CacheAssembly,
# and a filehandle, write out the predicted proteins as a fasta file with
# header lines >orgId:locusId sysName desc
# (orgId = gdb__genomeid will be used as the organism identifier)
# to the file handle
# Returns the orgId
sub WriteAssemblyAsOrgProteins($$) {
  my ($assembly, $fhOut) = @_;
  die "Invalid assembly" unless $assembly->{gdb} && $assembly->{gid};
  die ": or whitespace not allowed in assembly identifier $assembly->{gid}"
    if $assembly->{gid} =~ m/\s/ || $assembly->{gid} =~ m/:/;
  my $orgId = join("__", $assembly->{gdb}, $assembly->{gid});
  my $faaFile = $assembly->{faafile} || die "No faa file";
  open(my $fhIn, "<", $faaFile) || die "Cannot read $faaFile";
  my $state = {};
  my %seqSeen = ();
  while (my ($header, $seq) = ReadFastaEntry($fhIn, $state)) {
    die "Duplicate sequence for $header in $faaFile"
      if exists $seqSeen{$header};
    $seqSeen{$header} = 1;
    my ($locusId, $sysName, $desc) = ("", "", "");
    my @words = split / /, $header, -1;
    if ($assembly->{gdb} eq "FitnessBrowser" || $assembly->{gdb} eq "MicrobesOnline") {
      $locusId = shift @words;
      $sysName = shift @words;
      $desc = join(" ", @words);
    } elsif ($assembly->{gdb} eq "NCBI") {
      $locusId = shift @words;
      $desc = join(" ", @words);
      $desc =~ s/^MULTISPECIES: *//;
      # remove trailing organism descriptor
      $desc =~ s/ *\[[^\]]+\]$//;
      if (exists $assembly->{prot}{$locusId}) {
        # use locus tag and description from the feature table instead
        my $g = $assembly->{prot}{$locusId};
        $sysName = $g->{locus_tag};
        $desc = $g->{name} if $g->{name} ne "";
      }
    } elsif ($assembly->{gdb} eq "IMG") {
      $locusId = shift @words;
      $sysName = shift @words;
      $desc = join(" ", @words);
      # remove trailing organism descriptor
      $desc =~ s/ *\[[^\]]+\]$//;
    } elsif ($assembly->{gdb} eq "UniProt") {
      my @ids = split /[|]/, shift @words;
      shift @ids;
      ($locusId, $sysName) = @ids if @ids == 2;
      $desc = join(" ", @words);
      if ($desc =~ m/^(.*) OS=(.*)/) {
        $desc = $1;
        my $rest = $2;
        my $gn = $1 if $rest =~ m/ GN=(\S+) /;
        $desc .= " ($gn)" if $gn;
      }
    } else {
      die "Unknown genome database $assembly->{gdb}";
    }
    print $fhOut ">$orgId:$locusId $sysName $desc\n$seq\n";
  }
  close($fhIn) || die "Error reading $faaFile";
  return $orgId;
}

# Given a fetched assembly, such as from FetchAssembly::CacheAssembly,
# and a filehandle, write out the predicted proteins as a fasta file with
# header lines >orgId:locusId sysName desc
# Returns the orgId
sub WriteSixFrameAsOrgProteins {
  my ($assembly, $fhOut, $usearch) = @_;
  die "Invalid assembly" unless $assembly->{gdb} && $assembly->{gid};
  die ": or whitespace not allowed in assembly identifier $assembly->{gid}"
    if $assembly->{gid} =~ m/\s/ || $assembly->{gid} =~ m/:/;
  my $orgId = join("__", $assembly->{gdb}, $assembly->{gid});
  die "Not an executable: $usearch\n" unless -x $usearch;
  die "No such file: $assembly->{fnafile}\n" unless -e $assembly->{fnafile};
  my $xfile = "/tmp/$$.steps.aa6";
  # orfstyle 7 means allow ORF at beginning of sequence,
  # or immediately after a stop codon (no proper start codon),
  # and allow orfs to end at the end of a scaffold
  system($usearch, "-fastx_findorfs", $assembly->{fnafile}, "-aaout", $xfile,
         "-orfstyle", 7, "-mincodons", 30) == 0
           || die "usearch failed: $!";
  open(my $fhIn, "<", $xfile) || die "Cannot read $xfile from usearch\n";
  my $state = {};
  while (my ($header, $seq) = ReadFastaEntry($fhIn, $state)) {
    # The extra space is for the empty description
    my @parts = split / /, $header;
    my $scaffoldId = $parts[0];
    die "No scaffold in header $header from usearch -fastx_findorfs for $orgId\n"
      unless defined $scaffoldId && $scaffoldId ne "";
    @parts = split /[|]/, $header;
    my $partLast = $parts[-1];
    die "No position information in header $header from usearch -fastx_findorfs for $orgId\n"
      unless $partLast;
    die "Space in position information $partLast form usearch -fastx_findorgs for $orgId\n"
      if $partLast =~ m/ /;
    my $locusId = join("|", $scaffoldId, $partLast);
    print $fhOut ">$orgId:$locusId $locusId \n$seq\n";
  }
  unlink($xfile);
  return $orgId;
}

# Read a protein from a fasta file produced by WriteAssemblyAsOrgProteins
# Returns a hash of orgId, locusId, sysName, desc, aaseq
# If orgId is in the correct form, it also includes gdb and gid
# Set state to {} before the first invocation (as for ReadFastaEntry)
sub ReadOrgProtein($$) {
  my ($fh, $state) = @_;
  my ($header, $aaseq) = ReadFastaEntry($fh, $state);
  return undef unless defined $header;
  my %out = ( aaseq => $aaseq );

  # keep blank items at end in case sysName and/or desc is missing
  my @words = split / /, $header, -1;
  die "Not enough words in fasta header $header\n"
    unless @words >= 3; # locusSpec sysName desc
  my $locusSpec = shift @words;
  my ($orgId, $locusId) = ParseOrgLocus($locusSpec);
  $out{orgId} = $orgId;
  $out{locusId} = $locusId;
  $out{sysName} = shift @words;
  $out{desc} = join(" ", @words);
  my @pieces = split /__/, $out{orgId};
  $out{qw{gdb gid}} = @pieces if @pieces == 2;
  return \%out;
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

sub ReadOrgTable($) {
  my ($file) = @_;
  my @required_fields = qw{orgId genomeName};
  return ReadTable($file, \@required_fields);
}

# filename and list of required fields => list of hashes, each with field->value
# The list can be a single name or a reference to a list
# (It used to be an actual list; not sure how that went wrong with perl prototypes?)
sub ReadTable($*) {
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

# from orgId:locusId to orgId and locusId
# Allows ":" within locusId
sub ParseOrgLocus($) {
  my ($locusSpec) = @_;
  my @pieces = split /:/, $locusSpec;
  die "Not enough pieces in locus specifier $locusSpec"
    unless @pieces >= 2;
  my $orgId  = shift @pieces;
  my $locusId = join(":", @pieces);
  die "Invalid locus specifier $locusSpec" unless $orgId ne "" && $locusId ne "";
  return ($orgId, $locusId);
}

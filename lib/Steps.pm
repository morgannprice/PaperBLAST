package Steps;
require Exporter;
use strict;
our (@ISA,@EXPORT);
@ISA = qw(Exporter);
@EXPORT = qw(ReadSteps FetchUniProtSequence AssemblyToProt
             WriteAssemblyAsOrgProteins WriteSixFrameAsOrgProteins
             ReadOrgTable ReadOrgProtein ParseOrgLocus
             ReadReqs CheckReqs);
use LWP::Simple qw{get};
use LWP::UserAgent;

sub ReadSteps($); # stepfile => hash of steps, rules, ruleOrder
sub ReadSteps2($$); # stepfile, doimportflag => same result as ReadSteps

# The returned object includes:
# steps, a hash of step name to a hash that includes
#	i, name, desc, comment, and search, and optionally import
#	search is a list of queryType/value pairs (documented in bin/gapquery.pl)
#	use i to sort them by the same order as in the stepfile
#	import is set if this step was imported
# rules, a hash of rule name to a list of lists, each element being a step name or another rule name
#	The requirement is met if any of the sublists are met (OR at the top level)
#	The sublist is met if all of its components are met (AND at the 2nd level)
# ruleOrder, a list of rule names in order
# ruleComments, a hash of ruleId to comments.
# topComment, the top-level comment at the beginning of the steps file.
sub ReadSteps($) {
  my ($stepsFile) = @_;
  return ReadSteps2($stepsFile, 1); # 1 means allow imports
}

sub ReadSteps2($$) {
  # The purpose of the doimport flag is to prevent cycles A imports from B imports from A
  # When called with doimport=0, ignores import statements, and allows rules to be incomplete
  my ($stepsFile, $doimport) = @_;
  my %imports = (); # filepathname to step object
  my @stepTypes = qw{term hmm EC uniprot ignore ignore_hmm ignore_other curated};
  my %stepTypes = map { $_ => 1 } @stepTypes;
  my $steps = {}; # step name to hash of name, desc, and search, which is a list of pairs
  my $rules = {}; # rule name to list of lists, each element is a step name or another rule name
  my @ruleOrder = (); # list of rules in input order
  my %ruleComments = ();
  my $topComment = undef;
  open(my $fhSteps, "<", $stepsFile) || die "Cannot read $stepsFile";
  my $nSteps = 0;
  my @commentLines = (); # preceding comment lines. Emptied on lines that are not comments.
  while (my $line = <$fhSteps>) {
    $line =~ s/[\r\n]+$//;
    $line =~ s/\s+$//;
    if ($line =~ m/^#\s*(.*)$/) {
      push @commentLines, $1;
    } elsif ($line eq "") {
      $topComment = join(" ", @commentLines) if !defined $topComment;
      @commentLines = ();
    } elsif ($line =~ m/^(\S+):\s+(\S.*)$/) {
      my ($rulename, $pieces) = ($1, $2);
      my $newComment = join(" ", @commentLines);
      if (exists $ruleComments{$rulename}) {
        $ruleComments{$rulename} .= " " . $newComment;
      } else {
        $ruleComments{$rulename} = $newComment;
      }
      die "Cannot use step name $rulename as rule name in $stepsFile\n"
        if exists $steps->{$rulename};
      my @pieces = split /\s+/, $pieces;
      foreach my $piece (@pieces) {
        # In !doimport mode, some pieces could be missing
        die "Rule $rulename includes $piece which has not been defined yet in $stepsFile\n"
          unless exists $rules->{$piece} || exists $steps->{$piece} || ! $doimport;
        die "Rule $rulename cannot include itself in $stepsFile\n"
          if $piece eq $rulename;
      }
      push @ruleOrder, $rulename unless exists $rules->{$rulename};
      push @{ $rules->{$rulename} }, \@pieces;
    } elsif ($line =~ m/^import\s/) {
      my $toget = $line;
      $toget =~ s/^import\s+//;
      $toget =~ s/#.*//;
      $toget =~ s/\s+$//;
      my @pieces = split /:/, $toget;
      die "Invalid import line in $stepsFile\n$line\n(should be of the form\nimport stepsfile:step_or_rules\n"
        unless @pieces == 2;
      my ($importFileName, $importNames) = @pieces;
      die "Invalid import file $importFileName in $stepsFile\n"
        unless $importFileName =~ m/^[a-zA-Z0-9._-]+$/;
      my $importFilePathName = $importFileName;
      if ($stepsFile =~ m!/!) {
        $importFilePathName = $stepsFile;
        $importFilePathName =~ s!/[^/]+$!!;
        $importFilePathName = $importFilePathName . "/" . $importFileName;
      }
      die "Unknown import file $importFilePathName\nfrom line $line\nin steps file $stepsFile\n"
        unless -e $importFilePathName;
      my @importNames = split / /, $importNames;
      if ($doimport) {
        $imports{$importFilePathName} = ReadSteps2($importFilePathName, 0)
          if !exists $imports{$importFilePathName};
        my $stepsImport = $imports{$importFilePathName};

        # Add all the dependencies to the import list
        my %importNames = map { $_ => 1 } @importNames;
        my $nMaxRounds = 1000;
        my $nRound;
        for ($nRound = 0; $nRound < $nMaxRounds; $nRound++) {
          my $nImportOld = scalar(@importNames);
          foreach my $importName (@importNames) {
            if (exists $stepsImport->{rules}{$importName}) {
              foreach my $ruleInstance (@{ $stepsImport->{rules}{$importName}; }) {
                foreach my $rulePart (@$ruleInstance) {
                  if (!exists $importNames{$rulePart}) {
                    push @importNames, $rulePart;
                    $importNames{$rulePart} = 1;
                  }
                }
              }
            }
          }
          last if scalar(@importNames) == $nImportOld;
        }
        die "Too many dependencies" unless $nRound < $nMaxRounds;

        foreach my $importName (@importNames) {
          die "Steps file $stepsFile imports $importName from $importFileName, but $importName already exists\n"
            if exists $steps->{$importName} || exists $rules->{$importName};
          if (exists $stepsImport->{steps}{$importName}) {
            # Load the step.
            my $stepObj = $stepsImport->{steps}{$importName};
            $stepObj->{imported} = 1;
            $stepObj->{i} = $nSteps++; # renumber
            $steps->{$importName} = $stepObj;
          } elsif (exists $stepsImport->{rules}{$importName}) {
            # Load the rule, and, all dependencies
            my $ruleObj = $stepsImport->{steps}{$importName};
            $rules->{$importName} = $ruleObj;
            push @ruleOrder, $importName;
            $ruleComments{$importName} = $stepsImport->{ruleComments}{$importName}
              if exists $stepsImport->{ruleComments}{$importName};
          } else {
            die "Steps file $stepsFile imports unknown $importName from $importFilePathName\n";
          }
        }
      } else { # do not actually import
        foreach my $importName (@importNames) {
          $steps->{$importName} = {}; # empty object
        }
      }
    } elsif ($line =~ m/\t/) {
      my @F = split /\t/, $line;
      die "Invalid rule line: $line\n" unless @F >= 3;
      my $stepname = shift @F;
      die "Duplicate definition of step $stepname in $stepsFile\n"
        if exists $steps->{$stepname};
      die "Spaces or other whitespace within step name $stepname in $stepsFile are not allowed\n"
        if $stepname =~ m/\s/;
      my $desc = shift @F;
      my @search = ();
      foreach my $search (@F) {
        $search =~ m/^([a-zA-Z_]+):(.*)$/
          || die "Do not recognize search specifier $search for step $stepname in $stepsFile\n";
        my ($type, $value) = ($1,$2);
        die "Do not recognize the type in search specifier $search for step $stepname in $stepsFile\n"
          unless exists $stepTypes{$type};
        push @search, [$type, $value];
      }
      die "Invalid step has no search items: $line\n"
        unless @search > 0;
      $steps->{$stepname} = { 'name' => $stepname, 'desc' => $desc, 'search' => \@search,
                              'comment' => join(" ", @commentLines),
                              'i' => $nSteps++ };
    } else {
      die "Do not recognize line $line in $stepsFile as either rule or step\n";
    }
    @commentLines = () unless $line =~ m/^#/;
  }
  close($fhSteps) || die "Error reading $stepsFile";

  # Check for cyclic dependencies within the rules. Do a
  # a depth-first traversal starting at the rule all and
  # verify that it terminates.
  # (An older version of this check depended on the rule order instead.)
  die "No rule named all" unless exists $rules->{all};
  my @work = ("all");
  my $nMaxRounds = 10000;
  my $nRound = 0;
  while (@work > 0) {
    $nRound++;
    die "Depth-first traversal did not terminate -- cyclic depednencies!"
      if $nRound > $nMaxRounds;
    my $rule = shift @work;
    die unless exists $rules->{$rule};
    foreach my $instance (@{ $rules->{$rule} }) {
      foreach my $piece (@$instance) {
        next if exists $steps->{$piece};
        if (!exists $rules->{$piece}) {
          die "Unknown piece $piece in rule $rule" if $doimport;
        } else { # known rule
          push @work, $piece;
        }
      }
    }
  }

  return { 'steps' => $steps, 'rules' => $rules, 'ruleOrder' => \@ruleOrder,
           'topComment' => $topComment, 'ruleComments' => \%ruleComments };
}

# Returns sequence and description (just joining the DE lines)
sub FetchUniProtSequence($) {
  my ($id) = @_;
  # fetch the full entry because this allows either ids or accessions (i.e., Q72CU7_DESVH or Q72CU7) to be used
  my $URL = "https://www.uniprot.org/uniprot/${id}.txt";

  my $content;
  for (my $i = 0; $i < 3; $i++) {
    $content = get($URL);
    if (defined $content) {
      last;
    } else {
      sleep(1);
    }
  }
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
    } elsif ($assembly->{gdb} eq "local") {
      $locusId = shift @words;
      $desc = join(" ", @words);
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

# Given a requirements file (requires.tsv) and a hash of pathwayId => step object,
# returns a reference to a list of requirements, each of the form:
#   pathway, rule (defaults to "all"),
#   requiredPath, requiredRule (defaults to "all") or requiredStep,
#   not, comment, and reqSpec
sub ReadReqs($$) {
  my ($requiresFile, $pathways) = @_;
  my @header = qw{rule requires comment};
  my @requires = ReadTable($requiresFile, \@header);

  my @out = ();
  foreach my $req (@requires) {
    my $ruleSpec = $req->{rule}; # pathway or pathway:rule
    my ($pathwayId, $rule) = split /:/, $ruleSpec;
    $rule = "all" if !defined $rule;
    die "Unknown pathway $pathwayId in requirements $requiresFile\n"
      unless exists $pathways->{$pathwayId};
    die "Unknown rule $rule for pathway $pathwayId in requirements $requiresFile\n"
      if !exists $pathways->{$pathwayId}{rules}{$rule};
    my $reqSpec = $req->{requires};
    my $not = $reqSpec =~ m/^[!]/ ? 1 : 0;
    my $reqSpec2 = $reqSpec;
    $reqSpec2 =~ s/^[!]//;
    my ($reqPath, $reqPart) = split /:/, $reqSpec2;
    my $out = { 'pathway' => $pathwayId, 'rule' => $rule,
                'reqSpec' => $reqSpec, 'requiredPath' => $reqPath,
                'not' => $not,
                'comment' => $req->{comment}
              };
    $out->{comment} =~ s/^# *//;
    $reqPart = "all" if !defined $reqPart;
    if (exists $pathways->{$reqPath}{rules}{$reqPart}) {
      $out->{requiredRule} = $reqPart;
    } elsif (exists $pathways->{$reqPath}{steps}{$reqPart}) {
      $out->{requiredStep} = $reqPart;
      die "Not is not supported for required steps -- see $pathwayId $rule and $reqSpec\n"
        if $not;
    } else {
      die "Unknown required part $reqPart for pathway $reqPath in requirements $requiresFile\n"
    }
    push @out, $out;
  }
  return \@out;
}

# Given the list of sum.rules and sum.steps rows for an organism,
# and the requirements from the Requirement table, returns a list of warnings
# that are relevant.
sub CheckReqs($$$) {
  my ($sumRules, $sumSteps, $reqs) = @_;
  die unless $sumRules && $sumSteps && defined $reqs;

  my %ruleScores = (); # pathway => rule => score row
  foreach my $row (@$sumRules) {
    die "Duplicate score for rule $row->{pathway} $row->{rule}\n"
      if exists $ruleScores{ $row->{pathway} }{ $row->{rule} };
    $ruleScores{ $row->{pathway} }{ $row->{rule} } = $row;
  }

  my %stepScores = (); # pathway => step => score row
  foreach my $row (@$sumSteps) {
    die "Duplicate score for step $row->{pathway} $row->{step}\n"
      if exists $stepScores{ $row->{pathway} }{ $row->{step} };
    $stepScores{ $row->{pathway} }{ $row->{step} } = $row;
  }

  # List all the rules that are on best paths, working down from all
  my %onBestPath = (); # pathway => rule or step => 1
  while (my ($pathway, $ruleHash) = each %ruleScores) {
    my @work = ("all");
    while (@work > 0) {
      my $rule = shift @work;
      next if exists $onBestPath{$pathway}{$rule};
      $onBestPath{$pathway}{$rule} = 1;
      if (exists $ruleHash->{$rule}) {
        push @work, split / /, $ruleHash->{$rule}{path};
      }
    }
  }
  my @warn = ();
  foreach my $req (@$reqs) {
    my $pathwayId = $req->{pathwayId} || die;
    my $ruleId = $req->{ruleId} || die;
    # Ignore warnings that are not about the best path
    next unless exists $onBestPath{$pathwayId}{$ruleId};

    # Do not warn if already low confidence
    my $ruleS = $ruleScores{$pathwayId}{$ruleId} || die;
    next if $ruleS->{nLo} > 0;
    my $reqPathwayId = $req->{requiredPathwayId};
    if ($req->{requiredStepId} ne "") {
      die "Step requirements with isNot are not supported"
        if $req->{isNot};
      # Do not warn if the required step is high confidence
      my $stepS = $stepScores{$reqPathwayId}{ $req->{requiredStepId} } || die;
      next if $stepS->{score} eq "2";
    } elsif ($req->{requiredRuleId} ne "") {
      my $reqRuleId = $req->{requiredRuleId};
      my $reqS = $ruleScores{$reqPathwayId}{$reqRuleId} || die;
      if ($req->{isNot}) {
        # ! means warn if requiredRule has no low-confidence steps and is on the best path for $reqPath
        next unless $reqS->{nLo} == 0 && exists $onBestPath{$reqPathwayId}{$reqRuleId};
      } else {
        # Warn if requiredRule has any low- or medium-confidence steps
        next if $reqS->{nLo} == 0 && $reqS->{nMed} == 0;
      }
    } else {
      die "No requiredStep or requiredRule";
    }
    push @warn, $req;
  }
  return \@warn;
}

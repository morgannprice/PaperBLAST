#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw{$Bin};
use File::Which;
use IO::Handle;                 # for autoflush

my $usage = <<END
buildLitDb.pl -dir dir -snippets snippetsFile -prefix prefix1 ... prefixN
        [ -rif generif_tab.rif ]
	[ -curated curated_files ]
        [ -blastdir blast ]

Each prefix should be the output of parseEuropePMCHits.

The snippets file should be the output of buildSnippets.pl or combineSnippets.pl
	[the snippets.access file must exist as well]

The rif file should be from generifTables.pl (generif_tab should also
be one of the prefixes)

The (optional) curated arguments must be in curated_parsed format,
which is tab-delimited with fields dbId, proteinId, secondary_id,
short_name, description, organism, sequence, comment, and pubmed ids
(comma separated).

The -blast argument specifies the path to the blastall and formatdb executables
(or else there must be a blast/ directory of the working directory or the executable directory).

Creates dir/litsearch.db, dir/litsearch.faa, and formats the BLAST database
END
  ;

# sqlite3 expects CVS format, not exactly tab delimited format
# So, need to replace any " with "" and surround the field with quotes.
sub csv_quote($);

{
  my ($dir, $snippetsFile, $rifFile);
  my $blastdir;
  my @prefixes = ();
  my @curated = ();
  die $usage
    unless GetOptions('dir=s' => \$dir,
                      'snippets=s' => \$snippetsFile,
                      'rif=s' => \$rifFile,
                      # the {1,} syntax means 1 or more values
                      'curated=s{1,}' => \@curated,
                      'prefix=s{1,}' => \@prefixes)
      && @ARGV == 0;
  die "Must specify -dir\n" unless defined $dir;
  die "Not a directory: $dir\n" unless -d $dir;
  die "Must specify -prefix\n" unless @prefixes > 0;
  die "No such file: $snippetsFile\n" if ! -e $snippetsFile;
  die "No such file: $snippetsFile.access\n" if ! -e "$snippetsFile.access";
  die "No such file: $rifFile" if defined $rifFile && ! -e $rifFile;

  if (!defined $blastdir) {
    if (-d "blast") {
      $blastdir = "blast";
    } elsif (-d "$Bin/blast") {
      $blastdir = "$Bin/blast";
    } else {
      die "No -blastdir argument and not in . or $Bin\n";
    }
  }
  die "No such executable: $blastdir/blastall\n" unless -x "$blastdir/blastall";
  die "No such executable: $blastdir/formatdb\n" unless -x "$blastdir/formatdb";

  my $schema = "$Bin/litsearch.sql";
  die "Cannot find schema: not in $schema" unless -e $schema;

  # Check that all prefix files exist
  foreach my $prefix (@prefixes) {
    die "No such file: $prefix.papers\n" unless -e "$prefix.papers";
    die "No such file: $prefix.queries\n" unless -e "$prefix.queries";
    die "No such file: $prefix.faa\n" unless -e "$prefix.faa";
  }

  # Check that all curated files exist
  foreach my $curated (@curated) {
    die "No such file: $curated\n" unless -e $curated;
  }

  my $sqldb = "$dir/litsearch.db";
  unlink($sqldb);

  my $faafile = "$dir/litsearch.faa";

  system("sqlite3 $sqldb < $schema");

  # Read in the genes with text-mined/GeneRIF links and write just the Gene table for now.
  # These are to prevent duplicates
  my %locusLen = (); # locusId => protein length
  my %aaseq = (); # geneId => sequence if already saw a sequence
  open(GENE, ">", "$dir/Gene") || die "Cannot write to $dir/Gene";
  foreach my $prefix (@prefixes) {
    open (QUERYIN, "<", "$prefix.queries") || die "Cannot read $prefix.queries";
    while (my $line = <QUERYIN>) {
      chomp $line;
      my ($geneId, $organism, $protlen, $desc) = split /\t/, $line;
      die "Invalid input" unless defined $protlen && $protlen =~ m/^\d+$/;
      die "Mismatching length for $geneId in $prefix"
        if exists $locusLen{$geneId} && $locusLen{$geneId} ne $protlen;
      next if exists $locusLen{$geneId};
      $locusLen{$geneId} = $protlen;
      $desc = "" if !defined $desc;
      print GENE join("\t", $geneId, $organism, $protlen, $desc) . "\n";
    }
    close(QUERYIN) || die "Error reading $prefix.queries";

    open(FAAIN, "<", "$prefix.faa");
    while (my $header = <FAAIN>) {
      chomp $header;
      die "Invalid faa file $prefix.faa" unless $header =~ m/^>(.*)$/;
      my $geneId = $1;
      die "Unexpected gene $geneId in $prefix.faa" unless exists $locusLen{$geneId};
      my $seq = <FAAIN>;
      die "Invalid faa file $prefix.faa" unless defined $seq;
      chomp $seq;
      die "Incorrect length for $geneId in $prefix.faa" unless length($seq) == $locusLen{$geneId};
      if (exists $aaseq{$geneId}) {
        print STDERR "Warning: non-matching sequences for $geneId in $prefix.faa\n"
          if $aaseq{$geneId} ne $seq;
      } else {
        $aaseq{$geneId} = $seq;
      }
    }
    close(FAAIN) || die "Error reading $prefix.faa";
  } # end loop over prefixes
  close(GENE) || die "Error writing to $dir/Gene";

  # %linkSupport records which gene-paper links are supported by snippet or by GeneRIF
  # and should be kept
  my %linkSupport = ();          # pmcId::pmId => queryId => 1

  # Read in the snippets and write the Snippet table
  open(OUT, ">", "$dir/Snippet") || die "Cannot write to $dir/Snippet";
  my $nSuppressedByCase = 0;
  open(IN, "<", $snippetsFile) || die "Cannot read $snippetsFile";
  while (my $line = <IN>) {
    chomp $line;
    my ($pmcId, $pmId, $queryTerm, $queryId, $snippet) = split /\t/, $line;
    die "Not enough fields in $snippetsFile\n$line" unless defined $snippet;
    # if query term is like c2233 then consider it risky and require case to match
    if ($queryTerm =~ m/^[a-zA-Z]\d+$/ && $snippet !~ m/$queryTerm/) {
      $nSuppressedByCase++;
    } else {
      print OUT join("\t", $queryId, $queryTerm, $pmcId, $pmId, &csv_quote($snippet))."\n";
      $linkSupport{ join("::", $pmcId, $pmId) }{ $queryId } = 1;
    }
  }
  print STDERR "Wrote $dir/Snippet, suppressed $nSuppressedByCase snippets with risky locus tags and the wrong case\n";
  close(OUT) || die "Error writing to $dir/Snippet";

  # Read in the GeneRIFs and write the GeneRIF table
  # and set up RIFs
  open(OUT, ">", "$dir/GeneRIF") || die "Cannot write to $dir/GeneRIF";
  if (defined $rifFile) {
    open(IN, "<", $rifFile) || die "Cannot read $rifFile";
    while(my $line = <IN>) {
      chomp $line;
      my ($pmcId, $pmId, $queryIdShort, $queryId, $rif) = split /\t/, $line;
      die "Invalid inputin $rifFile" unless defined $rif && $rif;
      $linkSupport{ join("::", $pmcId, $pmId) }{ $queryId } = 1;
      print OUT join("\t", $queryId, $pmcId, $pmId, &csv_quote($rif))."\n";
    }
    close(IN) || die "Error reading $rifFile";
  }
  close(OUT) || die "Error writing to $dir/GeneRIF";

  # Write the PaperAccess table
  my %full = ();                # paperId (pmcId::pmId) => full text
  open(OUT, ">", "$dir/PaperAccess") || die "Cannot write to $dir/PaperAccess";
  open(IN, "<", "$snippetsFile.access") || die "Cannot read $snippetsFile.access";
  my %seen = ();
  while (my $line = <IN>) {
    chomp $line;
    my ($pmcId, $pmId, $access) = split /\t/, $line;
    die "Not enough fields in $snippetsFile.access\n$line" unless defined $access;
    my $key = $pmcId."::".$pmId;
    die "Duplicate row in $snippetsFile.access for $pmcId $pmId"
      if exists $seen{$key};
    $seen{$key} = 1;
    print OUT join("\t", $pmcId, $pmId, $access)."\n";
    $full{$key} = 1 if $access eq "full";
  }
  close(IN) || die "Error reading $snippetsFile.access";
  close(OUT) || die "Error writing to $dir/PaperAccess";

  # Write the GenePaper table with the gene/paper combinations that were kept
  my %geneHasPaper = ();
  my $nIgnoreNoSnippet = 0;
  open(GENEPAPER, ">", "$dir/GenePaper") || die "Cannot write to $dir/GenePaper";
  foreach my $prefix (@prefixes) {
    open(PAPERIN, "<", "$prefix.papers") || die "Cannot read $prefix.papers";
    while (my $line = <PAPERIN>) {
      chomp $line;
      my ($queryId, $queryTerm, $pmcId, $pmId, $doi, $title, $authors, $journal, $year, $isOpen) = split /\t/, $line;
      die "Invalid papers line in $prefix.papers:\n$line\n" unless defined $isOpen;
      die "Unexpected gene $queryId in $prefix.papers" unless exists $aaseq{$queryId};
      # Do not bother to remove duplicates
      my $key = join("::", $pmcId, $pmId);
      if (exists $full{$key} && !exists $linkSupport{$key}{$queryId}) {
        $nIgnoreNoSnippet++;
      } else {
        $geneHasPaper{$queryId} = 1;
        print GENEPAPER join("\t", $queryId, $queryTerm, $pmcId, $pmId, $doi, &csv_quote($title), $authors, $journal, $year, $isOpen)."\n";
      }
    }
    close(PAPERIN) || die "Error reading $prefix.papers";
  }
  close(GENEPAPER) || die "Error writing to $dir/GenePaper";
  print STDERR "Read paper links: ignored $nIgnoreNoSnippet cases with full text and no snippet. " . scalar(keys %geneHasPaper) . " cases left\n";

  # Save the sequences for the genes that still have a link
  my $nSkipGene = 0;
  open(FAA, ">", $faafile) || die "Cannot write to $faafile";
  while (my ($geneId, $seq) = each %aaseq) {
    if (exists $geneHasPaper{$geneId}) {
      print FAA ">$geneId\n$seq\n";
    } else {
      $nSkipGene++;
    }
  }

  # Build the curated tables, and also write out those proteins' sequences
  open(CGENE, ">", "$dir/CuratedGene") || die "Cannot write to $dir/CuratedGene";
  open(CPAPER, ">", "$dir/CuratedPaper") || die "Cannot write to $dir/CuratedPaper";
  my %curatedId = (); # db => id => 1 -- there should not be any duplicates
  foreach my $cfile (@curated) {
    open(IN, "<", $cfile) || die "Cannot read $cfile";
    while(my $line = <IN>) {
      chomp $line;
      my ($db, $protId, $id2, $name, $desc, $org, $seq, $comment, $pmIds) = split /\t/, $line, -1;
      die "Invalid line\n$line\nin $cfile" unless defined $pmIds;
      die "No database in $line" if $db eq "";
      die "No protein id in $line" if $protId eq "";
      die "No protein sequence in $line" if $seq eq "";
      $seq = uc($seq); # allow lower case
      $seq =~ m/^[A-Z*]+$/ || die "Invalid sequence $seq in $cfile";
      die "Duplicate curated entry $db $protId in $cfile"
        if exists $curatedId{$db}{$protId};
      die "Invalid pubmed ids $pmIds in $cfile" unless $pmIds =~ m/^[0-9,]*$/;
      $curatedId{$db}{$protId} = 1;
      my $combid = join("::", $db, $protId);
      if ($desc eq "") {
        $desc = $name || $id2;
        die "No description, name, or secondary id in $line" if $desc eq "";
        print STDERR "Warning: no description for curated protein $protId from $db\n";
      }
      print FAA ">$combid\n$seq\n";
      print CGENE join("\t", $db, $protId, $id2, $name, &csv_quote($desc),
                       $org, length($seq), &csv_quote($comment))."\n";
      my @pmIds = split /,/, $pmIds;
      my %pmIds = (); # check for duplicates
      foreach my $pmId (@pmIds) {
        die "Invalid pubmed id $pmId for $db $protId in $cfile" unless $pmId =~ m/^\d+$/;
        die "Duplicate pubmed id for $db $protId in $cfile" if exists $pmIds{$pmId};
        $pmIds{$pmId} = 1;
        print CPAPER join("\t", $db, $protId, $pmId)."\n";
      }
    }
    close(IN) || die "Error reading $cfile";
  }
  close(CGENE) || die "Error writing $dir/CuratedGene";
  close(CPAPER) || die "Error writing $dir/CuratedPaper";

  close(FAA) || die "Error writing to $faafile";
  print STDERR "Wrote $faafile, with $nSkipGene genes skipped because no remaining links to papers\n";

  open(SQLITE, "|-", "sqlite3", "$sqldb") || die "Cannot run sqlite3 on $sqldb";
  autoflush SQLITE 1;
  print SQLITE ".mode tabs\n";
  my @tables = qw{GenePaper Gene Snippet PaperAccess GeneRIF CuratedGene CuratedPaper};
  foreach my $table (@tables) {
    print STDERR "Loading table $table\n";
    print SQLITE ".import $dir/$table $table\n";
  }

  my $reportCmds = <<END
SELECT 'nGenesWithPubs', COUNT(DISTINCT geneId) FROM Gene;
SELECT 'nPubMedIds', COUNT(DISTINCT pmId) FROM GenePaper;
SELECT 'nPaperLinks', COUNT(*) FROM GenePaper;
SELECT 'nSnippets', COUNT(*) FROM Snippet;
SELECT 'nGeneRIF', COUNT(*) FROM GeneRIF;
SELECT 'nFullPapers', COUNT(*) FROM PaperAccess WHERE access="full";
SELECT 'nCurated', COUNT(*) FROM CuratedGene;
SELECT 'nCuratedPapers', COUNT(DISTINCT pmId) FROM CuratedPaper;
END
    ;
  print SQLITE $reportCmds;
  close(SQLITE) || die "Error running sqlite3 commands\n";

  system("$Bin/derepSequences.pl", "-dir", $dir) == 0
    || die "derepSequences.pl failed: $!";

  system("$blastdir/formatdb", "-p", "T", "-i", "$dir/uniq.faa", "-o", "T") == 0
    || die "formatdb failed";

  print STDERR "Success\n";
}

sub csv_quote($) {
  my ($in) = @_;
  return $in unless $in =~ m/"/;
  $in =~ s/"/""/g;
  return '"' . $in . '"';
}

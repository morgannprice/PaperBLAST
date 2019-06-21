#!/usr/bin/perl -w
#######################################################
## genomeSearch.cgi -- proteins in a genome that are
##   similar to a curated protein whose description
##   matches a text query
##
## Copyright (c) 2018 University of California
##
## Authors: Morgan Price
#######################################################
#
# Optional CGI parameters:
# query -- what term to search for
# word -- if non-empty, report whole word matches only
# gdb -- one of the genome sources supported by FetchAssembly:
#	"NCBI", "IMG", "UniProt", "MicrobesOnline", "FitnessBrowser"
# gid -- which organism or assembly in that database
# gquery -- terms for searching for a genome
# file -- an uploaded file with protein sequences in FASTA format
#
# doupload -- set to show the upload page
# Search -- set if the search button was pressed
#	[used to handle failed uploads or if no query was entered]

use strict;
use CGI qw(:standard Vars start_ul end_ul);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use Time::HiRes qw{gettimeofday};
use DBI;
use URI::Escape; # for uri_escape()
use HTML::Entities;
use IO::Handle; # for autoflush
use List::Util qw{min max};
use lib "../lib";
use pbutils; # for ReadFastaEntry(), Curated functions
use FetchAssembly; # for GetMatchingAssemblies(), CacheAssembly(), warning(), fail(), etc.
use pbweb;

sub query_fields_html;

# maximum size of posted data, in bytes
my $maxMB = 100;
$CGI::POST_MAX = $maxMB*1024*1024;
my $maxseqsK = 100;
my $maxseqs = $maxseqsK * 1000;
my $maxNtLen = 30 * 1000 * 1000;
my $maxseqsComma = "$maxseqsK,000";
my $maxEvalue = 0.01;
my $maxHitsEach = 3; # additional hits are hidden until the expander is clicked
my $nCollapseSet = 0; # keeping track of which expander goes with what

my $minCodons = 30; # for reading frames

# from a defline to a brief description and a link to the source
sub ProteinLink($);
sub SixFrameLink($$);

# Given a ublast file, a hash with the valid subject ids, and idToChit to describe the queries and their lengths,
# returns a reference to a list of hashes. Each row includes
# input (the subject), hit (the query), identity, coverage, score, irange and hrange (each as begin:end),
# and chits (the chits object for this query)
sub ParseUblast($$$);

# given a sequence id, sequence, hits, and whether or not this is the 6-frame translation, print the hits
sub PrintHits($$$$);

my $base = "../data";
my $sqldb = "$base/litsearch.db";
my $dbh = DBI->connect("dbi:SQLite:dbname=$sqldb","","",{ RaiseError => 1 }) || die $DBI::errstr;
my $usearch = "../bin/usearch";
die "No such executable: $usearch" unless -x $usearch;
my $blastdir = "../bin/blast";
die "No such directory: $blastdir" unless -d $blastdir;
my $blastdb = "$base/uniq.faa";
die "No such file: $blastdb" unless -e $blastdb;

my $cgi=CGI->new;

my $gdb = $cgi->param('gdb');
my $gid = $cgi->param('gid');
die "Invalid genome identifier $gid\n" if defined $gid && $gid !~ m/^[0-9A-Za-z_:.-]*$/;
my $gquery = $cgi->param('gquery');

my $upfile = $gid ? undef : $cgi->param('file');
my $query = $cgi->param('query');
my $word = $cgi->param('word');

my @gdbs = ("NCBI", "IMG", "UniProt", "MicrobesOnline", "FitnessBrowser", "local");
my %gdb_labels1 = ("NCBI" => "NCBI assemblies",
                   "UniProt" => "UniProt proteomes",
                   "IMG" => "JGI/IMG genomes", "FitnessBrowser" => "Fitness Browser genomes",
                   "local" => "Uploaded proteome"); # from other tools
my %gdb_labels = map { $_ => exists $gdb_labels1{$_} ? $gdb_labels1{$_} : "$_ genomes"} @gdbs;
die "Unknown genome database: $gdb\n"
  if $gdb && !exists $gdb_labels{$gdb};

&start_page('title' => "Curated BLAST",
           'banner' => "Curated BLAST for Genomes",
           'bannerURL' => "genomeSearch.cgi");
autoflush STDOUT 1; # show preliminary results

# A symbolic link to the Fitness Browser data directory is used (if it exists)
# to allow access to fitness browser genomes.
# That directory must include feba.db (sqlite3 database) and aaseqs (in fasta format)
SetFitnessBrowserPath("../fbrowse_data");

# Should include .JGI.info
SetPrivDir("../private");

my $tmpDir = "../tmp";
my $procId = $$;
my $timestamp = int (gettimeofday() * 1000);
my $basefile = $tmpDir . "/" . $procId . $timestamp;

my $assembly; # the fetched (and cached) assembly
my $genomeName;
my $fh_up;

if ($gdb && $gquery) {
  print p("Searching", $gdb_labels{$gdb}, "for", "'" . HTML::Entities::encode($gquery) . "'"), "\n";
  my @rows = GetMatchingAssemblies($gdb, $gquery);
  my $limit = GetMaxNAssemblies();
  if (@rows > 0) {
    my $desc = "Found " . scalar(@rows) . " assemblies";
    if (@rows == 1) {
      $desc = "Found 1 assembly";
    } elsif (@rows >= $limit) {
      $desc = "Found the first " . scalar(@rows) . " matching assemblies";
    }
    $desc .= ", please choose one" if @rows > 1;
    print p($desc . ":"), "\n";
    print start_form(-method => 'get', -action => 'genomeSearch.cgi'),
      hidden(-name => 'gdb', -value => $gdb, -override => 1);
    foreach my $row (@rows) {
      my $checked = @rows == 1 ? "CHECKED" : "";
      print qq{<INPUT type="radio" NAME="gid" VALUE="$row->{gid}" $checked />},
        a({ -href => $row->{URL} }, $row->{genomeName} ),
        " ", small("(" . $row->{gid} . ")"), br();
    }
    print query_fields_html(),
      p(submit(-name => 'Search', -value => 'Search in selected genome')),
      end_form;
  } else {
    print p("Sorry, no matching genomes were found.");
  }
  print p("Try", a({ -href => "genomeSearch.cgi?gdb=$gdb" }, "another genome"));
  finish_page();
} elsif ($gdb && $gid && $query) {
  $assembly = CacheAssembly($gdb, $gid, "../tmp")
    || fail("Cannot fetch assembly $gid from database $gdb");
  $genomeName = $assembly->{genomeName};
  my $link2 = $assembly->{gid};
  $link2 = a({ -href => $assembly->{URL2} }, $link2) if exists $assembly->{URL2};
  print p("Searching in", a({-href => $assembly->{URL} }, $genomeName),
         small("(" . $link2 . ")")),
         "\n";
  # finish searching down below
} elsif ($gdb && $gid) {
  # assembly chosen but no query was entered
  $assembly = CacheAssembly($gdb, $gid, "../tmp")
    || fail("Cannot fetch assembly $gid from database $gdb");
  warning("Please enter a search term")
    if $cgi->param('Search');
  print start_form(-method => 'get', -action => 'genomeSearch.cgi'),
    hidden(-name => 'gid', -value => $gid, -override => 1),
    hidden(-name => 'gdb', -value => $gdb, -override => 1),
    p("Search in",
      a({ -href => $assembly->{URL} }, $assembly->{genomeName}),
      small("(" . $assembly->{gid} . ")"),
      "from",  $gdb_labels{$gdb},
      "or try",
      a({ -href => "genomeSearch.cgi" }, "another genome")),
    query_fields_html(),
    p(submit(-name => 'Search', -value => 'Search in selected genome')),
    end_form;
  finish_page();
} elsif ($upfile && $query) {
  $genomeName = "uploaded file";
  $fh_up = $upfile->handle;
  die unless $fh_up;
} elsif ($cgi->param('doupload') || $cgi->param('uploading')) {
  warning("Please choose a FASTA file to analyze")
    if $cgi->param('uploading') && ! $upfile;
  warning("Please enter a query")
    if $cgi->param('uploading') && ! $query;
  print
    start_form( -autocomplete => 'on', -name => 'upload', -method => 'POST', -action => 'genomeSearch.cgi'),
      p("1. Upload amino acid or nucleotide sequences in FASTA format",
        br(),
        "(up to $maxMB MB and up to $maxseqsComma sequences)"),
      filefield(-name=>'file', -size=>50),
      query_fields_html(2),
      p(submit(-name => 'uploading', -value => 'Search')),
      end_form,
      p("Or", a({ -href => "genomeSearch.cgi"}, "search"), "for a genome");
  finish_page();
} else {
  warning("Please enter an organism name")
    if $cgi->param('findgenome');
  print
    GetMotd(),
    p("Given a genome and a query, find characterized proteins whose descriptions match the query,",
      "and then search the genome for homologs of those proteins",
      "(" . a({ -href => "genomeSearch.cgi?gdb=FitnessBrowser&gid=PS&query=perchlorate"},
              "example") . ")."),
    start_form( -autocomplete => 'on', -name => 'input', -method => 'GET', -action => 'genomeSearch.cgi'),
    p("Genome database to search:", 
      popup_menu(-name => 'gdb', -values => \@gdbs, -labels => \%gdb_labels, -default => $gdbs[0])),
    p(textfield(-name => 'gquery', -value => '', -size => 50, -maxlength => 200)),
    p(small("Example:", a({-href => "genomeSearch.cgi?gdb=NCBI&gquery=Azospira"}, "Azospira"))),
    p(submit(-name => "findgenome", -value => 'Find Genome')),
    end_form,
    p(br(),
      "Or",
      b(a({-href => "genomeSearch.cgi?doupload=1"}, "upload")),
      "a genome or proteome in fasta format."),
    p("Curated BLAST relies on curated descriptions of over 100,000 experimentally-characterized proteins from",
      a({-href => "http://uniprot.org",
         -title => "UniProtKB/Swiss-Prot (the manually annotated part of UniProt)"}, "Swiss-Prot").",",
      a({-href => "http://www.brenda-enzymes.org/index.php",
         -title => "The Comprehensive Enzyme Information System"}, "BRENDA").",",
      a({-href => "http://www.cazy.org/",
         -title => "the Carbohydrate-Active EnZymes database"}, "CAZy").",",
      a({-href => "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3245046/",
         -title => "Characterized Protein Database"}, "CharProtDB").",",
      a({-href => "http://metacyc.org/",
        -title => "MetaCyc Metabolic Pathway Database"}, "MetaCyc").",",
      a({-href => "http://ecocyc.org",
         -title => "EcoCyc: Encyclopedia of E. coli Genes and Metabolic Pathways"}, "EcoCyc").",",
      a({-href => "http://rebase.neb.com/rebase/rebase.html",
        -title => "The Restriction Enzyme Database"}, "REBASE").",",
      "and the",
      a({-href => "http://fit.genomics.lbl.gov/",
        -title => "Reannotations from genome-wide fitness data"}, "Fitness Browser")."."),
    h3("Related tools"),
    start_ul(),
    li(a({ -href => "curatedSearch.cgi" }, "Search for curated proteins by keyword")),
    li(a({ -href => "litSearch.cgi" }, "PaperBLAST: Find papers about a protein or its homologs")),
    li(a({ -href => 'gapView.cgi', -title => "Annotate the amino acid biosynthesis pathways in a genome" },
         "GapMind for amino acid biosynthesis")),
    end_ul(),
    h3("How it works"),
    p("Curated BLAST finds curated proteins whose descriptions match the query",
      "and uses", a({-href => "https://www.drive5.com/usearch/"}, "ublast"),
      "to compare these to the annotated proteins in the genome.",
      "Then, it uses ublast to compare the curated proteins to all possible proteins",
      "in the genome (the six-frame translation)."),
    p("For details see the",
      a({ -href => "https://msystems.asm.org/content/4/2/e00072-19.full",
          -title => "Curated BLAST for Genomes, mSystems 2019"}, "paper from 2019"),
      "or the",
      a({ -href => "https://github.com/morgannprice/PaperBLAST"}, "code").".",
      "Changes since the paper:",
      start_ul(),
      li(qq{Shows hits to unannotated reading frames if they hit regions of the curated
            protein that were not hit by any annotated protein
            (or if they score more highly).}),
      end_ul() );
  finish_page();
}

# Do the actual query -- first, validate the input
fail("No genome fetched") unless defined $genomeName;
die "No query\n" unless $query;

my $isNuc; # Is the main sequence set protein?

my $faaFile;
my $faaFileDesc;
my $fh;
if ($assembly) {
  $faaFile = $assembly->{faafile} || die "No faa file";
  $faaFileDesc = "protein fasta file for $genomeName";
  open($fh, "<", $faaFile) || die "Cannot read $faaFile\n";
} else {
  $fh = $fh_up;
  $faaFileDesc = "uploaded file";
}

my %seqs = (); # The main sequence set
my $state = {};
my $totlen = 0;
my $nNucChar = 0;
while (my ($header, $sequence) = ReadFastaEntry($fh,$state,"return_error")) {
  fail("Duplicate sequence for $header in $faaFileDesc") if exists $seqs{$header};
  fail(". found in sequence for $header in $faaFileDesc") if $sequence =~ m/[.]/;
  fail("- found in sequence for $header in $faaFileDesc") if $sequence =~ m/-/;
  $sequence =~ s/[*]//g;
  fail("Invalid/empty sequence for $header in $faaFileDesc") if $sequence eq "";
  $seqs{$header} = $sequence;
  $totlen += length($sequence);
  $nNucChar += ($sequence =~ tr/ACGTUN//);
}
close($fh) || fail("Error reading $faaFileDesc");
if (exists $state->{error}) {
  warning("Invalid fasta upload");
  warning($state->{error});
  print p("Please try", a({-href => "genomeSearch.cgi?doupload=1"}, "another upload"));
  finish_page();
}

fail("Too many sequences in $faaFileDesc. The limit is $maxseqs")
  if scalar(keys %seqs) > $maxseqs;
fail("No sequences in $faaFileDesc") if scalar(keys %seqs) == 0;
fail("Sorry, $faaFileDesc is too long. Please choose a different genome or proteome. The limit is 30 MB")
    if ($totlen >= $maxNtLen);
my $fracNuc = $nNucChar / $totlen;
$isNuc = $fracNuc >= 0.9;
if ($assembly && $isNuc) {
  warning("$faaFileDesc seems to contain nucleotide sequences instead!");
  $isNuc = 0;
}

print p("Found", scalar(keys %seqs),
        $isNuc ? "nucleotide" : "protein",
        "sequences in $upfile")."\n"
  if $upfile;

# Find relevant sequences in CuratedGene
my $URLnoq = $upfile ? "genomeSearch.cgi?doupload=1" : "genomeSearch.cgi?gdb=$gdb&gid=" . uri_escape($gid);
my $maxhits = 1000;
my $chits = CuratedMatch($dbh, $query, $maxhits+1);
my $quotedquery = HTML::Entities::encode($query);
if (@$chits > $maxhits) {
  print p(qq{Sorry, too many curated entries match the query '$quotedquery'. Please try},
          a({ -href => $URLnoq }, "another query").".");
  finish_page();
}
if (@$chits == 0) {
  print p(qq{None of the curated entries in PaperBLAST's database match '$quotedquery'. Please try},
          a({ -href => $URLnoq }, "another query") . ".");
  finish_page();
}
if ($word) {
  # filter for whole-word hits
  $chits = CuratedWordMatch($chits, $query);
  if (@$chits == 0) {
    print p(qq{None of the curated entries in PaperBLAST's database match '$quotedquery' as complete words. Please try},
            a({ -href => $URLnoq }, "another query") . ".");
    finish_page();
  }
}

my $wordstatement = $word ? " as complete word(s)" : "";
my $csURL = "curatedSearch.cgi?query=$quotedquery"
  . "&word=" . ($word ? 1 : 0);
print p("Found",
        a({ -href => $csURL }, scalar(@$chits), "curated entries"),
        qq{in PaperBLAST's database that match '$quotedquery'${wordstatement}.\n});

my $chitsfaaFile = "$basefile.chits.faa";
my $seqFile = "$basefile.seq";
my $ublastFile = "$basefile.ublast";

# Fetch the curated sequences that match
my %idToChit = (); # sequence identifier to curated gene hit(s)
foreach my $hit (@$chits) {
  my $seqid = $hit->{db} . "::" . $hit->{protId};
  my $uniqid = IdToUniqId($dbh, $seqid);
  push @{ $idToChit{$uniqid} }, $hit;
}
my @uniqids = sort keys %idToChit;
print p("These curated entries have", scalar(@uniqids), "distinct",
        a({ -href => $csURL . "&faa=1" }, "sequences") . "."),
  "\n";
FetchSeqs($blastdir, $blastdb, \@uniqids, $chitsfaaFile);

my %byCurated = (); # curated to list of hits; used to save the maximum score below

# Run ublast against the protein sequences
unless($isNuc) {
  open(FAA, ">", $seqFile) || die "Cannot write to $seqFile";
  while (my ($header, $seq) = each %seqs) {
    print FAA ">$header\n$seq\n";
  }
  close(FAA) || die "Error writing to $seqFile";

  my %parsed = (); # input sequence to list of hits

  print p("Running ublast with E &le; $maxEvalue\n");
  system("$usearch -ublast $chitsfaaFile -db $seqFile -evalue $maxEvalue -blast6out $ublastFile >& /dev/null") == 0
    || die "usearch failed: $!";
  unlink($seqFile);
  my $uhits = ParseUblast($ublastFile, \%seqs, \%idToChit);
  unlink($ublastFile);
  my %maxScore = (); # maximum score for each curated sequence

  foreach my $row (@$uhits) {
    push @{ $parsed{$row->{input}} }, $row;
    push @{ $byCurated{$row->{hit}} }, $row;
  }
  foreach my $input (keys %parsed) {
    my @rows = sort { $b->{score} <=> $a->{score} } @{ $parsed{$input} };
    $parsed{$input} = \@rows;
    $maxScore{$input} = $rows[0]{score};
  }
  print p("Found", scalar(keys %parsed), "relevant proteins in $genomeName, or try",
          a({-href => $URLnoq}, "another query"))."\n";
  my @inputs = sort { $maxScore{$b} <=> $maxScore{$a} } (keys %maxScore);
  foreach my $input (@inputs) {
    &PrintHits($input, $seqs{$input}, $parsed{$input}, 0); # 0 for proteins
  }
  print p(small("The hits are sorted by %identity * %coverage (highest first)"))
    if @inputs > 0;
}

# Search the six frame translation
my $fnafile; # genome sequence file (fasta nucleotide)
if ($assembly) {
  finish_page() if !exists $assembly->{fnafile}; # i.e. for UniProt proteomes
  $fnafile = $assembly->{fnafile};
  fail("Skipping 6-frame search -- no nucleotide sequences") unless $fnafile;
  open(my $fh, "<", $fnafile) || die "Cannot read $fnafile";
  my $state = {};
  my %ntlen = ();
  while (my ($header, $sequence) = ReadFastaEntry($fh,$state)) {
    fail("Duplicate sequence for $header in nucleotide sequence of $genomeName")
      if exists $ntlen{$header};
    $ntlen{$header} = length($sequence);
  }
  close($fh) || die "Error reading $fnafile";
  if (scalar(keys %ntlen) > $maxseqs) {
    print p("Not searching the 6-frame translation because this genome has too many scaffolds");
    finish_page();
  }
  my $tot = 0;
  foreach my $value (values %ntlen) { $tot += $value; }
  if ($tot > $maxNtLen) {
    print p("Not searching the 6-frame translation because this genome is too large");
    finish_page();
  }
} else { # uploaded
  if ($isNuc) {
    $fnafile = $seqFile;
    open(FNA, ">", $seqFile) || die "Cannot write to $seqFile";
    while (my ($header, $seq) = each %seqs) {
      print FNA ">$header\n$seq\n";
    }
    close(FNA) || die "Error writing to $seqFile";
  } else {
    # Do not search the 6-frame translation if uploaded a.a. sequences
    finish_page();
  }
}

my $xfile = "$fnafile.aa6";
if (! -e $xfile) {
  system("$usearch -fastx_findorfs $fnafile -aaout $xfile -orfstyle 7 -mincodons $minCodons >& /dev/null") == 0
    || die "usearch findorfs failed: $!";
}
unlink($seqFile) if $upfile;

# And read the 6-frame translation
my %seqsx = ();
$state = {};
open(my $fhx, "<", $xfile) || die "Cannot read $xfile";
while (my ($header, $sequence) = ReadFastaEntry($fhx,$state)) {
  $seqsx{$header} = $sequence;
}
close($fhx) || die "Error reading $xfile";

print p("Running ublast against the 6-frame translation.",
        "All reading frames of at least $minCodons codons are included."), "\n";
system("$usearch -ublast $chitsfaaFile -db $xfile -evalue $maxEvalue -blast6out $ublastFile >& /dev/null") == 0
  || die "usearch failed: $!";

my $uhits = ParseUblast($ublastFile, \%seqsx, \%idToChit);
unlink($ublastFile);
unlink($xfile) if $upfile;

my %maxCuratedScore = (); # curated => [maxscore, minbegin, maxend]
while (my ($curated, $hits) = each %byCurated) {
  my @scores = map $_->{score}, @$hits;
  my @beg = ();
  my @end = ();
  foreach my $hit (@$hits) {
    my ($cbeg,$cend) = split /:/, $hit->{hrange};
    push @beg, $cbeg;
    push @end, $cend;
  }
  $maxCuratedScore{$curated} = [max(@scores), min(@beg), max(@end)];
}
# Parse the 6-frame hits, ignoring any hits unless they are better than the best
# hit against the gene models
my %parsedx = (); # reading frame to list of hits
foreach my $row (@$uhits) {
  push @{ $parsedx{$row->{input}} }, $row;
}

# filter each list -- unless a hit is noticeably better than the best hit to an annotated protein,
# or it includes a region of the curated protein that was not seen before,
# it should be ignored. (The second case is so that for genes with the start codon
# or with frameshift errors, the hit to the other part of the gene is reported.)
#
# Also, if the best hit for a frame is masked in that way, mask all
# hits for that frame.
my $nWithHits = scalar(keys %parsedx);
foreach my $input (keys %parsedx) {
  my @rows = sort { $b->{score} <=> $a->{score} } @{ $parsedx{$input} };
  my @out = ();
  foreach my $i (0..(scalar(@rows)-1)) {
    my $row = $rows[$i];
    my $query = $row->{hit};
    my ($cbeg,$cend) = split /:/, $row->{hrange};
    my $maxC = exists $maxCuratedScore{$query} ? $maxCuratedScore{$query} : undef;
    $row->{maxCurated} = $maxC if defined $maxC;
    my ($maxS, $minBeg, $maxEnd) = @$maxC if defined $maxC;
    if (!defined $maxC
        || $row->{score} >= 1.1 * $maxS
        || $cbeg <= $minBeg - 5
        || $cend >= $maxEnd + 5) {
      push @out, $row;
    } else {
      last if $i == 0; # best hit must be useful or else suppress the reading frame entirely
    }
  }
  if (@out > 0) {
    $parsedx{$input} = \@out;
  } else {
    delete $parsedx{$input};
  }
}
my $nKept = scalar(keys %parsedx);
if ($nWithHits > 0) {
  my @found = ();
  if ($upfile && $isNuc) {
    push @found, "Found hits to $nWithHits reading frames.",
      "Or try", a({-href => $URLnoq}, "another query");
  } else {
    push @found, "Found hits to $nWithHits reading frames.";
    if ($nKept > 0) {
      push @found,
        qq{Except for $nKept reading frames, these were redundant with annotated proteins.
             These remaining reading frames may be pseudogenes, omissions in the genome annotation,
             or N-terminal extensions of annotated proteins.};
    } elsif (! $upfile) {
      push @found, "These were all redundant with annotated proteins.";
    }
  }
  print p(@found);
} else {
  print p("Did not find any hits to reading frames.");
  print p("Try", a({-href => $URLnoq}, "another query"))
    if $upfile && $isNuc;
}
my @inputsX = sort { $parsedx{$b}[0]{score} <=> $parsedx{$a}[0]{score} } (keys %parsedx);
foreach my $input (@inputsX) {
  &PrintHits($input, $seqsx{$input}, $parsedx{$input}, 1); # 1 for 6-frame translation
}

unlink($chitsfaaFile);
finish_page();


sub PrintHits($$$$) {
  my ($input, $seq, $hits, $sixframe) = @_;
  my $inputlink = $sixframe ? SixFrameLink($input, $hits) : ProteinLink($input);
  my $seqtype = $sixframe ? "reading frame" : "protein";
  my $pblink = small(a({ -href => "litSearch.cgi?query=>${input}%0A${seq}",
                         -title => "full PaperBLAST results for this $seqtype"},
                       "PaperBLAST"));

  my @show = ();
  my $iRow = 0;
  $nCollapseSet++;
  my $indentStyle = "margin: 0em; margin-left: 3em;";
  foreach my $row (@$hits) {
    $iRow++;
    my $chits = $row->{chits};
    my @descs = ();
    foreach my $chit (@$chits) {
      AddCuratedInfo($chit); # for the URL and showName fields
      my @showOrgWords = split / /, $chit->{organism};
      @showOrgWords = @showOrgWords[0..1] if @showOrgWords > 2;
      push @descs, a({-href => $chit->{URL}, -title => "from " . $chit->{db},
                      -onmousedown => loggerjs("curated", $chit->{subjectId}) },
                     $chit->{showName}) . ": " . $chit->{desc}
        . " " . small("from " . i(join(" ", @showOrgWords)));
    }
    my $percentcov = int($row->{coverage} * 100 + 0.5);
    my $clen = $chits->[0]{protein_length};
    my %trattr = ();
    $trattr{"-bgcolor"} = $iRow % 2 ? "#F2F2F2" : "#FCF3CF";
    if ($iRow > $maxHitsEach) {
      $trattr{"-class"} = "collapse${nCollapseSet}";
      $trattr{"-style"} = "display: none;"
    }
    my $showIdentity = int($row->{identity} + 0.5);
    my $seqlen = length($seq);
    if (exists $row->{maxCurated}) {
      my (undef, $minB, $maxE) = @{ $row->{maxCurated} };
      push @descs, small(i({-title => "Annotated proteins are similar to a.a. $minB:$maxE/$clen; this frame is similar to $row->{hrange}/$clen"},
                         "Also see hits to annotated proteins above"));
    }
    push @show, Tr(\%trattr,
                   td({-align => "left", -valign => "top"},
                      p({-style => $indentStyle}, join("<BR>", @descs))),
                   td({-align => "right", -valign => "top"},
                      small(a({ -href => "showAlign.cgi?" . join("&", "def1=$input", "seq1=$seq", "acc2=$row->{hit}"),
                                -title => "$row->{irange}/$seqlen aligns to $row->{hrange}/$clen of characterized protein" },
                              "${showIdentity}% id")
                            . ",<BR>${percentcov}% cov")));
  }
  push @show, Tr(td({ -colspan => 2 },
                    p({-style => $indentStyle},
                      a({ -href => "javascript:void(0);", -onclick => "tblexpander(this)" }, small("More...")))))
    if @show > $maxHitsEach;
  unshift @show, Tr(td({-align => "left", -valign => "top"}, $inputlink . "<BR>" . small("is similar to:")),
                    td({-align => "right", -valign => "top"}, $pblink));
  print p(table({-cellspacing => 0, -cellpadding => 2, -width => "100%" }, @show)), "\n";
}

sub ProteinLink($) {
    my ($input) = @_;
    return $input if !defined $gdb;
    my $inputlink = $input;
    if ($gdb eq "FitnessBrowser") {
      my @words = split / /, $input;
      my $locusId = shift @words;
      my $sysName = shift @words;
      my $desc = join(" ", @words);
      $inputlink = ($sysName || $locusId) . ": $desc";
      $inputlink = a({ -href => "http://fit.genomics.lbl.gov/cgi-bin/singleFit.cgi?orgId=${gid}&locusId=${locusId}" },
                     $inputlink);
    } elsif ($gdb eq "Microbesonline") {
      my @words = split / /, $input;
      my $locusId = shift @words;
      my $sysName = shift @words;
      my $desc = join(" ", @words);
      $inputlink = a({ -href => "http://www.microbesonline.org/cgi-bin/fetchLocus.cgi?locus=$locusId"},
                     $sysName || "VIMSS$locusId") . ": $desc";
    } elsif ($gdb eq "UniProt") {
      # either sp|accession|identifier description OS=organism GN=genename (and other attributes)
      # or tr|...
      # or for non reference proteomes, a different format (do not try to modify)
      if ($input =~ m/^[a-z][a-z][|]([A-Z90-9_]+)[|]([A-Z90-9_]+) (.*)$/) {
        my ($acc, $id, $desc) = ($1,$2,$3);
        my $gn = "";
        if ($desc =~ m/^(.*) OS=(.*)/) {
          $desc = $1;
          my $rest = $2;
          $gn = $1 if $rest =~ m/ GN=(\S+) /;
        }
        $inputlink = a({ -href => "http://www.uniprot.org/uniprot/" . $acc },
                       $acc) . " $desc";
        $inputlink .= " ($gn)" if $gn ne "";
      }
    } elsif ($gdb eq "NCBI") {
      # remove trailing organism descriptor
      $inputlink =~ s/\[[^\]]+\]$//;
      # change the initial protein id into a link
      if ($inputlink =~ m/^([A-Z0-9._]+) (.*)$/) {
        my ($acc, $desc) = ($1,$2);
        $desc =~ s/^MULTISPECIES: *//;
        my @ids = ();
        if (exists $assembly->{prot}{$acc}) {
          my $g = $assembly->{prot}{$acc};
          $desc = $g->{name} if $g->{name} ne "";
          push @ids, $g->{symbol} if $g->{symbol};
          if ($g->{locus_tag}) {
            my $geneurl = undef;
            my $title = undef;
            if ($g->{GeneID}) {
              $geneurl = "https://www.ncbi.nlm.nih.gov/gene/?term=" . $g->{GeneID};
              $title = "NCBI Gene";
            } elsif ($g->{genomic_accession} && $g->{start} && $g->{end}) {
              my $center = int(($g->{start} + $g->{end})/2);
              my ($left,$right) = ($center-5000,$center+5000);
              # The NCBI sequence viewer is smart enough to clip to valid regions
              $geneurl = "https://www.ncbi.nlm.nih.gov/nuccore/$g->{genomic_accession}/scaffold?report=graph&v=$left:$right";
              $title = "Genome Browser";
            }
            push @ids, defined $geneurl ? a( { -href => $geneurl, -title => $title }, $g->{locus_tag}) : $g->{locus_tag};
            push @ids, $assembly->{oldid}{$g->{locus_tag}} if exists $assembly->{oldid}{$g->{locus_tag}};
          }
        }
        push @ids, a({ -href => "https://www.ncbi.nlm.nih.gov/protein/$acc", -title => "NCBI Protein" }, $acc);
        $inputlink = join(" ", @ids) . ": " . $desc;
      }
    } elsif ($gdb eq "MicrobesOnline") {
      my @words = split / /, $input;
      my $locusId = shift @words;
      my $sysName = shift @words;
      my $desc = join(" ", @words);
      $inputlink = a({ -href => "http://www.microbesonline.org/cgi-bin/fetchLocus.cgi?locus=$locusId"},
                     $sysName || "VIMSS$locusId") . ": $desc";
    } elsif ($gdb eq "IMG") {
      my @words = split / /, $input;
      if (@words > 2 && $words[0] =~ m/^\d+$/ && $words[1] =~ m/^[A-Za-z][0-9a-zA-Z_]+$/ && $words[1] =~ m/[0-9_]/) {
        # The first word is the IMG gene # and the second word is a locus tag
        my $locusId = shift @words;
        my $sysName = shift @words;
        my $desc = join(" ", @words);
        $desc =~ s/ \[.*\]$//; # remove organism description at end
        $inputlink = a({ -href => "https://img.jgi.doe.gov/cgi-bin/w/main.cgi?section=GeneDetail&page=geneDetail&gene_oid="
                         . $locusId }, $sysName) . ": $desc";
      }
    }
    return $inputlink;
}

sub SixFrameLink($$) {
  my ($input, $hits) = @_;
  $input =~ m/^(.*)[|]([-+]\d):(\d+)-(\d+)[(](\d+)[)]$/
    || die "Cannot parse reading frame name $input";
  my ($scaffoldId, $frame, $begin, $end, $sclen) = ($1,$2,$3,$4,$5);
  die "Bad coordinates $begin:$end from reading frame name $input"
    unless $begin <= $end && $begin >= 1 && $end <= $sclen;
  # beg2 and end2 have beg2 > end2 if on - strand
  my ($beg2, $end2) = ($begin, $end);
  ($beg2,$end2) = ($end2,$beg2) if $frame < 0;
  $scaffoldId =~ s/ .*$// if $gdb && $gdb eq "IMG"; # remove extra fields from scaffold name
  my $show = "${begin}-${end} (frame $frame) on " . ($assembly ? $scaffoldId : "scaffold $scaffoldId");

  # If linking to a genome browser (either MicrobesOnline or the Fitness Browser), show
  # two objects, one for the entire frame and one for the expanse that has hits.
  my $minBeg = (sort {$a<=>$b} map { my ($b,$e) = split /:/, $_->{irange}; $b;  } @$hits)[0];
  my $maxEnd = (sort {$b<=>$a} map { my ($b,$e) = split /:/, $_->{irange}; $e;  } @$hits)[0];
  my $sign = $frame < 0 ? -1 : 1;
  my $beginUse = $sign > 0 ? $begin + ($minBeg-1) * 3 : $end - ($maxEnd-1) * 3;
  my $endUse = $sign > 0 ? $begin + ($maxEnd-1) * 3 : $end - ($minBeg-1) * 3;
  if ($gdb && $gdb eq "FitnessBrowser") {
    my $objspec1 = join(":",
                       "b", $begin, "e", $end,
                        "n", uri_escape("frame $frame"),
                       "s", $sign);
    my $objspec2 = join(":",
                        "b", $beginUse,
                        "e", $endUse,
                        "n", uri_escape("region with similarity"),
                        "s", $sign);
    my $URL = "http://fit.genomics.lbl.gov/cgi-bin/genomeBrowse.cgi?orgId=$gid&scaffoldId=$scaffoldId&object=$objspec1&object=$objspec2";
    $input = a({ -href => $URL, -title => "Fitness browser"}, $show);
  } elsif ($gdb && $gdb eq "MicrobesOnline") {
    my $object1 = join(":",
                       "f" . $beg2,
                       "t" . $end2,
                       "n" . uri_escape("Reading Frame $frame"),
                       "d" . uri_escape("From $beg2 to $end2 (frame $frame)"));
    # want begin > end if on - strand
    my $object2 = join(":",
                       "f" . $beginUse,
                       "t" . $endUse,
                       "n" . uri_escape("region with similarity"),
                       "d" . uri_escape("region with similarity ($minBeg/$maxEnd)"));
    my $URL = "http://www.microbesonline.org/cgi-bin/browser?"
                       . "mode=4;data=s${scaffoldId}:$object1:$object2";
    $input = a({ -href => $URL, -title => "MicrobesOnline genome browser" }, $show);
  } elsif ($gdb && $gdb eq "NCBI") {
    my $scaffold = $scaffoldId;
    $scaffold =~ s/ .*$//;
    $input = a({ -href => "https://www.ncbi.nlm.nih.gov/nuccore/$scaffold?report=graph"
                 . "&from=$begin&to=$end"
                 . "&mk=$beginUse:$endUse|hit_region|00008f",
                 -title => "NCBI's viewer" },
               $show);
  } elsif ($gdb && $gdb eq "IMG") {
    # Do not have the scaffold object id so cannot link to something like
    # https://img.jgi.doe.gov/cgi-bin/mer/main.cgi?section=ScaffoldGraph&page=scaffoldGraph&scaffold_oid=637000443&start_coord=4000&end_coord=29000
    $input = $show;
  }
  return $input;
}

sub ParseUblast($$$) {
  my ($file, $subjects, $idToChit) = @_;
  open (my $fh, "<", $file) || die "Cannot read $file";
  my @hits = ();
  while(<$fh>) {
    chomp;
    my ($query, $subject, $identity, $alen, $mm, $gap, $qbeg, $qend, $sbeg, $send, $eval, $bits)
      = split /\t/, $_;
    $query =~ s/^lcl[|]//;
    $query =~ s/ unnamed protein product$//;
    die "Unrecognized subject $subject" unless exists $subjects->{$subject};
    die "Unrecognized query $query" unless exists $idToChit->{$query};
    my $chits = $idToChit->{$query};
    my $clen = $chits->[0]{protein_length};
    my $coverage = ($qend - $qbeg + 1) / $clen;
    my $row = { input => $subject, hit => $query,
                identity => $identity, coverage => $coverage, score => $coverage * $identity,
                irange => "$sbeg:$send", hrange => "$qbeg:$qend",
                chits => $chits };
    push @hits, $row;
  }
  close($fh) || die "Error reading $file";
  return \@hits;
}

sub query_fields_html {
  my ($prefix) = @_;
  $prefix = defined $prefix ? "$prefix. " : "";
  return join("\n",
              p("${prefix}Enter a search term:",
                textfield(-name => 'query', -value => '', -size => 50, -maxlength => 200)),
              p({-style => "margin-left: 3em;" },
                small("Examples:", i("perchlorate"), "or", i("1.2.1.88"))),
              p({-style => "margin-left: 3em;" },
                checkbox(-name => "word", -label => "Match whole words only?"))
             );
}

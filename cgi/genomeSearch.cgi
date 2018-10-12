#!/usr/bin/perl -w
#######################################################
## genomeSearch.cgi
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
# Search -- set if the search button was pressed [used to handle failed uploads]

use strict;
use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use Time::HiRes qw{gettimeofday};
use DBI;
use LWP::Simple; # for get()
use URI::Escape; # for uri_escape()
use IO::Handle; # for autoflush
use lib "../lib";
use pbutils; # for ReadFastaEntry()
use FetchAssembly; # for FetchAssemblyInfo() etc.
use pbweb; # for TopDivHtml

sub start_page($);
sub GetMatchingAssemblies($$);
sub CacheAssembly($$$);

# page must be started already; reports any number of errors or warning
sub fail;
sub warning;
sub finish;

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
my $fastacmd = "$blastdir/fastacmd";
die "No such executable: $fastacmd" unless -x $fastacmd;
my $blastdb = "$base/uniq.faa";
die "No such file: $blastdb" unless -e $blastdb;

my $cgi=CGI->new;

my $gdb = $cgi->param('gdb');
my $gid = $cgi->param('gid');
my $gquery = $cgi->param('gquery');

my $upfile = $gid ? undef : $cgi->param('file');
my $query = $cgi->param('query');
my $word = $cgi->param('word');

my @gdbs = ("NCBI", "IMG", "UniProt", "MicrobesOnline", "FitnessBrowser");
my %gdb_labels1 = ("NCBI" => "NCBI assemblies",
                   "UniProt" => "UniProt proteomes",
                   "IMG" => "JGI/IMG genomes", "FitnessBrowser" => "Fitness Browser genomes");
my %gdb_labels = map { $_ => exists $gdb_labels1{$_} ? $gdb_labels1{$_} : "$_ genomes"} @gdbs;
die "Unknown genome database: $gdb\n"
  if $gdb && !exists $gdb_labels{$gdb};

&start_page("Curated BLAST for Genomes");

# A symbolic link to the Fitness Browser data directory is used (if it exists)
# to allow access to fitness browser genomes.
# That directory must include feba.db (sqlite3 database) and aaseqs (in fasta format)
my $fbdata = "../fbrowse_data"; # path relative to the cgi directory
fail("Cannot access Fitness Browser database")
  if $gid && $gid eq "FitnessBrowser" && ! -e $fbdata;

my $tmpDir = "../tmp";
my $procId = $$;
my $timestamp = int (gettimeofday() * 1000);
my $basefile = $tmpDir . "/" . $procId . $timestamp;

my $assembly; # the fetched (and cached) assembly
my $genomeName;
my $fh_up;

if ($gdb && $gquery) {
  print p("Searching", $gdb_labels{$gdb}, "for", "'" . $gquery . "'"), "\n";
  my @rows = GetMatchingAssemblies($gdb, $gquery);
  if (@rows > 0) {
    print p("Found",scalar(@rows),"assemblies"), "\n";
    print start_form,
      hidden(-name => 'gdb', -value => $gdb, -override => 1);
    foreach my $row (@rows) {
      my $checked = @rows == 1 ? "CHECKED" : "";
      print qq{<INPUT type="radio" NAME="gid" VALUE="$row->{gid}" $checked />},
        a({ -href => $row->{URL} }, $row->{genomeName} ),
        br();
    }
    print p("Enter a search term:",
            textfield(-name => 'query', -value => '', -size => 50, -maxlength => 200)),
          p(submit('Search selected genome')),
          end_form;;
  } else {
    print p("Sorry, no matching genomes were found.");
  }
  print p("Try", a({ -href => "genomeSearch.cgi?gdb=$gdb" }, "another query"));
  finish();
} elsif ($gdb && $gid && $query) {
  $assembly = CacheAssembly($gdb, $gid, "../tmp")
    || fail("Cannot fetch assembly $gid from database $gdb");
  $genomeName = $assembly->{genomeName};
  print p("Searching in", a({-href => $assembly->{URL} }, $genomeName)), "\n";
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
      query_fields(2),
      p(submit(-name => 'uploading', -value => 'Search')),
      end_form,
      p("Or", a({ -href => "genomeSearch.cgi"}, "search"), "for a genome");
  finish();
} else {
  warning("Please enter an organism name")
    if $cgi->param('findgenome');
  print
    p("Given a query term and a genome, find characterized proteins whose descriptions match the query,",
      "and then search the genome for homologs of those proteins."),
    start_form( -autocomplete => 'on', -name => 'input', -method => 'GET', -action => 'genomeSearch.cgi'),
    p("Genome database to search:", 
      popup_menu(-name => 'gdb', -values => \@gdbs, -labels => \%gdb_labels, -default => $gdbs[0])),
    p(textfield(-name => 'gquery', -value => '', -size => 50, -maxlength => 200)),
    p(submit(-name => "findgenome", -value => 'Find Genome')),
    end_form,
    p(br(),
      "Or",
      a({-href => "genomeSearch.cgi?doupload=1"}, "upload"),
      "a genome or proteome in fasta format");
  # Check $cgi->param('Search') in case of failed uploads ??
  finish();
}

fail("No genome fetched") unless defined $genomeName;
die "No query\n" unless $query;

# This should really be in pbweb
my %sourceToURL = ( "SwissProt" => "http://www.uniprot.org/uniprot/",
                    "SwissProt/TReMBL" => "http://www.uniprot.org/uniprot/",
                    "BRENDA" => "http://www.brenda-enzymes.org/sequences.php?AC=",
                    "MicrobesOnline" => "http://www.microbesonline.org/cgi-bin/fetchLocus.cgi?locus=",
                    "RefSeq" => "http://www.ncbi.nlm.nih.gov/protein/",
                    "metacyc" => "https://metacyc.org/gene?orgid=META&id=",
                    "ecocyc" => "https://ecocyc.org/gene?orgid=ECOLI&id=",
                    "CAZy" => "http://www.cazy.org/search?page=recherche&lang=en&tag=4&recherche="
                  );

my %seqs = (); # The main sequence set
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

my %seqs = ();
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
  finish();
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
my $maxhits = 1000;
my $limit = $maxhits + 1;
my $chits = $dbh->selectall_arrayref("SELECT * FROM CuratedGene WHERE desc LIKE ? LIMIT $limit",
                                     { Slice => {} }, "%" . $query . "%");
my $URLnoq = $upfile ? "genomeSearch.cgi?doupload=1" : "genomeSearch.cgi?gdb=$gdb&gid=" . uri_escape($gid);
if (@$chits > $maxhits) {
  print p(qq{Sorry, too many curated entries match the query '$query'. Please try},
          a({ -href => $URLnoq }, "another query").".");
  print end_html;
  exit(0);
}
if (@$chits == 0) {
  print p(qq{None of the curated entries in PaperBLAST's database match '$query'. Please try},
          a({ -href => $URLnoq }, "another query") . ".");
  print end_html;
  exit(0);
}
if ($word) {
  # filter for whole-word hits
  my $quoted = quotemeta($query); # this will quote % as well
  $quoted =~ s/\\%/\\b.*\\b/g; # turn % into a separation of words; note quoting of \\ so that it appears in the new string

  my @keep = grep { $_->{desc} =~ m/\b$quoted\b/i } @$chits;
  if (@keep == 0) {
    print p(qq{None of the curated entries in PaperBLAST's database match '$query' as complete words. Please try},
            a({ -href => $URLnoq }, "another query") . ".");
    print end_html;
    exit(0);
  }
  $chits = \@keep;
}

my $wordstatement = $word ? " as complete word(s)" : "";
print p("Found", scalar(@$chits), qq{curated entries in PaperBLAST's database that match '$query'${wordstatement}.\n});

my $listFile = "$basefile.list";
my $chitsfaaFile = "$basefile.chits.faa";
my $seqFile = "$basefile.seq";
my $ublastFile = "$basefile.ublast";

# Make the input file for fastacmd
my %idToChit = (); # sequence identifier to curated gene hit(s)
foreach my $hit (@$chits) {
  my $seqid = $hit->{db} . "::" . $hit->{protId};
  my $uniqRef = $dbh->selectcol_arrayref("SELECT sequence_id FROM SeqToDuplicate WHERE duplicate_id = ?",
                                         {}, $seqid);
  $seqid = $uniqRef->[0] if scalar(@$uniqRef) > 0;
  push @{ $idToChit{$seqid} }, $hit;
}
print p("These curated entries have", scalar(keys %idToChit), "distinct sequences.\n");
open(LIST, ">", $listFile) || die "Cannot write to $listFile";
foreach my $id (sort keys %idToChit) {
  print LIST "$id\n";
}
close(LIST) || die "Error writing to $listFile";
system("$fastacmd -i $listFile -d $blastdb -p T > $chitsfaaFile") == 0
  || die "fastacmd failed: $!";
unlink($listFile);

my %byCurated = (); # curated to list of hits; used to save the maximum score below

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
}

# And search the six frame translation
my $ntfile; # genome sequence file (fasta nucleotide)
if ($assembly) {
  $ntfile = $assembly->{ntFile};
  finish() unless $ntfile;
  if ($ntfile) {
    open(my $fh, "<", $ntfile) || die "Cannot read $ntfile";
    my $state = {};
    my %ntlen = ();
    while (my ($header, $sequence) = ReadFastaEntry($fh,$state)) {
      fail("Duplicate sequence for $header in nucleotide sequence of $genomeName")
        if exists $ntlen{$header};
      $ntlen{$header} = length($sequence);
    }
    close($fh) || die "Error reading $ntfile";
    if (scalar(keys %ntlen) > $maxseqs) {
      print p("Not searching the 6-frame translation because this genome has too many scaffolds");
      finish();
    }
    my $tot = 0;
    foreach my $value (values %ntlen) { $tot += $value; }
    if ($tot > $maxNtLen) {
      print p("Not searching the 6-frame translation because this genome is too large");
      finish();
    }
  }
} else { # uploaded
  if ($isNuc) {
    $ntfile = $seqFile;
    open(FNA, ">", $seqFile) || die "Cannot write to $seqFile";
    while (my ($header, $seq) = each %seqs) {
      print FNA ">$header\n$seq\n";
    }
    close(FNA) || die "Error writing to $seqFile";
  } else {
    # Do not search the 6-frame translation if uploaded a.a. sequences
    finish();
  }
}

my $xfile = "$ntfile.aa6";
if (! -e $xfile) {
  system("$usearch -fastx_findorfs $ntfile -aaout $xfile -orfstyle 7 -mincodons $minCodons >& /dev/null") == 0
    || die "usearch findorfs failed: $!";
}
unlink($seqFile) if $upfile;

# And read the 6-frame translation
my %seqsx = ();
my $state = {};
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

my %maxCuratedScore = ();
while (my ($curated, $hits) = each %byCurated) {
  my $maxscore = 0;
  foreach my $hit (@$hits) {
    $maxscore = $hit->{score} if $hit->{score} > $maxscore;
  }
  $maxCuratedScore{$curated} = $maxscore;
}
# Parse the 6-frame hits, ignoring any hits unless they are better than the best
# hit against the gene models
my %parsedx = (); # reading frame to list of hits
foreach my $row (@$uhits) {
  push @{ $parsedx{$row->{input}} }, $row;
}
# filter each list -- unless a hit is noticeably better than the best hit to an annotated protein,
# it should be ignored. Also, if the best hit for a frame is masked in that way, mask all
# hits for that frame.
my $nWithHits = scalar(keys %parsedx);
foreach my $input (keys %parsedx) {
  my @rows = sort { $b->{score} <=> $a->{score} } @{ $parsedx{$input} };
  my @out = ();
  foreach my $i (0..(scalar(@rows)-1)) {
    my $row = $rows[$i];
    my $query = $row->{hit};
    if (!exists $maxCuratedScore{$query} || $row->{score} >= 1.1 * $maxCuratedScore{$query}) {
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
  my @found = ("Found hits to $nWithHits reading frames.");
  if ($nKept > 0) {
    push @found,
      qq{Except for $nKept reading frames, these were redundant with annotated proteins.
             These remaining reading frames may be pseudogenes, omissions in the genome annotation,
             or N-terminal extensions of annotated proteins.};
  } elsif (! $upfile) {
    push @found, "These were all redundant with annotated proteins.";
  }
  print p(@found);
} else {
  print p("Did not find any hits to reading frames.");
}
print p("Or try", a({-href => $URLnoq}, "another query"))."\n"
  if $upfile && $isNuc;
my @inputsX = sort { $parsedx{$b}[0]{score} <=> $parsedx{$a}[0]{score} } (keys %parsedx);
foreach my $input (@inputsX) {
  &PrintHits($input, $seqsx{$input}, $parsedx{$input}, 1); # 1 for 6-frame translation
}
unlink($chitsfaaFile);
finish();

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
      # build a URL for this sequence, and, show an identifier
      my $db = $chit->{db};
      my $protId = $chit->{protId};
      my $URL = "";
      if ($db eq "reanno") {
        my ($orgId, $locusId) = split /:/, $protId;
        $URL = "http://fit.genomics.lbl.gov/cgi-bin/singleFit.cgi?orgId=$orgId&locusId=$locusId";
      }  elsif ($db eq "REBASE") {
        $URL = "http://rebase.neb.com/rebase/enz/${protId}.html";
      } elsif ($db eq "CharProtDB") {
        $URL = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3245046/"; # the CharProtDB paper
      } else {
        die "No URL for $db" unless exists $sourceToURL{ $db };
        $URL = $sourceToURL{$db} . $protId if $sourceToURL{$db};
      }
      my $showId = $protId;
      $showId =~ s/^.*://;
      $showId = $chit->{id2} if $db eq "reanno" && $chit->{id2};
      $showId = "VIMSS$showId" if $showId =~ m/^\d+/ && $db eq "reanno";
      $showId = "$showId / $chit->{id2}" if $db eq "SwissProt" && $chit->{id2};

      $showId = "$showId / $chit->{name}" if $chit->{name};
      $showId = $chit->{name} if $db eq "ecocyc" && $chit->{name};
      my @showOrgWords = split / /, $chit->{organism};
      @showOrgWords = @showOrgWords[0..1] if @showOrgWords > 2;
      push @descs, a({-href => $URL, -title => "from $db"}, $showId) . ": " . $chit->{desc}
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
                      a({ -href => "javascript:void(0);", -onclick => "expander(this,$nCollapseSet)" }, small("More...")))))
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
      my $orginfo = {}; #XXX
      my @words = split / /, $input;
      my $locusId = shift @words;
      my $sysName = shift @words;
      my $desc = join(" ", @words);
      $inputlink = ($sysName || $locusId) . ": $desc";
      $inputlink = a({ -href => "http://fit.genomics.lbl.gov/cgi-bin/singleFit.cgi?gid=${gid}&locusId=${locusId}" },
                     $inputlink)
        if exists $orginfo->{$gid};
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
        my %assemblyProt = (); #XXX
        my %assemblyOldLocusTag = (); #XXX
        if (exists $assemblyProt{$acc}) {
          my $g = $assemblyProt{$acc};
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
            push @ids, $assemblyOldLocusTag{$g->{locus_tag}} if exists $assemblyOldLocusTag{$g->{locus_tag}};
          }
        }
        push @ids, a({ -href => "https://www.ncbi.nlm.nih.gov/protein/$acc", -title => "NCBI Protein" }, $acc);
        $inputlink = join(" ", @ids) . ": " . $desc;
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

sub start_page($) {
  my ($title) = @_;
#XXX move to pbwb
my $style = <<END
.autocomplete {
  /*the container must be positioned relative:*/
  position: relative;
  display: inline-block;
}

.autocomplete-items {
  position: absolute;
  border: 1px solid #d4d4d4;
  border-bottom: none;
  border-top: none;
  z-index: 99;
  /*position the autocomplete items to be the same width as the container:*/
  top: 100%;
  left: 0;
  right: 0;
}
.autocomplete-items div {
  padding: 10px;
  cursor: pointer;
  background-color: #fff;
  border-bottom: 1px solid #d4d4d4;
}
.autocomplete-items div:hover {
  /*when hovering an item:*/
  background-color: #e9e9e9;
}
.autocomplete-active {
  /*when navigating through the items using the arrow keys:*/
  background-color: DodgerBlue !important;
  color: #ffffff;
}
END
;

print
  header(-charset => 'utf-8'),
  start_html(-head => Link({-rel => "shortcut icon", -href => "../static/favicon.ico"}),
             -style => { -code => $style },
             -script => [{ -type => "text/javascript", -src => "../static/autocomplete_uniprot.js" }],
             -title => $title),
  TopDivHtml(),
  h2($title);
  autoflush STDOUT 1; # show preliminary results
  # Autoexpander script
  print <<END
<SCRIPT>
function expander(o,n) {
  var x = document.getElementsByClassName("collapse"+n);
  console.log("Expander " + n + " count " + x.length);
  var i;
  for (i = 0; i < x.length; i++) {
    x[i].style.display = "table-row";
  }
  o.parentElement.style.display = "none";
}
</SCRIPT>
END
    ;
  print "\n";
}

sub warning {
  print p({ -style => "color: red;" }, @_), "\n";
}

sub fail {
  warning(@_);
  finish();
}

sub finish {
print <<END
<P>
<small>
<center>by <A HREF="http://morgannprice.org/">Morgan Price</A>,
<A HREF="http://genomics.lbl.gov/">Arkin group</A><BR>
Lawrence Berkeley National Laboratory
</center>
</small>
</P>
END
;
print end_html;
exit(0);
}

# Fetch matching genomes, and return list of hashes
# Each hash includes gid, genomeName, and URL
sub GetMatchingAssemblies($$) {
  my ($gdb,$gquery) = @_;
  return unless $gdb && $gquery;
  if ($gdb eq "NCBI") {
    my @hits = FetchNCBIInfo($gquery);
    foreach my $hit (@hits) {
      $hit->{gid} = $hit->{id};
      $hit->{genomeName} = $hit->{org};
      $hit->{URL} = "https://www.ncbi.nlm.nih.gov/assembly/" . $hit->{id};
    }
    return @hits;
  }
  # else
  fail("Database $gdb is not supported");
}

sub query_fields($) {
  my $prefix = "";
  $prefix = "2. " if @_;
  return join("\n",
              p("${prefix}Enter a search term:",
                textfield(-name => 'query', -value => '', -size => 50, -maxlength => 200)),
              p({-style => "margin-left: 5em;" },
                checkbox(-name => "word", -checked => 0, -label => "Match whole words only?"))
             );
}

sub CacheAssembly($$$) {
  my ($gdb, $gid, $dir) = @_;
  return undef unless $gdb && $gid;
  die "Not a directory: $dir\n" unless -d $dir;
  if ($gdb eq "NCBI") {
    my @hits = FetchNCBIInfo($gid);
    fail("Do not recognize NCBI assembly $gid")
      unless @hits;
    my $assembly = $hits[0];
    # note redundancy with code above
    $assembly->{gid} = $assembly->{id};
    $assembly->{genomeName} = $assembly->{org};
    $assembly->{URL} = "https://www.ncbi.nlm.nih.gov/assembly/" . $assembly->{id};
    my $faafile = "$dir/refseq_" . $assembly->{gid} . ".faa";
    my $ntfile = "$dir/refseq_" . $assembly->{gid} . ".fna";
    my $featurefile = "$tmpDir/refseq_" . $assembly->{id} . ".features.tab";
    unless (-e $faafile) {
      print "<P>Fetching protein fasta file for $assembly->{gid}\n";
      fail("Sorry, failed to fetch the protein fasta file for this assembly ($!). This assembly might not have any predicted proteins.")
        unless FetchNCBIFaa($assembly, $faafile);
      fail("Sorry, failed to fetch the nucleotide assembly for this assembly: $!")
        unless FetchNCBIFna($assembly, $ntfile);
      fail("Sorry, failed to fetch the feature file for this assembly: $!")
        unless &FetchNCBIFeatureFile($assembly, $featurefile);
    }
    my $features = ParseNCBIFeatureFile($featurefile);
    $assembly->{prot} = {};
    foreach my $row (@$features) {
      next unless $row->{class} eq "with_protein";
      $assembly->{prot}{$row->{product_accession}} = $row;
      $assembly->{prot}{$row->{"non-redundant_refseq"}} = $row;
    }
    foreach my $row (@$features) {
      if ($row->{class} eq "protein_coding"
          && $row->{locus_tag}
          && $row->{attributes} =~ m/old_locus_tag=([A-Z0-9_]+)/) {
        $assembly->{oldlocustag}{$row->{locus_tag}} = $1;
      }
    }
    $assembly->{faafile} = $faafile;
    $assembly->{ntfile} = $ntfile;
    return $assembly;
  }
  # else
  fail("Database $gdb is not supported");
}

#!/usr/bin/perl -w
#######################################################
## genomeSearch.cgi
##
## Copyright (c) 2018 University of California
##
## Authors: Morgan Price
#######################################################
#
# Optional CGI garameters:
# query -- what to search for
# file -- an uploaded file with protein sequences in FASTA format
# Search -- set if the search button was pressed

use strict;
use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use Time::HiRes qw{gettimeofday};
use DBI;
use IO::Handle; # for autoflush
use lib "../lib";
use pbutils; # for ReadFastaEntry()

# maximum size of posted data, in bytes
my $maxMB = 100;
$CGI::POST_MAX = $maxMB*1024*1024;
my $maxseqsK = 100;
my $maxseqs = $maxseqsK * 1000;
my $maxseqsComma = "$maxseqsK,000";
my $maxEval = 0.01;

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
my $upfile = $cgi->param('file');
my $query = $cgi->param('query');

# This should really be in pbutils
my %sourceToURL = ( "SwissProt" => "http://www.uniprot.org/uniprot/",
                    "SwissProt/TReMBL" => "http://www.uniprot.org/uniprot/",
                    "BRENDA" => "http://www.uniprot.org/uniprot/",
                    "MicrobesOnline" => "http://www.microbesonline.org/cgi-bin/fetchLocus.cgi?locus=",
                    "RefSeq" => "http://www.ncbi.nlm.nih.gov/protein/",
                    "metacyc" => "https://metacyc.org/gene?orgid=META&id=",
                    "ecocyc" => "https://ecocyc.org/gene?orgid=ECOLI&id=",
                    "CharProtDB" => "",
                    "CAZy" => "http://www.cazy.org/search?page=recherche&lang=en&tag=4&recherche="
                  );

my $title = "Curated BLAST for Genomes";
print
  header(-charset => 'utf-8'),
  start_html(-head => Link({-rel => "shortcut icon", -href => "../static/favicon.ico"}),
             -title => $title),
  h2($title);

if ($upfile && $query) {
  # Sufficient information to execute a query
  # First validate the fasta input
  my $up = $cgi->upload('file');
  die "Cannot upload $upfile" unless $up;
  my $fh = $up->handle;
  my %seqs = (); # header (with > removed) to sequence
  my %seqlen = ();
  my $state = {};
  while (my ($header, $sequence) = ReadFastaEntry($fh,$state)) {
    die "Duplicate sequence for $header\n" if exists $seqs{$header};
    die ". found in sequence for $header\n" if $sequence =~ m/[.]/;
    die "- found in sequence for $header\n" if $sequence =~ m/-/;
    $sequence =~ s/[*]//g;
    die "Invalid/empty sequence for $header\n" if $sequence eq "";
    $seqs{$header} = $sequence;
    $seqlen{$header} = length( $sequence );
  }
  die "Too many sequences: limit $maxseqs\n" if scalar(keys %seqs) > $maxseqs;
  die "No sequences in uploaded file\n" if scalar(keys %seqs) == 0;
  print p("Found", scalar(keys %seqs), "sequences in $upfile.");
  autoflush STDOUT 1; # show preliminary results
  print "\n";

  # Find relevant sequences in CuratedGene
  my $maxhits = 1000;
  my $limit = $maxhits + 1;
  my $chits = $dbh->selectall_arrayref("SELECT * FROM CuratedGene WHERE desc LIKE ? LIMIT $limit",
                                       { Slice => {} }, "%" . $query . "%");
  if (@$chits > $maxhits) {
    print p(qq{Sorry, too many curated entries match the query '$query'. Please choose a more specific query.});
    print end_html;
    exit(0);
  }
  if (@$chits == 0) {
    print p(qq{None of the curated entries in PaperBLAST's database match '$query'. Please try another query.});
    print end_html;
    exit(0);
  }
  print p("Found", scalar(@$chits), qq{curated entries in PaperBLAST's database that match '$query'.\n});

  my $tmpDir = "../tmp";
  my $procId = $$;
  my $timestamp = int (gettimeofday() * 1000);
  my $filename = $procId . $timestamp;
  my $listFile = "$tmpDir/$filename.list";
  my $chitsfaaFile = "$tmpDir/$filename.chits.faa";
  my $genomefaaFile = "$tmpDir/$filename.genome.faa";
  my $ublastFile = "$tmpDir/$filename.ublast";

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

  open(FAA, ">", $genomefaaFile) || die "Cannot write to $genomefaaFile";
  while (my ($header, $seq) = each %seqs) {
    print FAA ">$header\n$seq\n";
  }
  close(FAA) || die "Error writing to $genomefaaFile";

  print p("Running ublast with E &le; $maxEval\n");
  system("$usearch -ublast $chitsfaaFile -db $genomefaaFile -evalue $maxEval -blast6out $ublastFile >& /dev/null") == 0
           || die "usearch failed: $!";
  unlink($genomefaaFile);
  unlink($chitsfaaFile);

  my $uhits = [];
  open(HITS, "<", $ublastFile) || die "Cannot read $ublastFile";
  while(<HITS>) {
    chomp;
    my @F = split /\t/, $_;
    push @$uhits, \@F;
  }
  close(HITS) || die "Error reading $ublastFile";
  unlink($ublastFile);
  print p("Found", scalar(@$uhits), "hits\n");

  # Parse the hits. Query is from the curated hits; subject is from the input genome; output
  # will be sorted by subject, with subjects sorted by their top hit (as %identity * %coverage)
  my %parsed = (); # input sequence to list of hits
  foreach my $uhit (@$uhits) {
    my ($query, $subject, $identity, $alen, $mm, $gap, $qbeg, $qend, $sbeg, $send, $eval, $bits) = @$uhit;
    $query =~ s/^lcl[|]//;
    $query =~ s/ unnamed protein product$//;
    die "Unrecognized subject $subject" unless exists $seqs{$subject};
    die "Unrecognized query $query" unless exists $idToChit{$query};
    my $chits = $idToChit{$query};
    my $clen = $chits->[0]{protein_length};
    my $coverage = ($qend - $qbeg + 1) / $clen;
    my $row = { input => $subject, hit => $query,
                identity => $identity, coverage => $coverage, score => $coverage * $identity,
                irange => "$sbeg:$send", hrange => "$qbeg:$qend",
                chits => $chits };
    push @{ $parsed{$subject} }, $row;
  }
  my %maxScore = ();
  while (my ($input, $rows) = each %parsed) {
    my @rows = sort { $b->{score} <=> $a->{score} } @$rows;
    $parsed{$input} = \@rows;
    $maxScore{$input} = $rows[0]{score};
  }
  my @inputs = sort { $maxScore{$b} <=> $maxScore{$a} } (keys %maxScore);
  foreach my $input (@inputs) {
    my $header = "Hits for " .
      a({ -href => "litSearch.cgi?query=>${input}%0A$seqs{$input}", -title => "run PaperBLAST on $input" }, $input);
    my @show = ();
    my $rows = $parsed{$input};
    foreach my $row (@$rows) {
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
        } else {
          die "No URL for $db" unless exists $sourceToURL{ $db };
          $URL = $sourceToURL{$db} . $protId if $sourceToURL{$db};
        }
        my $showId = $protId;
        $showId =~ s/^.*://;
        $showId = "VIMSS$showId" if $showId =~ m/^\d+/ && $db eq "reanno";
        push @descs, a({-href => $URL, -title => "from $db"}, $showId) . ": " . $chit->{desc} . " from " . i($chit->{organism});
      }
      my $percentcov = int($row->{coverage} * 100 + 0.5);
      my $clen = $chits->[0]{protein_length};
      push @show, li(join("; ", @descs),
                     a({ -href => "showAlign.cgi?" . join("&", "def1=$input", "seq1=$seqs{$input}", "acc2=$row->{hit}"),
                         -title => "$row->{irange}/$seqlen{$input} versus $row->{hrange}/$clen" },
                       "($row->{identity}% identity, ${percentcov}% coverage)"));
    }
    print p({ -style => "margin-top: 1em; margin-bottom: 0em;" },
            $header, qq{<UL style="margin-top: 0em; margin-bottom: 0em;">}, @show, "</UL>")."\n";;
  }
} else {
  # Show the query form.
  # Note that failed uploads reach here because CGI sets upfile to undef
  if ($cgi->param('Search')) {
    print p("Cannot search without an input genome.") unless $upfile;
    print p("Please enter a query.") unless $query;
  }
  print
    p("Given a query term, find characterized proteins that are relevant and find their homologs in a genome."),
    start_form( -name => 'input', -method => 'POST', -action => 'genomeSearch.cgi'),
    p("1. Query:", textfield(-name => "query", -value => '', -size => 50, -maxlength => 200)),
    p({-style => "margin-left: 5em;" }, "use % as a wild card that matches any substring"),
    p("2. Genome: upload proteins in FASTA format", filefield(-name=>'file', -size=>50)),
    p({-style => "margin-left: 5em;" },"up to $maxseqsComma amino acid sequences or $maxMB MB"),
    p(submit('Search')),
    end_form;
}

print end_html;
exit(0);

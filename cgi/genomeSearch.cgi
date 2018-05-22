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
# alternative ways to specify which genome to search in:
#	orgId -- an organism identifier in the fitness browser
#	mogenome -- an genone name in MicrobesOnline
#	file -- an uploaded file with protein sequences in FASTA format
#
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
my $orgId = $cgi->param('orgId');
my $mogenome = $cgi->param('mogenome');
my $upfile = $cgi->param('file');
my $query = $cgi->param('query');

# A symbolic link to the Fitness Browser data directory is used (if it exists)
# to allow quick access to fitness browser genomes.
# That directory must include feba.db (sqlite3 database) and aaseqs (in fasta format)
my $fbdata = "../fbrowse_data"; # path relative to the cgi directory
my $fbdbh;
my $orginfo;
if (-e $fbdata) {
  $fbdbh = DBI->connect("dbi:SQLite:dbname=$fbdata/feba.db","","",{ RaiseError => 1 }) || die $DBI::errstr;
  die "Missing aaseqs file for fitness browser: $fbdata/aaseqs\n"
    unless -e "$fbdata/aaseqs";
  $orginfo = $fbdbh->selectall_hashref("SELECT * FROM Organism", "orgId");
  foreach my $hash (values %$orginfo) {
    $hash->{genome} = join(" ", $hash->{genus}, $hash->{species}, $hash->{strain});
  }
}

# The MicrobesOnline taxonomy table (must include the taxonomyId and name fields)
# This file is optional. If it exists, uses MicrobesOnline's public mysql server
# to fetch genomes for a taxonomyId
my $mofile = "../static/Taxonomy.mo";
my %moTax = (); # name to taxonomyId
if (-e $mofile) {
  my @tax = &ReadTable($mofile, ["taxonomyId","name"]);
  %moTax = map { $_->{name} => $_->{taxonomyId} } @tax;
}

my $tmpDir = "../tmp";
my $procId = $$;
my $timestamp = int (gettimeofday() * 1000);
my $basefile = $tmpDir . "/" . $procId . $timestamp;

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
  autoflush STDOUT 1; # show preliminary results
print "\n";

if (($mogenome || $orgId || $upfile) && $query) {
  # Sufficient information to execute a query
  my $fh; # read the fasta file
  if ($orgId) {
    die "Invalid organism id $orgId\n" unless exists $orginfo->{$orgId};
    print p("Loading the genome of $orginfo->{$orgId}{genome} from the",
            a({-href => "http://fit.genomics.lbl.gov/cgi-bin/orgAll.cgi"}, "Fitness Browser"));
    my $cachedfile = "$tmpDir/fbrowse_${orgId}.faa";
    my $nWrite = 0;
    if (! -e $cachedfile) {
      my $tmpfile = "$basefile.tmp";
      open(TMP, ">", $tmpfile) || die "Cannot write to $tmpfile";
      my $genes = $fbdbh->selectall_hashref("SELECT * from Gene WHERE orgId = ?", "locusId", {}, $orgId);
      open(my $aafh, "<", "$fbdata/aaseqs") || die "Cannot read $fbdata/aaeqs";
      my $state = {};
      while (my ($header, $sequence) = ReadFastaEntry($aafh,$state)) {
        my ($org2, $locusId) = split /:/, $header;
        die $header unless defined $locusId;
        next unless $org2 eq $orgId;
        my $gene = $genes->{$locusId}
          || die "Unrecognized locusId $locusId";
        print TMP ">$locusId $gene->{sysName} $gene->{desc}\n$sequence\n";
        $nWrite++;
      }
      close($aafh) || die "Error reading $fbdata/aaseqs";
      close(TMP) || die "Error writing to $tmpfile";
      rename($tmpfile, $cachedfile) || die "Rename $tmpfile to $cachedfile failed";
    }
    open($fh, "<", $cachedfile) || die "Cannot read $cachedfile";
  } elsif ($mogenome) {
    my $taxId = $moTax{$mogenome}
      || die "Invalid MicrobesOnline genome name $mogenome";
    print p("Loading the genome of $mogenome from", a({-href => "http://www.microbesonline.org/" }, "MicrobesOnline")), "\n";
    my $cachedfile = "$tmpDir/mogenome_${taxId}.faa";
    if (! -e $cachedfile) {
      my $mo_dbh = DBI->connect('DBI:mysql:genomics:pub.microbesonline.org', "guest", "guest")
        || die $DBI::errstr;
      my $genes = $mo_dbh->selectall_arrayref(qq{ SELECT locusId, sequence
                                                  FROM Locus JOIN Scaffold USING (scaffoldId)
                                                  JOIN AASeq USING (locusId, version)
                                                  WHERE taxonomyId = ? AND isActive=1 AND priority=1 AND type=1; },
                                              { Slice => {} }, $taxId);
      my $desc = $mo_dbh->selectall_hashref(qq{ SELECT locusId, description
                                                    FROM Locus JOIN Scaffold USING (scaffoldId)
                                                    JOIN Description USING (locusId,version)
                                                    WHERE taxonomyId = ? AND isActive=1 AND priority=1 AND Locus.type=1 },
                                                "locusId", {}, $taxId);
      my $sysNames = $mo_dbh->selectall_hashref(qq{ SELECT locusId, name
                                                    FROM Locus JOIN Scaffold USING (scaffoldId)
                                                    JOIN Synonym USING (locusId, version)
                                                    WHERE taxonomyId = ? AND isActive=1 AND priority=1 AND Locus.type=1 AND Synonym.type = 1 },
                                                "locusId", {}, $taxId);
      my $tmpfile = "$basefile.tmp";
      open(TMP, ">", $tmpfile) || die "Cannot write to $tmpfile";
      foreach my $gene (@$genes) {
        my $locusId = $gene->{locusId};
        my $sysName = $sysNames->{$locusId}{name} || "";
        my $desc = $desc->{$locusId}{description} || "";
        print TMP ">$locusId $sysName $desc\n$gene->{sequence}\n";
      }
      close(TMP) || die "Error writing to $tmpfile";
      rename($tmpfile, $cachedfile) || die "Rename $tmpfile to $cachedfile failed";
    }
    open($fh, "<", $cachedfile) || die "Cannot read $cachedfile";
  } elsif ($upfile) {
    my $up = $cgi->upload('file');
    die "Cannot upload $upfile" unless $up;
    $fh = $up->handle;
  } else {
    die "Unreachable";
  }

  # First validate the fasta input
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
  close($fh) || die "Error reading genome file";
  die "Too many sequences: limit $maxseqs\n" if scalar(keys %seqs) > $maxseqs;
  die "No sequences in genome file\n" if scalar(keys %seqs) == 0;
  print p("Found", scalar(keys %seqs), "sequences in $upfile.\n") if $upfile;

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

  my $listFile = "$basefile.list";
  my $chitsfaaFile = "$basefile.chits.faa";
  my $genomefaaFile = "$basefile.genome.faa";
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
    my $inputlink = $input;
    if ($orgId) {
      my @words = split / /, $input;
      my $locusId = shift @words;
      my $sysName = shift @words;
      my $desc = join(" ", @words);
      $inputlink = a({ -href => "http://fit.genomics.lbl.gov/cgi-bin/singleFit.cgi?orgId=${orgId}&locusId=${locusId}" },
                     $sysName) . ": $desc";
    } elsif ($mogenome) {
      my @words = split / /, $input;
      my $locusId = shift @words;
      my $sysName = shift @words;
      my $desc = join(" ", @words);
      $inputlink = a({ -href => "http://www.microbesonline.org/cgi-bin/fetchLocus.cgi?locus=$locusId"},
                     $sysName || "VIMSS$locusId") . ": $desc";
    }
    my $pblink = a({ -href => "litSearch.cgi?query=>${input}%0A$seqs{$input}" }, "PaperBLAST");
    my $header = join(" ", $inputlink, small($pblink));

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
  my @genomeSelectors = ();
  if ($fbdbh) {
    my @orginfo = sort { $a->{genome} cmp $b->{genome} } values(%$orginfo);
    my %orgLabels = ("" => "From the Fitness Browser:");
    my @orgOptions = ("");
    foreach my $hash (@orginfo) {
      my $orgId = $hash->{orgId};
      push @orgOptions, $orgId;
      $orgLabels{$orgId} = $hash->{genome};
    }
    push @genomeSelectors, p(popup_menu( -name => 'orgId', -values => \@orgOptions, -labels  => \%orgLabels,
                                       -default => ''));
  }
  if (keys(%moTax) > 0) {
    my @tax = sort keys %moTax;
    my %taxLabels = map { $_ => $_ } @tax;
    $taxLabels{""} = "From MicrobesOnline:";
    my @taxOptions = sort keys %taxLabels;
    push @genomeSelectors, p(popup_menu( -name => 'mogenome', -values => \@taxOptions, -labels => \%taxLabels,
                                         -default => ''));
  }
  push @genomeSelectors,
    p("or upload proteins in FASTA format", filefield(-name=>'file', -size=>50)),
    p({-style => "margin-left: 5em;" },"up to $maxseqsComma amino acid sequences or $maxMB MB");
  print
    p("Given a query term, find characterized proteins that are relevant and find their homologs in a genome."),
    start_form( -name => 'input', -method => 'POST', -action => 'genomeSearch.cgi'),
    p("1. Query:", textfield(-name => "query", -value => '', -size => 50, -maxlength => 200)),
    p({-style => "margin-left: 5em;" }, "use % as a wild card that matches any substring"),
    p("2. Select genome:"),
    @genomeSelectors,
    p(submit('Search'), reset()),
    end_form;
}

print end_html;
exit(0);

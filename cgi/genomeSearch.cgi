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
# word -- if non-empty, report whole word matches only
# alternative ways to specify which genome to search in:
#	orgId -- an organism identifier in the fitness browser
#	mogenome -- an genome name in MicrobesOnline
#	uniprotname -- a uniprot proteome's name
#	file -- an uploaded file with protein sequences in FASTA format
#
# Search -- set if the search button was pressed

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

# maximum size of posted data, in bytes
my $maxMB = 100;
$CGI::POST_MAX = $maxMB*1024*1024;
my $maxseqsK = 100;
my $maxseqs = $maxseqsK * 1000;
my $maxNtLen = 30 * 1000 * 1000;
my $maxseqsComma = "$maxseqsK,000";
my $maxEval = 0.01;
my $maxHitsEach = 3; # additional hits are hidden until the expander is clicked
my $nCollapseSet = 0; # keeping track fo which expander goes with what

my $minCodons = 30; # for reading frames

# from a defline to a brief description and a link to the source
sub ProteinLink($);
sub SixFrameLink($);

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
my $orgId = $cgi->param('orgId');
my $mogenome = $orgId ? undef : $cgi->param('mogenome');
my $uniprotname = $orgId || $mogenome ? undef : $cgi->param('uniprotname');
my $upfile = $orgId || $mogenome || $uniprotname ? undef : $cgi->param('file');
my $query = $cgi->param('query');
my $word = $cgi->param('word');

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
                    "BRENDA" => "http://www.brenda-enzymes.org/sequences.php?AC=",
                    "MicrobesOnline" => "http://www.microbesonline.org/cgi-bin/fetchLocus.cgi?locus=",
                    "RefSeq" => "http://www.ncbi.nlm.nih.gov/protein/",
                    "metacyc" => "https://metacyc.org/gene?orgid=META&id=",
                    "ecocyc" => "https://ecocyc.org/gene?orgid=ECOLI&id=",
                    "CAZy" => "http://www.cazy.org/search?page=recherche&lang=en&tag=4&recherche="
                  );

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

my $title = "Curated BLAST for Genomes";
print
  header(-charset => 'utf-8'),
  start_html(-head => Link({-rel => "shortcut icon", -href => "../static/favicon.ico"}),
             -style => { -code => $style },
             -script => [{ -type => "text/javascript", -src => "../static/autocomplete_uniprot.js" }],
             -title => $title),
    qq{<div style="background-color: #40C0CB; display: block; position: absolute; top: 0px; left: -1px; width: 100%; padding: 0.25em; z-index: 400;"><H2 style="margin: 0em;"><A HREF="litSearch.cgi" style="color: gold; font-family: 'Montserrat', sans-serif; font-style:italic; text-shadow: 1px 1px 1px #000000; text-decoration: none;">PaperBLAST &ndash; <small>Find papers about a protein or its homologs</small></A></H2></div><P style="margin: 0em;">&nbsp;</P>\n},
  h2($title);
  autoflush STDOUT 1; # show preliminary resul
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

my $hasGenome = ($mogenome || $orgId || $uniprotname || $upfile);
if ($hasGenome && $query) {
  # Sufficient information to execute a query
  # First read the protein fasta
  my $fh; # reads the fasta file
  my $genomeName;
  my $cachedfile = "$tmpDir/fbrowse_${orgId}.faa";
  if ($orgId) {
    # Allow use of orgId that is not in the database if the cached file has already been set up.
    # This is a way to support use of "local" genomes (although they are still not available
    # in the selector).
    if (!exists $orginfo->{$orgId} && ! -e $cachedfile) {
      print p("Sorry, organism nickname $orgId is not known. Try",
              a({-href => "genomeSearch.cgi?query=$query"}, "selecting it from UniProt or uploading it"),
              "instead."),
            end_html;
      exit(0);
    }
    $genomeName = exists $orginfo->{$orgId} ? $orginfo->{$orgId}{genome} : "local genome $orgId";
    print p("Loading the proteome of $genomeName from the",
            a({-href => "http://fit.genomics.lbl.gov/cgi-bin/orgAll.cgi"}, "Fitness Browser"))
      if exists $orginfo->{$orgId};
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
    $genomeName = $mogenome;
    print p("Loading the proteome of $mogenome from", a({-href => "http://www.microbesonline.org/" }, "MicrobesOnline")), "\n";
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
  } elsif ($uniprotname) {
    $genomeName = $uniprotname;
    my $dbfile = "../static/uniprot_proteomes.db";
    die "No such file: $dbfile\n" unless -e $dbfile;
    my $dbh = DBI->connect("dbi:SQLite:dbname=$dbfile","","",{ RaiseError => 1 })
      || die $DBI::errstr;
    my $uniprotList = $dbh->selectcol_arrayref("SELECT upid FROM Proteome WHERE name = ?", {}, $uniprotname);
    unless (scalar(@$uniprotList) == 1) {
      print p("Proteome name $uniprotname is not valid.");
      exit(0);
    }
    my ($upid) = @$uniprotList;
    print p("Loading the proteome of $uniprotname from UniProt",
            a({-href => "https://www.uniprot.org/proteomes/$upid"}, $upid)), "\n";
    my $cachedfile = "$tmpDir/uniprot_${upid}.faa";
    unless (-e $cachedfile) {
      # First try getting from https://www.uniprot.org/uniprot/?query=proteome:UP000002886&format=fasta
      # This is a bit slow, so I considered using links like
      # ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Archaea/UP000000242_399549.fasta.gz
      # but those are also surprisingly slow, and would need to figure out which section to look in.
      my $refURL = "https://www.uniprot.org/uniprot/?query=proteome:${upid}&format=fasta";
      my $uniparcURL = "https://www.uniprot.org/uniparc/?query=proteome:${upid}&format=fasta";
      my $faa = get($refURL);
      $faa = get($uniparcURL) if $faa eq "";
      if ($faa eq "") {
        print p("Proteome $uniprotname seems to be empty, see",
                a({href => $uniparcURL}, "here"));
        exit(0);
      }
      my $tmpfile = "$basefile.tmp";
      open($fh, ">", $tmpfile) || die "Cannot write $tmpfile";
      print $fh $faa;
      close($fh) || die "Error writing to $tmpfile";
      rename($tmpfile, $cachedfile) || die "Rename $tmpfile to $cachedfile failed";
    }
    open($fh, "<", $cachedfile) || die "Cannot read $cachedfile";
  } elsif ($upfile) {
    $genomeName = "uploaded file";
    my $up = $cgi->upload('file');
    die "Cannot upload $upfile" unless $up;
    $fh = $up->handle;
  } else {
    die "Unreachable";
  }

  # Validate the fasta input
  my %seqs = (); # header (with > removed) to sequence
  my $state = {};
  while (my ($header, $sequence) = ReadFastaEntry($fh,$state)) {
    die "Duplicate sequence for $header\n" if exists $seqs{$header};
    die ". found in sequence for $header\n" if $sequence =~ m/[.]/;
    die "- found in sequence for $header\n" if $sequence =~ m/-/;
    $sequence =~ s/[*]//g;
    die "Invalid/empty sequence for $header\n" if $sequence eq "";
    $seqs{$header} = $sequence;
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
  if ($word) {
    # filter for whole-word hits
    my $quoted = quotemeta($query); # this will quote % as well
    $quoted =~ s/\\%/\\b.*\\b/g; # turn % into a separation of words; note quoting of \\ so that it appears in the new string

    my @keep = grep { $_->{desc} =~ m/\b$quoted\b/i } @$chits;
    if (@keep == 0) {
      print p(qq{None of the curated entries in PaperBLAST's database match '$query' as complete words. Please try another query.});
      print end_html;
      exit(0);
    }
    $chits = \@keep;
  }

  my $wordstatement = $word ? " as complete word(s)" : "";
  print p("Found", scalar(@$chits), qq{curated entries in PaperBLAST's database that match '$query'${wordstatement}.\n});

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

  my $uhits = ParseUblast($ublastFile, \%seqs, \%idToChit);
  unlink($ublastFile);

  my $URLq = "genomeSearch.cgi";
  if ($orgId) {
    $URLq .= "?orgId=$orgId";
  } elsif ($mogenome) {
    $URLq .= "?mogenome=$mogenome";
  } elsif ($uniprotname) {
    $URLq .= "?uniprotname=" . uri_escape($uniprotname);
  }

  my %parsed = (); # input sequence to list of hits
  my %byCurated = (); # curated to list of hits
  foreach my $row (@$uhits) {
    push @{ $parsed{$row->{input}} }, $row;
    push @{ $byCurated{$row->{hit}} }, $row;
  }
  my %maxScore = ();
  foreach my $input (keys %parsed) {
    my @rows = sort { $b->{score} <=> $a->{score} } @{ $parsed{$input} };
    $parsed{$input} = \@rows;
    $maxScore{$input} = $rows[0]{score};
  }
  print p("Found", scalar(keys %parsed), "relevant proteins in $genomeName, or try",
          a({-href => $URLq}, "another query"))."\n";
  my @inputs = sort { $maxScore{$b} <=> $maxScore{$a} } (keys %maxScore);
  foreach my $input (@inputs) {
    &PrintHits($input, $seqs{$input}, $parsed{$input}, 0); # 0 for proteins
  }

  # And search the six frame translation
  my $ntfile; # genome sequence file (fasta nucleotide)
  my $sc = []; # list of [scaffoldId, sequence], if ntfile does not already exist
  if ($orgId) {
    $ntfile = "$tmpDir/fbrowse_${orgId}.fna";
    $sc = $fbdbh->selectall_arrayref("SELECT scaffoldId, sequence FROM ScaffoldSeq WHERE orgId = ?",
                                     {}, $orgId)
      unless -e $ntfile;
  } elsif ($mogenome) {
    my $taxId = $moTax{$mogenome}
      || die "Invalid MicrobesOnline genome name $mogenome";
    $genomeName = $mogenome;
    print p("Loading the genome of $mogenome from", a({-href => "http://www.microbesonline.org/" }, "MicrobesOnline")), "\n";
    $ntfile = "$tmpDir/mogenome_${taxId}.fna";
    if (! -e $ntfile) {
      my $mo_dbh = DBI->connect('DBI:mysql:genomics:pub.microbesonline.org', "guest", "guest")
        || die $DBI::errstr;
      $sc = $mo_dbh->selectall_arrayref(qq{ SELECT scaffoldId, sequence
                                               FROM Scaffold JOIN ScaffoldSeq USING (scaffoldId)
                                               WHERE taxonomyId = ? AND isActive = 1 },
                                           {}, $taxId);
      die "Cannot fetch genome sequence for $taxId from MicrobesOnline\n"
        unless @$sc > 0;
    }
  }

  # create ntfile if necessary
  if (@$sc > 0) {
    my $totlen = 0;
    foreach my $row (@$sc) {
      $totlen += length($row->[1]);
    }
    if ($totlen > $maxNtLen) {
      print p("Skipping the 6-frame translation because the genome is too large");
    } else {
      my $tmpfile = "$basefile.fna";
      open(FNA, ">", $tmpfile) || die "Cannot write to $tmpfile";
      foreach my $row (@$sc) {
        my ($scaffoldId, $seq) = @$row;
        print FNA ">${scaffoldId}\n$seq\n";
      }
      close(FNA) || die "Error writing to $tmpfile";
      rename($tmpfile, $ntfile) || die "Rename $tmpfile to $ntfile failed";
    }
  }

  # If we successfully created a genome file then
  if (-e $ntfile) {
    my $xfile = "$ntfile.aa6";
    if (! -e $xfile) {
      system("$usearch -fastx_findorfs $ntfile -aaout $xfile -orfstyle 7 -mincodons $minCodons >& /dev/null") == 0
        || die "usearch findorfs failed: $!";
    }
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
    system("$usearch -ublast $chitsfaaFile -db $xfile -evalue $maxEval -blast6out $ublastFile >& /dev/null") == 0
      || die "usearch failed: $!";

    my $uhits = ParseUblast($ublastFile, \%seqsx, \%idToChit);
    unlink($ublastFile);

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
      } else {
        push @found, "These were all redundant with annotated proteins.";
      }
      print p(@found),"\n";
    } else {
      print p("Did not find any hits to reading frames.");
    }
    my @inputsX = sort { $parsedx{$b}[0]{score} <=> $parsedx{$a}[0]{score} } (keys %parsedx);
    foreach my $input (@inputsX) {
      &PrintHits($input, $seqsx{$input}, $parsedx{$input}, 1); # 1 for 6-frame translation
    }
  }
  unlink($chitsfaaFile);


} else {
  # Show the query form.
  # Note that failed uploads reach here because CGI sets upfile to undef
  if ($cgi->param('Search')) {
    print p("Cannot search without an input genome.") if $query;
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
  my $uniprot_default = $uniprotname ? qq{ value="$uniprotname" } : qq{ placeholder="Genome name" };
  my $uniprot_selector = <<END
<div class="autocomplete" style="width:100%;">
From UniProt:<BR>
<input id="uniprotname" type="text" name="uniprotname" $uniprot_default style="width:90%;" >
</div>
<SCRIPT>
autocomplete(document.getElementById("uniprotname"));
</SCRIPT>
END
;
  push @genomeSelectors,
    $uniprot_selector,
    p("Or upload proteins in FASTA format (up to $maxseqsComma amino acid sequences or $maxMB MB)",
      br(),
      filefield(-name=>'file', -size=>50));
  print
    p("Given a query term, find characterized proteins whose descriptions match the query. Then, search a genome for homologs of those proteins."),

    start_form( -autocomplete => 'off', -name => 'input', -method => 'POST', -action => 'genomeSearch.cgi'),
    p(b("1. Enter a query:"), textfield(-name => "query", -value => '', -size => 50, -maxlength => 200)),
    p({-style => "margin-left: 5em;" }, "use % as a wild card that matches any substring"),
    p({-style => "margin-left: 5em;" }, checkbox(-name => "word", -checked => 0, -label => "Match whole words only?")),
    p(b("2. Select a genome:")),
    @genomeSelectors,
    p(submit('Search'), reset()),
    end_form,
}

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

sub PrintHits($$$$) {
  my ($input, $seq, $hits, $sixframe) = @_;
  my $inputlink = $sixframe ? SixFrameLink($input) : ProteinLink($input);
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
    my $inputlink = $input;
    if ($orgId) {
      my @words = split / /, $input;
      my $locusId = shift @words;
      my $sysName = shift @words;
      my $desc = join(" ", @words);
      $inputlink = ($sysName || $locusId) . ": $desc";
      $inputlink = a({ -href => "http://fit.genomics.lbl.gov/cgi-bin/singleFit.cgi?orgId=${orgId}&locusId=${locusId}" },
                     $inputlink)
        if exists $orginfo->{$orgId};
    } elsif ($mogenome) {
      my @words = split / /, $input;
      my $locusId = shift @words;
      my $sysName = shift @words;
      my $desc = join(" ", @words);
      $inputlink = a({ -href => "http://www.microbesonline.org/cgi-bin/fetchLocus.cgi?locus=$locusId"},
                     $sysName || "VIMSS$locusId") . ": $desc";
    } elsif ($uniprotname) {
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
    }
    return $inputlink;
}

sub SixFrameLink($) {
  my ($input) = @_;
  $input =~ m/^(.*)[|]([-+]\d):(\d+)-(\d+)[(](\d+)[)]$/
    || die "Cannot parse reading frame name $input";
  my ($scaffoldId, $frame, $begin, $end, $sclen) = ($1,$2,$3,$4,$5);
  die "Bad coordinates $begin:$end from reading frame name $input"
    unless $begin <= $end && $begin >= 1 && $end <= $sclen;
  # beg2 and end2 have beg2 > end2 if on - strand
  my ($beg2, $end2) = ($begin, $end);
  ($beg2,$end2) = ($end2,$beg2) if $frame < 0;
  my $show = "${begin}-${end} (frame $frame) on scaffold $scaffoldId";
  if ($orgId) {
    my $objspec = join(":",
                       "b", $begin, "e", $end, "n", uri_escape("frame $frame"),
                       "s", $frame < 0 ? -1 : 1);
    my $URL = "http://fit.genomics.lbl.gov/cgi-bin/genomeBrowse.cgi?orgId=$orgId&scaffoldId=$scaffoldId&object=$objspec";
    $input = a({ -href => $URL, -title => "Fitness browser"}, $show);
  } elsif ($mogenome) {
    my $URL = "http://www.microbesonline.org/cgi-bin/browser?"
      . "mode=4;data=s${scaffoldId}:nReading Frame $frame:f${beg2}t${end2}:dFrom $beg2 to $end2 (frame $frame)";
    $input = a({ -href => $URL, -title => "MicrobesOnline genome browser" }, $show);
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

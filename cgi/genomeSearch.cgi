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
my $maxseqsComma = "$maxseqsK,000";
my $maxEval = 0.01;

my $maxHitsEach = 3;

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
my $uniprotname = $cgi->param('uniprotname');
my $upfile = $cgi->param('file');
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
                    "BRENDA" => "http://www.uniprot.org/uniprot/",
                    "MicrobesOnline" => "http://www.microbesonline.org/cgi-bin/fetchLocus.cgi?locus=",
                    "RefSeq" => "http://www.ncbi.nlm.nih.gov/protein/",
                    "metacyc" => "https://metacyc.org/gene?orgid=META&id=",
                    "ecocyc" => "https://ecocyc.org/gene?orgid=ECOLI&id=",
                    "CharProtDB" => "",
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
    x[i].style.display = "list-item";
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
  } elsif ($uniprotname) {
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

  my $URLq = "genomeSearch.cgi";
  if ($orgId) {
    $URLq .= "?orgId=$orgId";
  } elsif ($mogenome) {
    $URLq .= "?mogenome=$mogenome";
  } elsif ($uniprotname) {
    $URLq .= "?uniprotname=" . uri_escape($uniprotname);
  }
  print p("Found", scalar(@$uhits), "hits, or try",
          a({-href => $URLq}, "another query"))."\n";

  # Parse the hits. Query is from the curated hits; subject is from the input genome; output
  # will be sorted by subject, with subjects sorted by their top hit (as %identity * %coverage)
  my %parsed = (); # input sequence to list of hits
  my $nCollapseSet = 0;
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
    } elsif ($uniprotname) {
      # either sp|accession|identifier description OS=organism GN=genenname (and other attributes)
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
    if (@show > $maxHitsEach) {
      my @first = ();
      my @rest = @show;
      while(@first < $maxHitsEach) {
        push @first, (pop @rest);
      }
      @show = @first;
      $nCollapseSet++;
      push @show, li(a({ -href => "javascript:void(0);", -onclick => "expander(this,$nCollapseSet)" }, "More..."));
      foreach my $show (@rest) {
        $show =~ s/^<LI>/<LI class="collapse${nCollapseSet}" style="display: none;">/i;
        push @show, $show;
      }
    }
    print p({ -style => "margin-top: 1em; margin-bottom: 0em;" },
            $header, qq{<UL style="margin-top: 0em; margin-bottom: 0em;">}, @show, "</UL>")."\n";;
  }
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
    p("1. Query:", textfield(-name => "query", -value => '', -size => 50, -maxlength => 200)),
    p({-style => "margin-left: 5em;" }, "use % as a wild card that matches any substring"),
    p({-style => "margin-left: 5em;" }, checkbox(-name => "word", -checked => 0, -label => "Match whole words only?")),
    p("2. Select genome:"),
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

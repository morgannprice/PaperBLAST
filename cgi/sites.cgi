#!/usr/bin/perl -w
use strict;

# sites.cgi -- compare a protein sequence to known functional sites
#
# Optional CGI parameters:
# query -- this should be the protein sequence in FASTA or UniProt
# hmmId -- an hmm identifier from PFam or TIGRFam, or of the form hex.MD5 (created when uploading an hmm)
# hmmfile -- an uploaded hmm to search
# format -- set to "tsv" to download a tab-delimited table
#
# If query is not specified, shows a query box

use strict;
use CGI qw(:standard Vars start_ul);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use Time::HiRes qw{gettimeofday};
use DBI;
use LWP::Simple qw{get};
use HTML::Entities;
use IO::Handle; # for autoflush
use lib "../lib";
use pbutils qw{NewerThan};
use pbweb qw{start_page finish_page GetMotd loggerjs parseSequenceQuery sequenceToHeader HmmToFile analysisLinks};
use Bio::SearchIO;
use URI::Escape;

sub FormatAlnString($$$$$$$);
sub FormatAlnSites($$$$$$);

sub AddMouseOver($$$);

my %charToLong = ("A" => "Ala", "C" => "Cys", "D" => "Asp", "E" => "Glu",
                  "F" => "Phe", "G" => "Gly", "H" => "His", "I" => "Ile",
                  "K" => "Lys", "L" => "Leu", "M" => "Met", "N" => "Asn",
                  "P" => "Pro", "Q" => "Gln", "R" => "Arg", "S" => "Ser",
                  "T" => "Thr", "V" => "Val", "W" => "Trp", "Y" => "Tyr",
                  "O" => "Pyrrolysine", "U" => "Selenocysteine");
# These indicate which amino acid goes with which style, i.e. see .aa1 in pbweb::start_page
my %charSets = ();
if (defined param('color') && param('color') eq "clustal") {
  %charSets = ("AILMFWV" => 1,
               "KR" => 2,
               "ED" => 3,
               "NQST" => 4,
               "C" => 5,
               "G" => 6,
               "P" => 7,
               "HY" => 8);
} else {
  %charSets = map { $_ => $_ } split //, "ACDEFGHIKLMNPQRSTVWY";
  $charSets{"-"} = "Gap";
}

my %charToSet = ();

my $borderCol = "#EEEEEE"; # also in pb.js alnHighlight()

my $tmpDir = "../tmp";
my $blastall = "../bin/blast/blastall";
my $nCPU = 4;
my $base = "../data";
my $blastdb = "$base/hassites.faa";
my $sqldb = "$base/litsearch.db";
my $fbdata = "../fbrowse_data"; # path relative to the cgi directory

my $procId = $$;
my $timestamp = int (gettimeofday() * 1000);
my $filename = $procId . $timestamp;
my $seqFile = "$tmpDir/$filename.fasta";

my $docstring = <<END
SitesBLAST's database includes
(1) <A HREF="https://web.expasy.org/docs/swiss-prot_guideline.html">SwissProt</A>
entries with experimentally-supported functional features;
and (2) protein structures with bound ligands, from the
<A HREF="https://zhanglab.ccmb.med.umich.edu/BioLiP">BioLip</A> database.
END
;

my $maxMB = 10; # maximum size of uploaded HMM
$CGI::POST_MAX = $maxMB*1024*1024;

my $cgi=CGI->new;
my $query = $cgi->param('query') || "";
my $hmmSpec = $cgi->param('hmmId') || "";
$hmmSpec =~ s/^\s+//;
$hmmSpec =~ s/\s+$//;
my $hmmUp = $cgi->upload('hmmfile');
my $maxHits = 20;
my $maxE = 0.001;

my $format = $cgi->param('format') || "";
$format = "" unless $format eq "tsv";

my $dbh = DBI->connect("dbi:SQLite:dbname=$sqldb","","",{ RaiseError => 1 }) || die $DBI::errstr;
autoflush STDOUT 1; # show preliminary results

my ($header, $seq, $hmmFile, $queryLen);
if ($hmmSpec) {
  die "Invalid hmmSpec" unless $hmmSpec =~ m/^[a-zA-Z0-9_.-]+$/;
  $hmmFile = HmmToFile($hmmSpec);
  # Compute query length
  open(my $fh, "<", $hmmFile) || die "Cannot read $hmmFile";
  while(my $line = <$fh>) {
    $queryLen = $1 if $line =~ m/^LENG\s+(\d+)$/;
  }
  close($fh) || die "Error reading $hmmFile";
  die "No LENG line in $hmmFile\n" unless $queryLen;
} elsif ($hmmUp) {
  # Check the uploaded file, save it, set up an MD5-based hmmid, and redirect
  my $fh = $hmmUp->handle;
  my @lines = <$fh>;
  my @errors = ();
  push @errors, "Uploaded file does not start with HMMER3"
    unless ($lines[0] =~ m/^HMMER3/i);
  push @errors, "Uploaded file does not end with a // line"
    unless (@lines > 0 && $lines[-1] =~ m!^//\r?\n?!);
  my @endlines = grep m!^//!, @lines;
  push @errors, "Uploaded file has more than one record-ending // line"
    unless @endlines == 1;
  my @namelines = grep m/^NAME /, @lines;
  push @errors, "Uploaded file has no NAME field"
    unless @namelines > 0;
  push @errors, "Uploaded file has more than one NAME"
    unless @namelines <= 1;
  if (@errors > 0) {
    print header,
      start_html(-title => "HMM Upload failed"),
      h2("Not a valid HMM file"),
      p(join(". ", @errors)),
      a({ -href => "sites.cgi"}, "Try another search");
    finish_page();
  }
  my $name = $namelines[0];
  chomp $name;
  $name =~ s/^NAME +//;

  my $hex = Digest::MD5::md5_hex(@lines);
  my $hmmId = "hex.$hex";
  my $file = "../tmp/$hmmId.hmm";
  unless (-e $file) {
    open(my $sfh, ">", $file) || die "Cannot write to $file";
    print $sfh @lines;
    close($sfh) || die "Error writing to $file";
  }
  print $cgi->redirect("sites.cgi?hmmId=$hmmId");
  exit(0);
}
# else it is a sequence, but, must call parseSequenceQuery after page is started

if ($format ne "tsv") {
  my $title = "SitesBLAST";
  if ($hmmSpec) {
    my $showId = $hmmSpec;
    if (defined $hmmFile && $hmmSpec =~ m/^hex[.]/) {
      $showId = "Uploaded HMM $hmmSpec";
    } else {
      my $hmmName = `egrep '^NAME' $hmmFile`;
      $hmmName =~ s/[\r\n].*//;
      $hmmName =~ s/^NAME +//;
      $showId = "$hmmSpec ($hmmName)" unless $hmmSpec eq $hmmName;
    }
    $title = "SitesHMM for $showId";
  }
  start_page('title' => $title,
             'banner' => 'SitesBLAST &ndash; <small>Find functional sites</small>',
             'bannerURL' => 'sites.cgi');
  print <<END
<SCRIPT src="../static/pb.js"></SCRIPT>
END
    ;
}

if (!defined $hmmFile) {
  ($header, $seq) = parseSequenceQuery(-query => $query,
                                       -dbh => $dbh,
                                       -blastdb => $blastdb,
                                       -fbdata => $fbdata);
  $queryLen = length($seq);
}

if ($format eq "tsv") {
  die "No query found\n" if !defined $seq && !defined $hmmFile;
  print "Content-Type:text/tab-separated-values\n";
  my $title = $hmmSpec ? "SitesHMM" : "SitesBLAST";
  print "Content-Disposition: attachment; filename=$title.tsv\n\n";
}

if ($format ne "tsv" && !defined $hmmFile && !defined $seq) {
  print
    GetMotd(),
    p("Given a protein sequence, SitesBLAST finds homologs that have known functional residues and",
      "shows whether the functional residues are conserved.",
      small("(" . a({-href => "sites.cgi?query=VIMSS590795"}, "example" ) . ")")),
    p($docstring),
    start_form( -name => 'input', -method => 'GET', -action => 'sites.cgi'),
    p(br(),
      b("Enter a protein sequence in FASTA format, or an identifier from UniProt, RefSeq, or MicrobesOnline"),
      br(),
      textarea( -name  => 'query', -value => '', -cols  => 70, -rows  => 10 )),
    p("Or an HMM accession:", textfield(-name => "hmmId", -value => "", -size => 12),
      br(),
      small("like",
            a({-href => "sites.cgi?hmmId=TIGR00110"}, "TIGR00110"),
            "or",
            a({-href => "sites.cgi?hmmId=PF01817"}, "PF01817"))),
    p(submit('Search'), reset()),
    end_form,
    start_form(-name => 'upload', -method => 'POST', -action => 'sites.cgi'),
    p("Or upload an HMM:",
      filefield(-name => 'hmmfile', -size => 50),
      submit('Go')),
    p("Or try", a({-href => "treeSites.cgi", -title =>"Sites on a Tree: view functional residues in an alignment"},
                  "Sites on a Tree"),
     "or", a({-href => "litSearch.cgi", -title => "PaperBLAST: Find papers about a protein or its homologs"},
             "PaperBLAST"));
  finish_page();
  exit(0);
}

my @hits;
if (defined $hmmFile) {
  $header = $hmmSpec =~ m/^hex/ ? "uploaded HMM" : $hmmSpec;
  unless ($format eq "tsv") {
    print
      p("Comparing $header to proteins with known functional sites using HMMer"),
      p("Or try",
        a({ -href => "hmmSearch.cgi?hmmId=$hmmSpec"}, "Family Search vs. Papers")),
      "\n";
  }
  my $resultsFile = "$tmpDir/$hmmSpec.hmmer.sites";
  unless (-e $resultsFile
          && NewerThan($resultsFile, $hmmFile)
          && NewerThan($resultsFile,  $blastdb)) {
    my $hmmsearch = "../bin/hmmsearch";
    die "No such executable: $hmmsearch" unless -x $hmmsearch;
    # If it has no trusted cutoff line, try running it with -E 0.01
    my $tmpResultsFile = "$resultsFile.$$.tmp";
    if (system("$hmmsearch", "--cut_tc", "-o", $tmpResultsFile, $hmmFile, $blastdb) != 0) {
      print p("Rerunning HMMer with -E 0.01 instead of trying to use the trusted cutoff")."\n";
      system($hmmsearch, "-E", 0.01, "-o", $tmpResultsFile, $hmmFile, $blastdb) == 0
        || die "Error runniung hmmsearch: $!";
    }
    rename($tmpResultsFile, $resultsFile)
      || die "Error renaming to $resultsFile";
  }
  # The hmmer parser sometimes outputs lines like
  # "Missed this line: [No individual domains that satisfy reporting thresholds (although complete target did)]"
  # So, temporarily redirect STDOUT
  open(my $old_stdout, ">&STDOUT") || die "Cannot duplicate STDOUT: $!";
  open(STDOUT, ">", "/dev/null") || die "Cannot redirect STDOUT: $!";
  my $searchio = Bio::SearchIO->new(-format => 'hmmer', -file => $resultsFile)
    || die "Failed to read $resultsFile";
  while (my $result = $searchio->next_result) {
    while (my $hit = $result->next_hit) {
      # In some cases, a hit has no high-scoring domains, but meets the overall threshold; skip those
      my @hsps = $hit->hsps;
      if (@hsps > 0) {
        push @hits, $hit;
        last if @hits >= $maxHits;
      }
    }
  }
  open(STDOUT, ">&", $old_stdout) || die "Cannot duplicate filehandle: $!";
} else {
  my $hasDef = $header ne "";
  $header = sequenceToHeader($seq) unless $hasDef;
  $query = ">$header\n$seq\n";

  my $query2 = $query; $query2 =~ s/[|]/./g;

  unless ($format eq "tsv") {
    my @links = analysisLinks('desc' => $header, 'seq' => $seq,
                              'skip' => { 'SitesBLAST' => 1 },
                              'fbLoad' => 1);
    print
      qq{<DIV style="float:right; padding-left: 10px; width:25%; background-color: #EEEEEE; font-size: 95%;">},
      map(p({-style => "margin-top: 0.5em; margin-bottom: 0.5em;"}, $_),
          "Other sequence analysis tools:", @links),
            qq{</DIV>},
      p("Comparing $header to proteins with known functional sites using BLASTp with E &le; $maxE."),
      p("Or try",
        a({-href => "treeSites.cgi?query=".uri_escape($query),
           -title => "Sites on a Tree: view functional residues in an alignment"},
          "Sites on a Tree")),
        "\n";
  }
  open(my $fhFaa, ">", $seqFile) || die "Cannot write to $seqFile\n";
  print $fhFaa ">$header\n$seq\n";
  close($fhFaa) || die "Error writing to $seqFile\n";
  die "No such executable: $blastall\n" unless -x $blastall;
  # m S means mask complex sequences for lookup but not for alignment
  system("$blastall -F 'm S' -p blastp -i $seqFile -d $blastdb -e $maxE -a $nCPU -o $seqFile.out -b $maxHits -v $maxHits -a $nCPU > /dev/null 2>&1") == 0
    || die "$blastall failed: $!\n";
  unlink($seqFile);
  my $searchio = Bio::SearchIO->new(-format => 'blast', -file => "$seqFile.out")
    || die "Failed to read $seqFile.out\n";
  while (my $result = $searchio->next_result) {
    while (my $hit = $result->next_hit) {
      push @hits, $hit;
    }
  }
  unlink("$seqFile.out");
}

my $nHits = scalar(@hits) || "no";
$nHits .= " (the maximum)" if $nHits eq $maxHits;
my $baseURL = "sites.cgi?" . (defined $hmmFile ? "hmmId=$hmmSpec" : "query=".uri_escape($query));
my $baseTitle = defined $hmmFile ? "Family Search vs. Sites" : "SitesBLAST";
print p("Found $nHits hits to proteins with known functional sites",
        "(" . a({ -href => "$baseURL&format=tsv",
                  -title => "$baseTitle results in a tab-delimited format" }, "download") . ")"), "\n"
  unless $format eq "tsv";
print join("\t", qw{queryId subjectDb subjectId identity queryAlnBegin queryAlnEnd sujectAlnBegin subjectAlnEnd
                    subjectDescription evalue bits
                    subjectSite subjectSiteRange subjectSiteSeq querySiteRange querySiteSeq pubMedIds})."\n"
  if $format eq "tsv";

my $nHit = 0;
foreach my $hit (@hits) {
  $nHit++;
  my $hitId = $hit->name;
  my $bits = $hit->bits;
  my $hsp = $hit->hsp; # will consider only the best one
  my ($queryBeg, $hitBeg) = $hsp->start('list');
  my ($queryEnd, $hitEnd) = $hsp->end('list');
  my $eval = $hit->significance(); # for an HMM, the e-value for the whole sequence
  # The aligned sequences for query and subject over the given range, with dashes
  # (For an HMM query, the sequence is partly lower-case)
  my ($alnQ, $alnS);
  if ($hmmSpec) {
    $bits = $hsp->score();
    # Non-matching positions are lower-case in HMMer output, but matching positions will be highlighted anyway
    $alnQ = uc( $hsp->query_string );
    $alnS = uc( $hsp->hit_string );
    # Convert . to - (HMMer alignments use . for gaps, but BLAST use -)
    $alnQ =~ s/[.]/-/g;
    $alnS =~ s/[.]/-/g;
  } else {
    $bits = $hit->bits;
    my $aln = $hsp->get_aln();
    $alnQ = $aln->get_seq_by_pos(1)->seq;
    $alnS = $aln->get_seq_by_pos(2)->seq;
  }

  my $seqQ = $alnQ; $seqQ =~ s/-//g;
  unless ($hmmSpec) {
    die "BLAST parsing error, query sequences do not line up:\n$seqQ\n"
      . substr($seq,$queryBeg-1,$queryEnd-$queryBeg+1)
      unless $seqQ eq substr($seq, $queryBeg-1, $queryEnd-$queryBeg+1);
  }

  # Get information about the subject
  my @hitIds = split /:/, $hitId;
  my $db = shift @hitIds;
  my $id = shift @hitIds;
  my $chain = @hitIds > 0 ? shift @hitIds : "";
  die "Invalid hit id $hitId" unless @hitIds == 0;
  my $desc = "No description";
  my $id2 = "";
  my $info =  $dbh->selectrow_hashref("SELECT * FROM HasSites WHERE db = ? AND id = ? AND chain = ?",
                                      {}, $db, $id, $chain);
  if ($info) {
    $desc = $info->{desc};
    $id2 = $info->{id2};
  }

  my $subjectLen = $info->{length};

  my $sites = $dbh->selectall_arrayref("SELECT * FROM Site WHERE db = ? AND id = ? AND chain = ? ORDER BY posFrom",
                                       { Slice => {} }, $db, $id, $chain);

  # Save alignment of positions
  my %sposToAlnPos = ();
  my %alnPosToQpos = ();
  my $atQ = $queryBeg;
  my $atS = $hitBeg;
  my $alnLen = length($alnQ);
  foreach my $alnPos (1..$alnLen) {
    my $valQ = substr($alnQ, $alnPos-1, 1);
    if ($valQ ne "-") {
      $alnPosToQpos{$alnPos} = $atQ;
      $atQ++;
    }
    my $valS = substr($alnS, $alnPos-1, 1);
    if ($valS ne "-") {
      $sposToAlnPos{$atS} = $alnPos;
      $atS++;
    }
  }

  # Save posAlnFrom and posAlnTo, with 0 for off-the-left and alnLen+1 for off-the-right
  foreach my $site (@$sites) {
    foreach my $at ("From","To") {
      my $pos = $site->{"pos".$at};
      my $posAln;
      if ($pos < $hitBeg) {
        $posAln = 0;
      } elsif ($pos > $hitEnd) {
        $posAln = $alnLen+1;
      } else {
        $posAln = $sposToAlnPos{$pos};
        die "Alignment error: no position for $pos in subject $hitId"
          unless defined $posAln;
      }
      $site->{"posAln".$at} = $posAln;
    }
  }

  # Compute union of pubmed ids
  my %pmIds = ();
  foreach my $site (@$sites) {
    foreach my $pmId (split /,/, $site->{pmIds}) {
      $pmIds{$pmId} = 1 unless $pmId eq "";
    }
  }
  my @pmIds = sort { $b <=> $a } keys %pmIds;
  my $paperLink = "";
  $paperLink = "(" . a({ -href => "http://pubmed.ncbi.nlm.nih.gov/pubmed/?term=" . join(",",@pmIds) . "&sort=pubdate",
                         -onmousedown => loggerjs($db eq "PDB" ? "PDB" : "curatedsite", $hitId) },
                       "see",
                       @pmIds > 1 ? scalar(@pmIds) . " papers" : "paper")
    . ")"
      if @pmIds > 0;

  # Header region (2 lines)
  my $hitURL = "";
  $hitURL = "https://www.rcsb.org/structure/" . $id if $db eq "PDB";
  $hitURL = "https://www.uniprot.org/uniprot/" . $id if $db eq "SwissProt";
  my $identityString = $hmmSpec ? "" : int(0.5 + 100 * $hsp->frac_identical);
  my $coverageString = int(0.5 + 100 * ($queryEnd-$queryBeg+1) / $queryLen);
  # subject length is always known right now, but not 100% sure this will always be the case
  my $fromSubjectString = $subjectLen ne "" ? "/" . $subjectLen : "";
  print p({-style => "margin-bottom: 0.5em;"},
          a({-title => $db eq "PDB" ? "PDB entry $id chain $chain" : "Swiss-Prot identifier $id2",
             -href => $hitURL},
            $id.$chain),
          $desc,
          small($paperLink),
          br(),
          small(a({-title => "$bits bits, E = $eval"},
                  $hmmSpec ? "$bits bits," : "${identityString}% identity,",
                  "${coverageString}% coverage:",
                  "${queryBeg}:${queryEnd}/${queryLen} of " . ($hmmSpec?"HMM":"query") . " aligns to",
                  "${hitBeg}:${hitEnd}${fromSubjectString} of ${id}${chain}")))
    if $format ne "tsv";

  # Fetch relevant ligand information
  my %ligands = map { $_->{ligandId} => 1 } @$sites;
  my %ligandName = ();        # only if it is information in the table
  foreach my $ligandId (keys %ligands) {
    next if $ligandId eq "" || $ligandId eq "NUC";
    my $ligInfo = $dbh->selectrow_hashref("SELECT * FROM PdbLigand WHERE ligandId = ?",
                                          {
                                          }, $ligandId);
    if (defined $ligInfo) {
      my $name = $ligInfo->{ligandName};
      $name = lc($name) unless $name =~ m/[a-z]/;
      $ligandName{$ligandId} = $name;
    }
  }

  my %sposToSite = ();
  my %alnposToSite = ();
  foreach my $site (@$sites) {
    foreach my $i ($site->{posFrom} .. $site->{posTo}) {
      push @{ $sposToSite{$i} }, $site;
      push @{ $alnposToSite{$sposToAlnPos{$i}} }, $site
        if exists $sposToAlnPos{$i};
    }
  }
  my $iSite = 0;
  foreach my $site (@$sites) {
    $site->{iSite} = $iSite++;
  }

  # For each site, compute isAligned, shortDesc, and longDesc
  # (the descriptions do not include anything about the position)
  foreach my $site (@$sites) {
    $site->{isAligned} = $site->{posAlnFrom} >= 1 && $site->{posAlnFrom} <= $alnLen
      && $site->{posAlnTo} >= 1 && $site->{posAlnTo} <= $alnLen;
    if ($db eq "PDB") {
      if ($site->{ligandId} eq "") {
        $site->{shortDesc} = "active site";
      } elsif ($site->{ligandId} eq "NUC") {
        $site->{shortDesc} = "binding DNA/RNA";
      } elsif (exists $ligandName{$site->{ligandId}}) {
        $site->{shortDesc} = "binding " . $ligandName{$site->{ligandId}};
      } else {
        $site->{shortDesc} = "binding " . $site->{ligandId};
      }
      $site->{longDesc} = $site->{shortDesc};
    } elsif ($db eq "SwissProt") {
      my $type = $site->{type};
      my $comment = $site->{comment};
      my $desc = "$type $comment";
      if ($type eq "mutagenesis") {
        if ($comment =~ m/^([A-Z]+)->([A-Z]+): (.*)$/) {
          if ($site->{isAligned}) {
            $desc = "mutation to $2: $3";
          } else {
            $desc = "$1&rightarrow;$2: $3";
          }
        } else {
          $desc = "mutation $comment";
        }
      } elsif ($type eq "functional") {
        $desc = $comment;
      } elsif ($type eq "modified") {
        if ($comment =~ m/^signal peptide/i
            || $comment =~ m/^disulfide link/) {
          $desc = $comment;
        } elsif ($comment =~ m/^natural variant/) {
          $desc = $comment;
          if ($desc =~ m/^natural variant: ([A-Z]) *-> *([A-Z]) (.*)$/) {
            my ($original, $variant, $comment2) = ($1, $2,$3);
            $comment2 =~ s/^[(]//;
            $comment2 =~ s/[)]$//;
            if ($site->{isAligned}) {
              $desc = "to $variant: $comment2";
            } else {
              $desc = "$original &rightarrow; $variant: $comment2";
            }
          }
        } else {
          $desc = "modified: $comment";
        }
      }
      $site->{longDesc} = $desc;
      my @words = split / /, $desc;
      my $maxWords = 9;
      if (@words > $maxWords) {
        $site->{shortDesc} = join(" ", splice(@words, 0, $maxWords)) . "...";
      } else {
        $site->{shortDesc} = $desc;
      }
    } else {
      die "Unknown db $db";
    }
  }

  # Alignment region -- there's 3 lines, for the query, to highlight functional sites, and for the subject
  # There's a header at the left, and then 1 column for each alignment position
  # Each column (or header) is a floating div with 3 lines (implemented using br()) and
  # with a 1em margin at top. This is necessary so that there is some space if the
  # alignment wraps (as it usually does)
  # The div for that column has id = Aln${nHit}P${i}
  # Setting the main div to be display:flex would cause it to be on one line with scrolling,
  # but would not work in opera and perhaps some other common browsers, not sure.

  # for sequences outside the bounds of the query)
  my $queryText = $hmmSpec? "HMM" : "query";
  my $queryColumns = FormatAlnString($alnQ, $queryBeg, {},
                                     $queryText,
                                     $hitBeg, $alnS, $id.$chain);
  my $siteColumns = FormatAlnSites($alnQ, $alnS, $queryBeg, $hitBeg, $id.$chain, \%sposToSite);
  my $hitColumns = FormatAlnString($alnS, $hitBeg, \%sposToSite, $id.$chain, $queryBeg, $alnQ, $queryText);
  my @alnColumns = ();
  foreach my $i (0..($alnLen-1)) {
    my @onMouseOver = ();
    my @onMouseOut = ();
    # use i+1 in the id because posAlnFrom and such are 1-based, with 0 for off-to-the-left
    if (exists $alnposToSite{$i+1}) {
      foreach my $site (@{ $alnposToSite{$i+1} }) {
        my $iSite = $site->{iSite};
        push @onMouseOver, qq{document.getElementById("Site${nHit}S${iSite}").style.backgroundColor="lightgrey";};
        push @onMouseOut, qq{document.getElementById("Site${nHit}S${iSite}").style.backgroundColor="white";};
      }
    }
    my %colArgs = (-class => "alnCol", -id => "Aln${nHit}P" . ($i+1) );
    $colArgs{"-onmouseover"} = join("", @onMouseOver) if @onMouseOver > 0;
    $colArgs{"-onmouseout"} = join("", @onMouseOut) if @onMouseOut > 0;
    push @alnColumns, div(\%colArgs,
                          join(br(),
                               $queryColumns->[$i],
                               $siteColumns->[$i],
                               $hitColumns->[$i]));
  }
  my $alnTitle = ($hmmSpec? "HMM consensus for positions" : "amino acids")
    . " $queryBeg to $queryEnd/$queryLen";
  my @alnLabels = (a({-title => $alnTitle}, "$queryText"),
                   "sites",
                   a({-title => "amino acids $hitBeg to $hitEnd$fromSubjectString"}, $id.$chain ));
  print "\n",
    div( div({-class => "alnLabel"},
             join(br(), @alnLabels)),
         "\n",
         @alnColumns),
    "\n", div({-style => "clear:both;"}, ""), "\n"
      if $format ne "tsv";

  my @alignedSites = grep $_->{isAligned}, @$sites;
  my @unalignedSites = grep ! $_->{isAligned}, @$sites;
  foreach my $isAligned (1,0) {
    my $sitesThisAlign = $isAligned ? \@alignedSites : \@unalignedSites;
    next unless @$sitesThisAlign > 0;
    print p("Sites not aligning to the query:") if !$isAligned && $format ne "tsv";
    my @bullets = ();

    if ($db eq "PDB") {
      my %byLigand = ();
      foreach my $site (@$sitesThisAlign) {
        push @{ $byLigand{$site->{ligandId} } }, $site;
      }
      foreach my $ligandId (sort keys %byLigand) {
        my $ligSites = $byLigand{$ligandId};
        my $ligShow = "";
        if ($ligandId eq "") {
          $ligShow = "active site:";
        } elsif ($ligandId eq "NUC") {
          $ligShow = "binding DNA/RNA:";
        } elsif (exists $ligandName{$ligandId}) {
          $ligShow = a({ -href => "http://www.rcsb.org/ligand/".$ligandId,
                         -title => "PDB ligand id " . $ligandId },
                       "binding " . $ligandName{$ligandId} . ":");
        } else {
          $ligShow = "binding ${ligShow}:";
        }
        my @posShow = ();
        foreach my $site (@$ligSites) {
          my $posShow = "";
          my $sChar = substr($alnS, $site->{posAlnFrom}-1, 1) if $isAligned;
          $posShow .= $sChar if $isAligned;
          $posShow .= $site->{posFrom};
          $posShow = a({-title => "$site->{pdbFrom} in PDB numbering for $id$chain"},
                       $posShow)
            if $site->{pdbFrom} ne "";
          my ($qChar, $qPos);
          if ($isAligned) {
            $qChar = substr($alnQ, $site->{posAlnFrom}-1, 1);
            if ($qChar eq "-") {
              $posShow .= " (vs. gap)";
            } else {
              $qPos = $alnPosToQpos{ $site->{posAlnFrom} };
              my $qShow = a({-title => "$qChar$qPos in " . ($hmmSpec?"HMM":"query") }, $qChar.$qPos);
              my $matchChar = $qChar eq $sChar ? "=" : "<span style='color:darkred'>&ne;</span>";
              $posShow .= " ($matchChar $qShow)";
            }
            $posShow = AddMouseOver($posShow, $nHit, $site->{posAlnFrom});
          }
          $posShow = span({-id => "Site${nHit}S" . $site->{iSite} }, $posShow);
          push @posShow, $posShow;
          print join("\t", $header, $db,
                     $db eq "PDB" ? "$id:$chain" : $id2,
                     $hmmSpec? "" : sprintf("%.1f", 100 * $hsp->frac_identical),
                     $queryBeg, $queryEnd, $hitBeg, $hitEnd,
                     $desc, $eval, $bits,
                     $site->{shortDesc},
                     $site->{posFrom}, $sChar || "",
                     $qPos || "", $qChar || "",
                    $site->{pmIds})."\n"
                       if $format eq "tsv";
        }
        push @bullets, li($ligShow, join(", ", @posShow));
      }
    } elsif ($db eq "SwissProt") {
      # Sometimes have multiple indications for a single position, so, try to collapse
      my %byPos = ();           # from => to => sites
      foreach my $site (@$sitesThisAlign) {
        push @{ $byPos{ $site->{posFrom} }{ $site->{posTo} } }, $site;
      }
      my @typeOrder = qw{functional binding modified mutagenesis};
      my %typeOrder = map { $typeOrder[$_] => $_ } 0..(scalar(@typeOrder)-1);
      foreach my $posFrom (sort {$a<=>$b} keys %byPos) {
        my $hashTo = $byPos{$posFrom};
        foreach my $posTo (sort {$a<=>$b} keys %$hashTo) {
          my $sitesHere = $hashTo->{$posTo};
          foreach my $site (@$sitesHere) {
            my $type = $site->{type};
            $site->{orderScore} = exists $typeOrder{$type} ? $typeOrder{$type} : scalar(@typeOrder);
          }
          my @sitesHere = sort { $a->{orderScore} <=> $b->{orderScore} } @$sitesHere;
          die unless @sitesHere > 0;
          my $site1 = $sitesHere[0];
          my $showPos = "";
          my ($sSeq,$qSeq);
          if ($isAligned) {
            $sSeq = substr($alnS, $site1->{posAlnFrom}-1, $site1->{posAlnTo}-$site1->{posAlnFrom}+1);
            $qSeq = substr($alnQ, $site1->{posAlnFrom}-1, $site1->{posAlnTo}-$site1->{posAlnFrom}+1);
          }
          # If it's a multi-position site, and either the 1st or last position aligns to a gap,
          # we can't just use the full range to get the range in the query. Instead, we need to
          # see if there is a corresponding range.
          my @qPos;
          foreach my $alnPos ($site1->{posAlnFrom} .. $site1->{posAlnTo}) {
            push @qPos, $alnPosToQpos{$alnPos} if exists $alnPosToQpos{$alnPos};
          }
          my ($qFirst, $qLast);
          if (@qPos > 0) {
            $qFirst = $qPos[0];
            $qLast = $qPos[-1];
          }
          if ($posTo - $posFrom <= 5) {
            $showPos .= $sSeq if $isAligned;
            $showPos .= $posFrom eq $posTo ? $posFrom : " $posFrom:$posTo";
            if ($isAligned) {
              if (!defined $qFirst) {
                $showPos .= " (vs. gap)";
              } else {
                if ($qSeq eq $sSeq) {
                  $showPos .= " (= ";
                } else {
                  $showPos .= " (<span style='color:darkred'>&ne;</span> ";
                }
                my $qShow = $qSeq;
                $qShow .= " " if $posTo ne $posFrom;
                if (@qPos > 0) {
                  $qShow .= $qFirst;
                  $qShow .= ":" . $qLast if $qFirst != $qLast;
                  $qShow = a({-title => "$qShow in " . ($hmmSpec?"HMM":"query")}, $qShow);
                }
                $showPos .= $qShow . ")";
              }
            }
          } else {
            # large region, report %conserved if aligned
            $showPos .= "$posFrom:$posTo";
            if ($isAligned && $qFirst) {
              # compute %conservation and then
              my $regionQ = substr($alnQ, $site1->{posAlnFrom}-1, $site1->{posAlnTo}-$site1->{posAlnFrom}+1);
              my $regionS = substr($alnS, $site1->{posAlnFrom}-1, $site1->{posAlnTo}-$site1->{posAlnFrom}+1);
              die unless length($regionQ) == length($regionS);
              my $n = length($regionQ);
              die if $n < 1;
              my $nMatch = 0;
              for (my $i = 0; $i < $n; $i++) {
                $nMatch++ if substr($regionQ,$i,1) eq substr($regionS,$i,1);
              }
              my $percMatch = int(0.5 + 100 * $nMatch/$n);
              $showPos .= " (vs. $qFirst:$qLast, ${percMatch}% identical)";
            } elsif ($isAligned ){
              $showPos .= " (vs. gap)";
            }
          }
          $showPos = AddMouseOver($showPos, $nHit, $site1->{posAlnFrom})
            if $isAligned;
          my @siteDesc = map span({-id => "Site${nHit}S" . $_->{iSite} }, $_->{longDesc}), @sitesHere;
          push @bullets, li($showPos, join("; ", @siteDesc));
          my %pmHere = (); # all pmids across @sitesHere
          foreach my $site (@sitesHere) {
            foreach my $pmId (split /,/, $site->{pmIds}) {
              $pmHere{$pmId} = 1 if $pmId ne "";
            }
          }
          print join("\t", $header, $db,
                     $db eq "PDB" ? $id.$chain : $id2,
                     $hmmSpec? "" : sprintf("%.1f", 100 * $hsp->frac_identical),
                     $queryBeg, $queryEnd, $hitBeg, $hitEnd,
                     $desc, $eval, $bits,
                     join("; ", map $_->{shortDesc}, @$sitesHere),
                     $posFrom eq $posTo ? $posFrom : "$posFrom..$posTo",
                     $sSeq || "",
                     defined $qFirst ? ($qFirst == $qLast ? $qFirst : "$qFirst..$qLast") : "",
                     $qSeq || "",
                     join(",", sort keys %pmHere))."\n"
                       if $format eq "tsv";
        }                       # end loop over PosTo
      }                         # end loop over PosFrom
    } else {                    # not PDB or SwissProt
      die "Unknown db $db";
    }
    print start_ul(), @bullets, end_ul(), "\n" if $format ne "tsv";
  }                             # end loop over isAligned
}                               # end loop over hits

exit(0) if $format eq "tsv";

if (! $hmmSpec) {
  my @pieces = $seq =~ /.{1,60}/g;
  print
    h3("Query Sequence"),
    p({-style => "font-family: monospace;"}, small(join(br(), ">" . HTML::Entities::encode($header), @pieces))),
    p("Or try a", a({-href => "sites.cgi"}, "new SitesBLAST search"));
}
print  h3("SitesBLAST's Database"),
  p($docstring);
finish_page();


sub FormatAlnSites($$$$$$) {
  my ($queryAln, $hitAln, $queryBeg, $hitBeg, $hitName, $hposSite) = @_;
  my @out = ();
  my $queryAt = $queryBeg;
  my $hitAt = $hitBeg;
  die unless length($hitAln) == length($queryAln);
  for (my $i = 0; $i < length($hitAln); $i++) {
    my $queryChar = substr($queryAln, $i, 1);
    my $hitChar = substr($hitAln, $i, 1);
    my $string = "&nbsp;";
    my $class = "alnS";
    my $title;
    if ($hitChar ne "-" && exists $hposSite->{$hitAt}) { # site in subject
      my $sitesHere = $hposSite->{$hitAt};
      # Old logic ste $hasSite1 only if there was a site that covered
      # this character only. Decided this was wierd.
      if ($queryChar eq $hitChar) {
        $string = "&vert;";
        $class = "alnS1";
      } else {
        $string = "x";
        $class = "alnS0";
      }
      my $queryLong = exists $charToLong{$queryChar} ? $charToLong{$queryChar} : $queryChar;
      my $hitLong = exists $charToLong{$hitChar} ? $charToLong{$hitChar} : $hitChar;
      $title = "${hitLong}${hitAt} in $hitName: "
        .  join("; ", map $_->{shortDesc}, @$sitesHere);
      $title .= " (${queryLong}${queryAt} in query)" if $queryChar ne "-";
    }
    if (defined $title) {
      push @out, a({-title => $title, -class => $class}, $string);
    } else {
      push @out, a({-class => $class}, $string);
    }
    $queryAt++ unless $queryChar eq "-";
    $hitAt++ unless $hitChar eq "-";
  }
  return \@out;
}

sub FormatAlnString($$$$$$$) {
  my ($alnSeq, $beg, $posToSite, $seqName, $begOther, $alnOther, $nameOther) = @_;

  if (keys(%charToSet) == 0) {
    while (my ($chars, $setno) = each %charSets) {
      foreach my $char (split //, $chars) {
        $charToSet{$char} = $setno;
      }
    }
  }

  my @out = ();
  my $at = $beg;
  my $atOther = $begOther;
  for (my $i = 0; $i < length($alnSeq); $i++) {
    my $char = substr($alnSeq, $i, 1);
    my $charOther = substr($alnOther, $i, 1) if defined $alnOther;

    my $class = exists $charToSet{$char} ? "aa" . $charToSet{$char} : "aa";
    my $longAA = exists $charToLong{$char} ? $charToLong{$char} : $char;
    my $title;

    if ($char ne "-") {
      $title = "${longAA}${at} in $seqName";
      if (defined $charOther && $charOther ne "-") {
        my $longOther = exists $charToLong{$charOther} ? $charToLong{$charOther} : $char;
        $title .= " (${longOther}${atOther} in $nameOther)";
      }
    }
    my %attrs = ( -class => $class );
    $attrs{title} = $title if defined $title;
    push @out, span(\%attrs, $char);

    $at++ if $char ne "-";
    $atOther++ if defined $charOther && $charOther ne "-";
  }
  return \@out;
}

sub AddMouseOver($$$) {
  my ($html, $nHit, $iAln) = @_;
  return $html unless $iAln >= 0;
  my $id = "Aln${nHit}P${iAln}";
  return span({ -onmouseover => "alnHighlight('$id', 1)",
                -onmouseout => "alnHighlight('$id', 0)" },
              $html);
}


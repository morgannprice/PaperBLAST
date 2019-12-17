#!/usr/bin/perl -w
use strict;

# sites.cgi -- compare a protein sequence to known functional sites
#
# Optional CGI parameters:
# query -- this should be the protein sequence in FASTA or UniProt or plain format
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
use pbweb qw{start_page finish_page GetMotd loggerjs UniProtToFasta RefSeqToFasta VIMSSToFasta FBrowseToFasta};
use Bio::SearchIO;

sub fail($);
sub FormatAlnString($$$$);

my %charToLong = ("A" => "Ala", "C" => "Cys", "D" => "Asp", "E" => "Glu",
                  "F" => "Phe", "G" => "Gly", "H" => "His", "I" => "Ile",
                  "K" => "Lys", "L" => "Leu", "M" => "Met", "N" => "Asn",
                  "P" => "Pro", "Q" => "Gln", "R" => "Arg", "S" => "Ser",
                  "T" => "Thr", "V" => "Val", "W" => "Trp", "Y" => "Tyr",
                  "O" => "Pyrrolysine", "U" => "Selenocysteine");
my %charColSets = ("lightblue" => "AILMFWV",
                   "red" => "KR",
                   "magenta" => "ED",
                   "lightgreen" => "NQST",
                   "pink" => "C",
                   "orange" => "G",
                   "yellow" => "P",
                   "cyan" => "HY");
my %charToColor = ();

my $tmpDir = "../tmp";
my $blastall = "../bin/blast/blastall";
my $nCPU = 4;
my $base = "../data";
my $blastdb = "$base/hassites.faa";
my $sqldb = "$base/sites.db";
my $fbdata = "../fbrowse_data"; # path relative to the cgi directory

my $procId = $$;
my $timestamp = int (gettimeofday() * 1000);
my $filename = $procId . $timestamp;
my $seqFile = "$tmpDir/$filename.fasta";

sub fail($) {
    my ($notice) = @_;
    print
        p($notice),
        p(a({-href => "litSearch.cgi"}, "New search")),
        end_html;
    exit(0);
}

my $cgi=CGI->new;
my $query = $cgi->param('query') || "";
my $maxHits = 20;
my $maxE = 0.001;

my $title = "SitesBLAST";
start_page('title' => $title);
print <<END
<SCRIPT src="../static/pb.js"></SCRIPT>
END
;

# remove leading and trailing whitespace
$query =~ s/^\s+//;
$query =~ s/\s+$//;
unless ($query) {
  print
    GetMotd(),
    start_form( -name => 'input', -method => 'GET', -action => 'sites.cgi'),
    p(br(),
      b("Enter a protein sequence in FASTA format, or an identifier from UniProt, RefSeq, or MicrobesOnline"),
      br(),
      textarea( -name  => 'query', -value => '', -cols  => 70, -rows  => 10 )),
    p(submit('Search'), reset()),
    end_form;
  finish_page();
} else {
  my ($header, $seq);

  # Single-line query sequence
  if ($query !~ m/\n/ && $query =~ m/^[A-Z*]+$/) {
    $header = substr($seq, 0, 10) . "..." . "(" . length($query) . " amino acids)";
    $seq = $query;
  } elsif ($query !~ m/\n/) {
    # single-line, identifier only
    my $short = $query;
    $query = undef;

    # Is it a VIMSS id?
    $query = &VIMSSToFasta($short)
      if $short =~ m/^VIMSS\d+$/i;
    $query = &FBrowseToFasta($fbdata, $short)
      if !defined $query && $short =~ m/^[0-9a-zA-Z_]+$/;
    $query = &UniProtToFasta($short)
      if !defined $query && $short =~ m/^[A-Z0-9_]+$/;
    $query = &RefSeqToFasta($short)
      if !defined $query && $short =~ m/^[A-Z0-9._]+$/;
    my $shortSafe = HTML::Entities::encode($short);
    &fail("Sorry -- we were not able to find a protein sequence for the identifier <b>$shortSafe</b>. We checked it against MicrobesOnline, UniProt, and the NCBI protein database (RefSeq and Genbank). Please use the sequence as a query instead.")
      unless defined $query;
  }

  if (!defined $header) {
    fail("Not in fasta format") unless $query =~ m/^>/;
    my @lines = split /\n/, $query;
    $header = shift @lines;
    $header =~ s/^>//;
    $seq = join("", @lines);
    $seq =~ s/\s//g;
    fail("Invalid sequence") unless $seq =~ m/^[A-Z*]+$/;
    $seq =~ s/[*]//g;
  }
  fail("Sequence is too short") unless length($seq) >= 10;

  my $dbh = DBI->connect("dbi:SQLite:dbname=$sqldb","","",{ RaiseError => 1 }) || die $DBI::errstr;
  autoflush STDOUT 1; # show preliminary results
  print p("Comparing $header to proteins with known functional sites using BLASTp with E &le; $maxE"), "\n";
  open(my $fhFaa, ">", $seqFile) || die "Cannot write to $seqFile\n";
  print $fhFaa ">$header\n$seq\n";
  close($fhFaa) || die "Error writing to $seqFile\n";
  die "No such executable: $blastall\n" unless -x $blastall;
  # m S means mask complex sequences for lookup but not for alignment
  system("$blastall -F 'm S' -p blastp -i $seqFile -d $blastdb -e $maxE -a $nCPU -o $seqFile.out -b $maxHits -v $maxHits >& /dev/null") == 0
    || die "$blastall failed: $!\n";
  unlink($seqFile);
  my $searchio = Bio::SearchIO->new(-format => 'blast', -file => "$seqFile.out")
    || die "Failed to read $seqFile.out\n";
  my @hits = ();
  while (my $result = $searchio->next_result) {
    while (my $hit = $result->next_hit) {
      push @hits, $hit;
    }
  }
  my $nHits = scalar(@hits) || "no";
  $nHits .= " (the maximum)" if $nHits eq $maxHits;
  print p("Found $nHits hits to proteins with known functional sites"),"\n";
  unlink("$seqFile.out");
  foreach my $hit (@hits) {
    my $hitId = $hit->name;
    my $bits = $hit->bits;
    my $eval = $hit->significance;
    my $hsp = $hit->hsp; # will consider only the best one
    my ($queryBeg, $hitBeg) = $hsp->start('list');
    my ($queryEnd, $hitEnd) = $hsp->end('list');
    my $queryLen = length($seq);
    my $aln = $hsp->get_aln();
    # The aligned sequences for query and subject over the given range, with dashes
    my $alnQ = $aln->get_seq_by_pos(1)->seq;
    my $alnS = $aln->get_seq_by_pos(2)->seq;
    my $seqQ = $alnQ; $seqQ =~ s/-//g;
    die "BLAST parsing error, query sequences do not line up:\n$seqQ\n"
      . substr($seq,$queryBeg-1,$queryEnd-$queryBeg+1)
      unless $seqQ eq substr($seq, $queryBeg-1, $queryEnd-$queryBeg+1);

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

    my $sites = $dbh->selectall_arrayref("SELECT * FROM Sites WHERE db = ? AND id = ? AND chain = ? ORDER BY posFrom",
                                  { Slice => {} }, $db, $id, $chain);

    # Save alignment of positions
    my %sposToAlnPos = ();
    my %alnPosToQpos = ();
    my $atQ = $queryBeg;
    my $atS = $hitBeg;
    my $alnLen = $aln->length();
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

    # Header region (2 lines)
    my $hitURL = "";
    $hitURL = "https://www.rcsb.org/structure/" . $id if $db eq "PDB";
    $hitURL = "https://www.uniprot.org/uniprot/" . $id if $db eq "SwissProt";
    my $identityString = int(100 * $hsp->frac_identical);
    my $coverageString = int(100 * ($queryEnd-$queryBeg+1) / length($seq));
    print p({-style => "margin-bottom: 0.5em;"},
            a({-title => $db eq "PDB" ? "entry $id chain $chain" : $id2,
               -href => $hitURL},
              $id.$chain),
            $desc,
            br(),
            small(a({-title => "$bits bits, E = $eval"},
                    "${identityString}% identity")
                  . ", "
                  . a({"-title" => "${queryBeg}:${queryEnd}/${queryLen} of query (${coverageString}%) aligns to ${hitBeg}:${hitEnd} of ${id}${chain}" },
                      "${coverageString}% coverage")));

    my %sposToSite = ();
    foreach my $site (@$sites) {
      foreach my $i ($site->{posFrom} .. $site->{posTo}) {
        push @{ $sposToSite{$i} }, $site;
      }
    }

    # Alignment region (no-wrap, fixed font, 1 extra space at beginning and end
    # for sequences outside the bounds of the query)
    my @lines = (i(a({-title => "$queryBeg to $queryEnd/$queryLen of query"}, "Query:"))
                 . " "
                 . FormatAlnString($alnQ, $queryBeg, {}, "query"),
                 i(a({-title => "$hitBeg to $hitEnd of $id$chain"}, "Sbjct:"))
                 . " "
                 . FormatAlnString($alnS, $hitBeg, \%sposToSite, $id.$chain) . " ");
    print div({-style => "font-family: monospace; white-space:nowrap;" },
              p({-style => "margin-top: 0em; margin-bottom: 2em;"},
                join(br(), @lines)));
    print "\n";

    my @siteShow = ();
    foreach my $site (@$sites) {
      my $showLigand = $site->{ligandId};
      if ($site->{ligandId} eq "NUC") {
        $showLigand = "DNA or RNA";
      } else {
        $showLigand = a({-href => "http://www.rcsb.org/ligand/".$site->{ligandId} },
                        $site->{ligandId}) if $site->{ligandId};
      }
      my @parts = ($site->{posFrom} . ":" . $site->{posTo},
                   $site->{type}, $showLigand, $site->{comment});
      if ($site->{posAlnFrom} >= 1
          && $site->{posAlnFrom} <= $alnLen
          && $site->{posAlnTo} >= 1
          && $site->{posAlnTo} <= $alnLen) {
        push @parts, "(" . substr($alnS, $site->{posAlnFrom}-1, $site->{posAlnTo}-$site->{posAlnFrom}+1)
          . " vs. " . substr($alnQ, $site->{posAlnFrom}-1, $site->{posAlnTo}-$site->{posAlnFrom}+1)
            . ")"
              if $site->{posAlnTo} - $site->{posAlnFrom} + 1 <= 10;
      } else {
        push @parts, "(outside the alignment)";
      }
      my @pmIds = split /,/, $site->{pmIds};
      my $note = @pmIds > 1 ? scalar(@pmIds) . " papers" : "paper";
      push @parts,
        small("see",
              a({ -href => "http://www.ncbi.nlm.nih.gov/pubmed/" . join(",",@pmIds),
            -onmousedown => loggerjs($db eq "PDB" ? "PDB" : "curatedsite", $hitId) },
          $note))
          if @pmIds > 0;
      push @siteShow, li(@parts);
    }
    print start_ul(), @siteShow, end_ul(), "\n";
  }
  my @pieces = $seq =~ /.{1,60}/g;
  print
    h3("Query Sequence"),
    p({-style => "font-family: monospace;"}, small(join(br(), ">" . HTML::Entities::encode($header), @pieces))),
    h3(a({-href => "sites.cgi"}, "New Search"));

  finish_page();
}

sub FormatAlnString($$$$) {
  my ($alnSeq, $beg, $posToSite, $seqName) = @_;

  if (keys(%charToColor) == 0) {
    while (my ($color, $chars) = each %charColSets) {
      foreach my $char (split //, $chars) {
        $charToColor{$char} = $color;
      }
    }
  }

  my @out = ();
  my $at = $beg;
  for (my $i = 0; $i < length($alnSeq); $i++) {
    my $char = substr($alnSeq, $i, 1);
    if ($char eq "-") {
      push @out, "-";
    } else {
      my @styleparts = ();
      push @styleparts, "background-color: " . $charToColor{$char} . ";"
        if exists $charToColor{$char};
      if (exists $posToSite->{$at}) {
        my $n = @{ $posToSite-> {$at} };
        push @styleparts, "font-weight: bold; text-decoration-line: underline; text-decoration-style: "
          . ($n > 1 ? "double" : "solid") . ";";
      }
      my $longAA = exists $charToLong{$char} ? $charToLong{$char} : $char;
      push @out, a({-title => "${longAA}${at} in $seqName",
                    -style => join(" ",@styleparts)}, $char);
      $at++;
    }
  }
  return join("", @out);
}

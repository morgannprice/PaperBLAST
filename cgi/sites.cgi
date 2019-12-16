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
use pbweb qw{start_page finish_page GetMotd loggerjs};
use Bio::SearchIO;

sub fail($);
my $tmpDir = "../tmp";
my $blastall = "../bin/blast/blastall";
my $nCPU = 4;
my $base = "../data";
my $blastdb = "$base/hassites.faa";
my $sqldb = "$base/sites.db";

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
      b("Enter a protein sequence in FASTA format"),
      br(),
      textarea( -name  => 'query', -value => '', -cols  => 70, -rows  => 10 )),
    p(submit('Search'), reset()),
    end_form;
  finish_page();
} else {
  my ($header, $seq);
  if ($query !~ m/\n/ && $query =~ m/^[A-Z*]+$/) {
    $header = substr($seq, 0, 10) . "..." . "(" . length($query) . " amino acids)";
    $seq = $query;
  } elsif ($query =~ m/^>/) {
    my @lines = split /\n/, $query;
    $header = shift @lines;
    $header =~ s/^>//;
    $seq = join("", @lines);
    $seq =~ s/\s//g;
    fail("Invalid sequence") unless $seq =~ m/^[A-Z*]+$/;
    $seq =~ s/[*]//g;
    fail("Sequence is too short") unless length($seq) >= 10;
  } else {
    fail("Not in fasta format");
  }
  my $dbh = DBI->connect("dbi:SQLite:dbname=$sqldb","","",{ RaiseError => 1 }) || die $DBI::errstr;
  autoflush STDOUT 1; # show preliminary results
  print p("Comparing $header to proteins with known functional sites using BLASTp with E &le; $maxE"), "\n";
  open(my $fhFaa, ">", $seqFile) || die "Cannot write to $seqFile\n";
  print $fhFaa ">$header\n$seq\n";
  close($fhFaa) || die "Error writing to $seqFile\n";
  die "No such executable: $blastall\n" unless -x $blastall;
  system("$blastall -p blastp -i $seqFile -d $blastdb -e $maxE -a $nCPU -o $seqFile.out -b $maxHits -v $maxHits >& /dev/null") == 0
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
    die "BLAST parsing error, query sequences do not line up"
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

    my @rows = ();
    foreach my $site (@$sites) {
      my $okRow;
      foreach my $iRow (0 .. (scalar(@rows) - 1)) {
        next unless $iRow >= 0 && @rows > $iRow;
        my $ok = 1;
        my $row = $rows[$iRow];
        foreach my $site2 (@$row) {
          if ($site2->{posAlnFrom} <= $site->{posAlnTo}
              && $site2->{posAlnTo} >= $site->{posAlnFrom}) {
            $ok = 0;
            last;
          }
        }
        if ($ok) {
          $okRow = $iRow;
          last;
        }
      }
      if (!defined $okRow) {
        push @rows, [];
        $okRow = scalar(@rows)-1;
      }
      push @{ $rows[$okRow] }, $site;
      $site->{iRow} = $okRow;
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

    # Alignment region (no-wrap, fixed font, 1 extra space at beginning and end for sequences outside the bounds of the query)
    my $formatQ = $alnQ;
    my $formatS = $alnS;
    my @lines = (i("Query:   ") . " " . $formatQ . " ",
                 i("Subject: ") . " " . $formatS . " ");
    print div({-style => "font-family: monospace; white-space:nowrap;" },
              pre({-style => "margin-top: 0em; margin-bottom: 2em;"},
                join(br(), @lines)));
    print "\n";

    my @siteShow = ();
    foreach my $site (@$sites) {
      my $showLigand = $site->{ligandId};
      $showLigand = a({-href => "http://www.rcsb.org/ligand/".$site->{ligandId} },
                      $site->{ligandId}) if $site->{ligandId};
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
  finish_page();
}


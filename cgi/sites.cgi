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
my %charSets = ("AILMFWV" => 1,
                "KR" => 2,
                "ED" => 3,
                "NQST" => 4,
                "C" => 5,
                "G" => 6,
                "P" => 7,
                "HY" => 8);
my %charToSet = ();

my $mygrey = "#EEEEEE";

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
  if ($query !~ m/\n/ && $query =~ m/^[A-Z*]+$/ && length($query) >= 5) {
    $header = substr($seq, 0, 10) . "..." . "(" . length($query) . " amino acids)";
    $seq = $query;
  } elsif ($query !~ m/\n/ && $query =~ m/^[a-zA-Z][a-zA-Z90-9_.]+$/) {
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
    fail("Sorry, input is not in fasta format") unless $query =~ m/^>/;
    my @lines = split /\n/, $query;
    $header = shift @lines;
    $header =~ s/^>//;
    $seq = join("", @lines);
    $seq =~ s/\s//g;
    fail("Invalid sequence, only uppercase letters and * (for stop codons) are allowed") unless $seq =~ m/^[A-Z*]+$/;
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
  system("$blastall -F 'm S' -p blastp -i $seqFile -d $blastdb -e $maxE -a $nCPU -o $seqFile.out -b $maxHits -v $maxHits -a $nCPU >& /dev/null") == 0
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
  my $nHit = 0;
  foreach my $hit (@hits) {
    $nHit++;
    my $hitId = $hit->name;
    my $bits = $hit->bits;
    my $eval = $hit->significance;
    my $hsp = $hit->hsp; # will consider only the best one
    my ($queryBeg, $hitBeg) = $hsp->start('list');
    my ($queryEnd, $hitEnd) = $hsp->end('list');
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

    my $queryLen = length($seq);
    my $subjectLen = $info->{length};

    my $sites = $dbh->selectall_arrayref("SELECT * FROM Site WHERE db = ? AND id = ? AND chain = ? ORDER BY posFrom",
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

    # Compute union of pubmed ids
    my %pmIds = ();
    foreach my $site (@$sites) {
      foreach my $pmId (split /,/, $site->{pmIds}) {
        $pmIds{$pmId} = 1 unless $pmId eq "";
      }
    }
    my @pmIds = sort { $b <=> $a } keys %pmIds;
    my $paperLink = "";
    $paperLink = "(" . a({ -href => "http://www.ncbi.nlm.nih.gov/pubmed/" . join(",",@pmIds),
                          -onmousedown => loggerjs($db eq "PDB" ? "PDB" : "curatedsite", $hitId) },
                         "see",
                         @pmIds > 1 ? scalar(@pmIds) . " papers" : "paper")
               . ")"
      if @pmIds > 0;

    # Header region (2 lines)
    my $hitURL = "";
    $hitURL = "https://www.rcsb.org/structure/" . $id if $db eq "PDB";
    $hitURL = "https://www.uniprot.org/uniprot/" . $id if $db eq "SwissProt";
    my $identityString = int(0.5 + 100 * $hsp->frac_identical);
    my $coverageString = int(0.5 + 100 * ($queryEnd-$queryBeg+1) / length($seq));
    # subject length is always known right now, but not 100% sure this will always be the case
    my $fromSubjectString = $subjectLen ne "" ? "/" . $subjectLen : "";
    print p({-style => "margin-bottom: 0.5em;"},
            a({-title => $db eq "PDB" ? "entry $id chain $chain" : $id2,
               -href => $hitURL},
              $id.$chain),
            $desc,
            small($paperLink),
            br(),
            small(a({-title => "$bits bits, E = $eval"},
                    "${identityString}% identity,",
                      "${coverageString}% coverage:",
                    "${queryBeg}:${queryEnd}/${queryLen} of query aligns to ${hitBeg}:${hitEnd}${fromSubjectString} of ${id}${chain}")));

    # Fetch relevant ligand information
    my %ligands = map { $_->{ligandId} => 1 } @$sites;
    my %ligandName = (); # only if it is information in the table
    foreach my $ligandId (keys %ligands) {
      next if $ligandId eq "" || $ligandId eq "NUC";
      my $ligInfo = $dbh->selectrow_hashref("SELECT * FROM PdbLigand WHERE ligandId = ?",
                                              {}, $ligandId);
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
    my $queryColumns = FormatAlnString($alnQ, $queryBeg, {}, "query", $hitBeg, $alnS, $id.$chain);
    my $siteColumns = FormatAlnSites($alnQ, $alnS, $queryBeg, $hitBeg, $id.$chain, \%sposToSite);
    my $hitColumns = FormatAlnString($alnS, $hitBeg, \%sposToSite, $id.$chain, $queryBeg, $alnQ, "query");
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
    my @alnLabels = (a({-title => "$queryBeg to $queryEnd/$queryLen of query"}, "Query:"),
                     "Sites:",
                     a({-title => "$hitBeg to $hitEnd$fromSubjectString of $id$chain"}, "Subject:"));
    print "\n",
      div({-style => qq{display:inline; font-family:"Courier New",monospace; } },
      div({-style => "vertical-align: top; float:left; display:inline-block; font-style: italic; margin-top: 1em;"},
          join(br(), @alnLabels)),
          "\n",
          @alnColumns),
      "\n",
      div({-style => "clear:both;"}, ""),
      "\n";

    my @alignedSites = grep $_->{isAligned}, @$sites;
    my @unalignedSites = grep ! $_->{isAligned}, @$sites;
    foreach my $isAligned (1,0) {
      my $sitesThisAlign = $isAligned ? \@alignedSites : \@unalignedSites;
      next unless @$sitesThisAlign > 0;
      print p("Sites not aligning to the query:") if !$isAligned;
      my @bullets = ();

      if ($db eq "PDB") {
        my %byLigand = ();
        foreach my $site (@$sitesThisAlign) {
          push @{ $byLigand{$site->{ligandId} } }, $site;
        }
        foreach my $ligandId (sort keys %byLigand) {
          my $ligSites = $byLigand{$ligandId};
          my $ligShow;
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
            my $qChar = substr($alnS, $site->{posAlnFrom}-1, 1) if $isAligned;
            $posShow .= $qChar if $isAligned;
            $posShow .= $site->{posFrom};
            $posShow = a({-title => "$site->{pdbFrom} in PDB numbering for $id$chain"},
                         $posShow)
              if $site->{pdbFrom} ne "";
            if ($isAligned) {
              my $qChar = substr($alnQ, $site->{posAlnFrom}-1, 1);
              $posShow .= " (vs. $qChar"
                . ($qChar eq "-" ? $alnPosToQpos{ $site->{posAlnFrom} } : "")
                . ")";
              $posShow = AddMouseOver($posShow, $nHit, $site->{posAlnFrom});
            }
            $posShow = span({-id => "Site${nHit}S" . $site->{iSite} }, $posShow);
            push @posShow, $posShow;
          }
          push @bullets, li($ligShow, join(", ", @posShow));
        }
      } elsif ($db eq "SwissProt") {
        # Sometimes have multiple indications for a single position, so, try to collapse
        my %byPos = (); # from => to => sites
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
            if ($posTo - $posFrom <= 5) {
              $showPos .= substr($alnS, $site1->{posAlnFrom}-1, $site1->{posAlnTo}-$site1->{posAlnFrom}+1)
                if $isAligned;
              $showPos .= $posFrom eq $posTo ? $posFrom : " $posFrom:$posTo";
              if ($isAligned) {
                $showPos .= " (vs. " . substr($alnQ, $site1->{posAlnFrom}-1, $site1->{posAlnTo}-$site1->{posAlnFrom}+1);
                $showPos .= " " if $posTo ne $posFrom;
                if (exists $alnPosToQpos{$site1->{posAlnFrom}} && exists $alnPosToQpos{$site1->{posAlnFrom}}) {
                  $showPos .= $alnPosToQpos{$site1->{posAlnFrom}};
                  $showPos .= ":" . $alnPosToQpos{$site1->{posAlnFrom}}
                    if $posTo ne $posFrom;
                }
                $showPos .= ")";
              }
            } else {
              # large region, report %conserved if aligned
              $showPos .= "$posFrom:$posTo";
              if ($isAligned) {
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
                $showPos .= " (vs. " . $alnPosToQpos{$site1->{posAlnFrom}}
                  . ":" . $alnPosToQpos{$site1->{posAlnTo}}
                    . ", ${percMatch}% identical)";
              }
            }
            $showPos = AddMouseOver($showPos, $nHit, $site1->{posAlnFrom})
              if $isAligned;
            my @siteDesc = map span({-id => "Site${nHit}S" . $_->{iSite} }, $_->{longDesc}), @sitesHere;
            push @bullets, li($showPos, join(br(), @siteDesc));
          } # end loop over PosTo
        } # end loop over PosFrom
      } else { # not PDB or SwissProt
        die "Unknown db $db";
      }
      print start_ul(), @bullets, end_ul(), "\n";
    } # end loop over isAligned
  } # end loop over hits
  my @pieces = $seq =~ /.{1,60}/g;
  print
    h3("Query Sequence"),
    p({-style => "font-family: monospace;"}, small(join(br(), ">" . HTML::Entities::encode($header), @pieces))),
    h3(a({-href => "sites.cgi"}, "New Search"));

  finish_page();
}

sub FormatAlnSites($$$$$$) {
  my ($queryAln, $hitAln, $queryBeg, $hitBeg, $hitName, $hposSite) = @_;
  my @out = ();
  my $queryAt = $queryBeg;
  my $hitAt = $hitBeg;
  die unless length($hitAln) == length($queryAln);
  for(my $i = 0; $i < length($hitAln); $i++) {
    my $queryChar = substr($queryAln, $i, 1);
    my $hitChar = substr($hitAln, $i, 1);
    my $string = "&nbsp;";
    if ($hitChar ne "-" && exists $hposSite->{$hitAt}) {
      my $sitesHere = $hposSite->{$hitAt};
      my $color = "black";

      my $hasSite1 = 0;
      if (exists $hposSite->{$hitAt}) {
        foreach my $site (@{ $hposSite->{$hitAt} }) {
          if ($site->{posFrom} eq $site->{posTo}) {
            $hasSite1 = 1;
          }
        }
      }
      my $agree;
      if ($queryChar eq $hitChar) {
        $string = "&vert;";
        $agree = 1;
      } else {
        $string = "x";
        $color = "darkred" if $hasSite1;
      }
      my $class = $hasSite1 ? ($agree ? "alnS1" : "alnS0") : "alnS";
      my $queryLong = exists $charToLong{$queryChar} ? $charToLong{$queryChar} : $queryChar;
      my $hitLong = exists $charToLong{$hitChar} ? $charToLong{$hitChar} : $hitChar;
      my $title = "${hitLong}${hitAt} in $hitName: "
        .  join("; ", map $_->{shortDesc}, @$sitesHere);
      $title .= " (${queryLong}${queryAt} in query)" if $queryChar ne "-";
      $string = a({-title => $title, -class => $class}, $string);
    }
    push @out, $string;
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

    if ($char eq "-") {
      push @out, "-";
    } else {
      my $class = exists $charToSet{$char} ? "aa" . $charToSet{$char} : "aa";
      my $longAA = exists $charToLong{$char} ? $charToLong{$char} : $char;
      my $title = "${longAA}${at} in $seqName";
      if (defined $charOther && $charOther ne "-") {
        my $longOther = exists $charToLong{$charOther} ? $charToLong{$charOther} : $char;
        $title .= " (${longOther}${atOther} in $nameOther)";
      }

      push @out, span({-title => $title, -class => $class }, $char);
      $at++;
    }
    $atOther++ if defined $charOther && $charOther ne "-";
  }
  return \@out;
}

sub AddMouseOver($$$) {
  my ($html, $nHit, $iAln) = @_;
  return $html unless $iAln >= 0;
  my $id = "Aln${nHit}P${iAln}";
  return span({ -onmouseover => qq{document.getElementById("$id").style.borderColor = "black"},
                -onmouseout => qq{document.getElementById("$id").style.borderColor = "$mygrey"} },
              $html);
}


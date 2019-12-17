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
    print p({-style => "margin-bottom: 0.5em;"},
            a({-title => $db eq "PDB" ? "entry $id chain $chain" : $id2,
               -href => $hitURL},
              $id.$chain),
            $desc,
            small($paperLink),
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
                 . FormatAlnString($alnQ, $queryBeg, {}, "query", $hitBeg, $alnS, $id.$chain),
                 i(a({-title => "$hitBeg to $hitEnd of $id$chain"}, "Sbjct:"))
                 . " "
                 . FormatAlnString($alnS, $hitBeg, \%sposToSite, $id.$chain, $queryBeg, $alnQ, "query"));
    print div({-style => "font-family: monospace; white-space:nowrap;" },
              p({-style => "margin-top: 0em; margin-bottom: 2em;"},
                join(br(), @lines)));
    print "\n";

    # First, remove the features that do not fully align, these will be rendered separately
    foreach my $site (@$sites) {
      $site->{isAligned} = $site->{posAlnFrom} >= 1 && $site->{posAlnFrom} <= $alnLen
        && $site->{posAlnTo} >= 1 && $site->{posAlnTo} <= $alnLen;
    }
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
          } elsif ($ligandId ne "NUC") {
            $ligShow = "binding " . a({-href => "http://www.rcsb.org/ligand/".$ligandId },
                                      $ligandId) . ":";
          } else {
            $ligShow = "binding DNA/RNA:";
          }
          my @posShow = ();
          foreach my $site (@$ligSites) {
            my $hitChar = substr($alnS, $site->{posAlnFrom}-1, 1);
            my $posShow = $hitChar . $site->{posFrom};
            $posShow = a({-title => "$site->{pdbFrom} in PDB numbering for $id$chain"},
                         $posShow)
              if $site->{pdbFrom} ne "";
            if ($site->{isAligned}) {
              my $posQ = $alnPosToQpos{ $site->{posAlnFrom} };
              $posShow .= " (vs. " . substr($alnQ, $site->{posAlnFrom}-1, 1) . $posQ . ")";
            }
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
                $showPos .= $alnPosToQpos{$site1->{posAlnFrom}};
                $showPos .= ":".$alnPosToQpos{$site1->{posAlnTo}}
                  if $posTo ne $posFrom;
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
            my @siteDesc = ();
            foreach my $site (@sitesHere) {
              my $type = $site->{type};
              my $comment = $site->{comment};
              my $desc = "$type $comment";
              if ($type eq "mutagenesis") {
                if ($comment =~ m/^([A-Z]+)->([A-Z]+): (.*)$/) {
                  if ($isAligned) {
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
                    if ($isAligned) {
                      $desc = "to $variant: $comment2";
                    } else {
                      $desc = "$original &rightarrow; $variant: $comment2";
                    }
                  }
                } else {
                  $desc = "modified: $comment";
                }
              }
              push @siteDesc, $desc;
            }
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

sub FormatAlnString($$$$) {
  my ($alnSeq, $beg, $posToSite, $seqName, $begOther, $alnOther, $nameOther) = @_;

  if (keys(%charToColor) == 0) {
    while (my ($color, $chars) = each %charColSets) {
      foreach my $char (split //, $chars) {
        $charToColor{$char} = $color;
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
      my @styleparts = ();
      push @styleparts, "background-color: " . $charToColor{$char} . ";"
        if exists $charToColor{$char};
      if (exists $posToSite->{$at}) {
        my $n = @{ $posToSite-> {$at} };
        push @styleparts, "font-weight: bold; text-decoration-line: underline; text-decoration-style: "
          . ($n > 1 ? "double" : "solid") . ";";
      }
      my $longAA = exists $charToLong{$char} ? $charToLong{$char} : $char;
      my $title = "${longAA}${at} in $seqName";
      if (defined $charOther && $charOther ne "-") {
        my $longOther = exists $charToLong{$charOther} ? $charToLong{$charOther} : $char;
        $title .= " (${longOther}${atOther} in $nameOther)";
      }

      push @out, a({-title => $title,
                    -style => join(" ",@styleparts)}, $char);
      $at++;
    }
    $atOther++ if defined $charOther && $charOther ne "-";
  }
  return join("", @out);
}

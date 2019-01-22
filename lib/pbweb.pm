# Utilities for PaperBLAST's web site
package pbweb;
use strict;
use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use Time::HiRes qw{gettimeofday};

our (@ISA,@EXPORT);
@ISA = qw(Exporter);
@EXPORT = qw(UniqToGenes SubjectToGene AddCuratedInfo GenesToHtml
             GetMotd FetchFasta HmmToFile loggerjs start_page finish_page);

# Returns a list of entries from SubjectToGene, 1 for each duplicate (if any),
# sorted by priority
sub UniqToGenes($$) {
  my ($dbh, $uniqId) = @_;
  my $dups = $dbh->selectcol_arrayref("SELECT duplicate_id FROM SeqToDuplicate WHERE sequence_id = ?",
                                      {}, $uniqId);
  my @subject_ids = ($uniqId);
  push @subject_ids, @$dups;

  my @genes = map { &SubjectToGene($dbh, $_) } @subject_ids;
  @genes = sort { $a->{priority} <=> $b->{priority} } @genes;
  return @genes;
}

# Given a row from the CuratedGene table, add fields for
# subjectId, source, curated, curateId, priority, URL, and showName,
# and clean up the comment field for presentation
sub AddCuratedInfo($) {
  my ($gene) = @_;
  my $db = $gene->{db};
  my $protId = $gene->{protId};
  die unless $db && $protId;
  $gene->{subjectId} = join("::", $db, $protId);
  $gene->{source} = $db;
  $gene->{curated} = 1;
  $gene->{curatedId} = $protId;

  if ($db eq "CAZy") {
    $gene->{source} = "CAZy via dbCAN";
    $gene->{URL} = "http://www.cazy.org/search?page=recherche&lang=en&recherche=$protId&tag=4";
    $gene->{priority} = 4;
  } elsif ($db eq "CharProtDB") {
    $gene->{priority} = 4;
    # their site is not useful, so just link to the paper
    $gene->{URL} = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3245046/";
    if ($gene->{comment}) {
      # label the comment as being from CharProtDB, as otherwise it is a bit mysterious.
      # And remove the Alias or Aliaess part.
      $gene->{comment} =~ s/Aliase?s?: [^ ;]+;? ?//;
      $gene->{comment} = i("CharProtDB") . " " . $gene->{comment};
    }
  } elsif ($db eq "SwissProt") {
    $gene->{URL} = "http://www.uniprot.org/uniprot/$protId";
    $gene->{priority} = 2;
    # and clean up the comments
    my @comments = split /_:::_/, $gene->{comment};
    @comments = map { s/[;. ]+$//; $_; } @comments;
    @comments = grep m/^SUBUNIT|FUNCTION|COFACTOR|CATALYTIC|ENZYME|DISRUPTION/, @comments;
    @comments = map {
      my @words = split / /, $_;
      my $word1 = $words[0];
      $words[0] = b(lc($words[0]));
      $words[1] = b(lc($words[1])) if @words > 1 && $words[1] =~ m/^[A-Z]+:$/;
      my $out = join(" ", @words);
      if ($word1 eq "COFACTOR:") {
        # Remove Evidence= and Xref= fields, as often found in the cofactor entry
        $out =~ s/ Evidence=[^;]*;?//g;
        $out =~ s/ Xref=[^;]*;?//g;
        # Transform Name=x; to x;
        $out =~ s/ Name=([^;]+);?/ $1/g;
        # Transform Note=note. to <small>(note.)</small>
        # require last char to be non-space to avoid " )" in output
        $out =~ s!Note=(.*)[.]!<small>($1.)</small>!g;
      } elsif ($word1 eq "CATALYTIC") {
        # Convert Xref=Rhea:RHEA:nnnn to a link to the Rhea entry, if it exists
        # Remove everthing else after the Reaction (i.e. EC: or ChEBI: entries)
        my $rheaId = $1 if $out =~ m/Xref=Rhea:RHEA:(\d+)/;
        $out =~ s/Reaction=([^;]+);.*$/$1/;
        $out .= " " . small("("
                      . a({ -href => "https://www.rhea-db.org/reaction?id=$rheaId" }, "RHEA:$rheaId")
                      . ")") if defined $rheaId;
      }
      $out;
    } @comments;
    my $comment = join("<BR>\n", @comments);
    $comment =~ s!{ECO:[A-Za-z0-9_:,.| -]+}!!g;
    $gene->{comment} = $comment;
  } elsif ($db eq "ecocyc") {
    $gene->{source} = "EcoCyc";
    $gene->{URL} = "https://ecocyc.org/gene?orgid=ECOLI&id=$protId";
    $gene->{priority} = 1;
  } elsif ($db eq "metacyc") {
    $gene->{source} = "MetaCyc";
    $gene->{URL} = "https://metacyc.org/gene?orgid=META&id=$protId";
    $gene->{priority} = 3;
  } elsif ($db eq "reanno") {
    $gene->{source} = "Fitness-based Reannotations";
    $gene->{comment} = "Mutant Phenotype: " . $gene->{comment};
    $gene->{priority} = 5;
    my ($orgId, $locusId) = split /:/, $protId;
    die "Invalid protId $protId" unless $locusId;
    $gene->{URL} = "http://fit.genomics.lbl.gov/cgi-bin/singleFit.cgi?orgId=$orgId&locusId=$locusId";
  } elsif ($db eq "REBASE") {
    $gene->{priority} = 4;
    $gene->{URL} = "http://rebase.neb.com/rebase/enz/$protId.html";
  } elsif ($db eq "BRENDA") {
    $gene->{priority} = 2.5; # just behind Swiss-Prot
    $gene->{source} = "BRENDA";
    $gene->{URL} = "http://www.brenda-enzymes.org/sequences.php?AC=" . $protId;
  } else {
    die "Unexpected curated database $db";
  }
  my @ids = ( $gene->{name}, $gene->{id2} );
  push @ids, $protId if $db eq "SwissProt";
  @ids = grep { $_ ne "" } @ids;
  $gene->{showName} = join(" / ", @ids) || $protId;
  $gene->{showName} = $protId if $db eq "REBASE";
}

# The returned entry will include:
# showName, URL, priority (for choosing what to show first), subjectId, desc, organism, protein_length, source,
# and other entries that depend on the type -- either papers for a list of GenePaper/PaperAccess items,
# or pmIds (a list of pubmed identifiers)
sub SubjectToGene($$) {
  my ($dbh, $subjectId) = @_;
  if ($subjectId =~ m/::/) { # curated gene
    my ($db, $protId) = split /::/, $subjectId;
    my $gene = $dbh->selectrow_hashref("SELECT * FROM CuratedGene WHERE db = ? AND protId = ?", {}, $db, $protId);
    die "Unrecognized subject $subjectId" unless defined $gene;
    AddCuratedInfo($gene);
    $gene->{pmIds} = $dbh->selectcol_arrayref("SELECT pmId FROM CuratedPaper WHERE db = ? AND protId = ?",
                                              {}, $db, $protId);
    return $gene;
  } else { # look in Gene table
    my $gene = $dbh->selectrow_hashref("SELECT * FROM Gene WHERE geneId = ?", {}, $subjectId);
    die "Unrecognized gene $subjectId" unless defined $gene;
    $gene->{subjectId} = $subjectId;
    $gene->{priority} = 6; # literature mined is lowest
    if ($subjectId =~ m/^VIMSS(\d+)$/) {
      my $locusId = $1;
      $gene->{source} = "MicrobesOnline";
      $gene->{URL} = "http://www.microbesonline.org/cgi-bin/fetchLocus.cgi?locus=$locusId";
    } elsif ($subjectId =~ m/^[A-Z]+_[0-9]+[.]\d+$/) { # refseq
      $gene->{URL} = "http://www.ncbi.nlm.nih.gov/protein/$subjectId";
      $gene->{source} = "RefSeq";
    } elsif ($subjectId =~ m/^[A-Z][A-Z0-9]+$/) { # SwissProt/TREMBL
      $gene->{URL} = "http://www.uniprot.org/uniprot/$subjectId";
      $gene->{source} = "SwissProt/TReMBL";
    } else {
      die "Cannot build a URL for subject $subjectId";
    }

    my $papers = $dbh->selectall_arrayref(qq{ SELECT DISTINCT * FROM GenePaper
                                              LEFT JOIN PaperAccess USING (pmcId,pmId)
                                              WHERE geneId = ?
                                              ORDER BY year DESC },
                                          { Slice => {} }, $subjectId);
    $gene->{papers} = $papers;

    # set up showName
    my @terms = map { $_->{queryTerm} } @$papers;
    my %terms = map { $_ => 1 } @terms;
    @terms = sort keys %terms;
    $gene->{showName} = join(", ", @terms) if !defined $gene->{showName};

    return $gene;
  }
}

my $li_with_style = qq{<LI style="list-style-type: none;" margin-left: 6em; >};
my $ul_with_style = qq{<UL style="margin-top: 0em; margin-bottom: 0em;">};

# Given the HTML for the coverage string, format the list of genes
sub GenesToHtml($$$$$) {
  my ($dbh, $uniqId, $genes, $coverage_html, $maxPapers) = @_;
  
  my @headers = ();
  my @content = ();
  my %paperSeen = (); # to avoid redundant papers -- pmId.pmcId.doi => term => 1
  my %paperSeenNoSnippet = (); # to avoid redundant papers -- pmId.pmcId.doi => 1
  my %seen_uniprot = (); # to avoid showing the curated Swiss-Prot entry and then the text-mined Swiss-Prot entry
  # (A metacyc entry could also mask a text-mined Swiss-Prot entry; that also seems ok)

  # Produce top-level and lower-level output for each gene (@headers, @content)
  # Suppress duplicate papers if no additional terms show up
  # (But, a paper could show up twice with two different terms, instead of the snippets
  # being merged...)
  foreach my $gene (@$genes) {
    die "No subjectId" unless $gene->{subjectId};
    $gene->{desc} = "No description" unless $gene->{desc}; # could be missing in MicrobesOnline or EcoCyc
    foreach my $field (qw{showName priority subjectId desc protein_length source}) {
      die "No $field for $gene->{subjectId}" unless $gene->{$field};
    }
    die "URL not set for $gene->{subjectId}" unless exists $gene->{URL};
    my $fromText = $gene->{organism} ? " from " . i($gene->{organism}) : "";
    my @pieces = ( a({ -href => $gene->{URL}, -title => $gene->{source},
                       -onmousedown => loggerjs("curated", $gene->{showName}) },
                     $gene->{showName}),
                   ($gene->{curated} ? b($gene->{desc}) : $gene->{desc}) . $fromText);
    # Formerly showed the coverage on the 1st line, but decided this was a bit wierd
    #push @pieces, $coverage_html if $gene == $genes->[0];
    # The alignment to show is always the one reported, not necessarily the one for this gene
    # (They are all identical, but only $subjectId is guaranteed to be in the blast database
    # and to be a valid argument for showAlign.cgi)
    if (exists $gene->{pmIds} && @{ $gene->{pmIds} } > 0) {
      my @pmIds = @{ $gene->{pmIds} };
      my %seen = ();
      @pmIds = grep { my $keep = !exists $seen{$_}; $seen{$_} = 1; $keep; } @pmIds;
      my $note = @pmIds > 1 ? scalar(@pmIds) . " papers" : "paper";
      push @pieces, "(see " .
        a({ -href => "http://www.ncbi.nlm.nih.gov/pubmed/" . join(",",@pmIds),
            -onmousedown => loggerjs("curatedpaper", $gene->{showName})},
          $note)
          . ")";
    }
    # For CAZy entries, add a link to the actual genbank entry because the CAZy entry is a bit mysterious
    if ($gene->{source} =~ m/^CAZy/) {
      my $id = $gene->{showName};
      $id =~ s/[.]\d+$//;
      if ($id =~ m/^[A-Z0-9_]+/) {
        push @pieces, "(see " .
          a({ -href => "https://www.ncbi.nlm.nih.gov/protein/$id",
              -title => "NCBI protein entry",
              -onmousedown => loggerjs("cazygenbank", $gene->{showName}) },
            "protein")
            . ")";
      }
    }
    # Skip the header if this is a UniProt entry that is redundant with a curated
    # (Swiss-Prot) entry
    push @headers, join(" ", @pieces)
      unless exists $seen_uniprot{$gene->{showName}} && !exists $gene->{curated};
    $seen_uniprot{ $gene->{curatedId} } = 1
      if exists $gene->{curatedId};
    
    push @content, $gene->{comment} if $gene->{comment};
    my $nPaperShow = 0;
    foreach my $paper (@{ $gene->{papers} }) {
      my @pieces = (); # what to say about this paper
      my $snippets = [];
      $snippets = $dbh->selectall_arrayref(
          "SELECT DISTINCT * from Snippet WHERE geneId = ? AND pmcId = ? AND pmId = ?",
          { Slice => {} },
          $gene->{subjectId}, $paper->{pmcId}, $paper->{pmId})
        if $paper->{pmcId} || $paper->{pmId};
      
      my $paperId = join(":::", $paper->{pmId}, $paper->{pmcId}, $paper->{doi});
      my $nSkip = 0; # number of duplicate snippets
      foreach my $snippet (@$snippets) {
        my $text = $snippet->{snippet};
        # In case of XML or HTML tags slipping into the snippet (which is rare)
        $text =~ s!<!&lt;!g;
        $text =~ s!/>!/&gt;!g;
        my $term = $snippet->{queryTerm};
        if (exists $paperSeen{$paperId}{$term}) {
          $nSkip++;
        } else {
          $text =~ s!($term)!<B><span style="color: red;">$1</span></B>!gi;
          push @pieces, "&ldquo;...$text...&rdquo;";
        }
      }
      # ignore this paper if all snippets were duplicate terms
      next if $nSkip == scalar(@$snippets) && $nSkip > 0;
      $nPaperShow++;
      if ($nPaperShow > $maxPapers) {
        push @content, a({-href => "litSearch.cgi?more=".$uniqId},
                         "More");
        last;
        next;
      }
      foreach my $snippet (@$snippets) {
        my $term = $snippet->{queryTerm};
        $paperSeen{$paperId}{$term} = 1;
      }
      
      # Add RIFs
      my $rifs = [];
      $rifs = $dbh->selectall_arrayref(qq{ SELECT DISTINCT * from GeneRIF
                                                        WHERE geneId = ? AND pmcId = ? AND pmId = ? },
                                       { Slice => {} },
                                       $gene->{subjectId}, $paper->{pmcId}, $paper->{pmId})
        if $paper->{pmcId} || $paper->{pmId};
      my $GeneRIF_def = a({ -title => "from Gene Reference into Function (NCBI)",
                            -href => "https://www.ncbi.nlm.nih.gov/gene/about-generif",
                            -style => "color: black; text-decoration: none; font-style: italic;" },
                          "GeneRIF");
      # just 1 snippet if has a GeneRIF
      pop @pieces if @$rifs > 0 && @pieces > 1;
      foreach my $rif (@$rifs) {
        # normally there is just one
        unshift @pieces, $GeneRIF_def . ": " . $rif->{ comment };
      }
      
      my $paper_url = undef;
      my $pubmed_url = "http://www.ncbi.nlm.nih.gov/pubmed/" . $paper->{pmId};
      if ($paper->{pmcId} && $paper->{pmcId} =~ m/^PMC\d+$/) {
        $paper_url = "http://www.ncbi.nlm.nih.gov/pmc/articles/" . $paper->{pmcId};
      } elsif ($paper->{pmid}) {
        $paper_url = $pubmed_url;
      } elsif ($paper->{doi}) {
        if ($paper->{doi} =~ m/^http/) {
          $paper_url = $paper->{doi};
        } else {
          $paper_url = "http://doi.org/" . $paper->{doi};
        }
      }
      my $title = $paper->{title};
      $title = a({-href => $paper_url, -onmousedown => loggerjs("pb", $gene->{showName})}, $title)
        if defined $paper_url;
      my $authorShort = $paper->{authors};
      $authorShort =~ s/ .*//;
      my $extra = "";
      $extra = "(" . a({ -href => $pubmed_url, -onmousedown => loggerjs("pb", $gene->{showName}) }, "PubMed") . ")"
        if !$paper->{pmcId} && $paper->{pmId};
      my $paper_header = $title . br() .
        small( a({ -title => $paper->{authors} }, "$authorShort,"),
               $paper->{journal}, $paper->{year}, $extra);
      
      if (@pieces == 0) {
        # Skip if printed already for this gene (with no snippet)
        next if exists $paperSeenNoSnippet{$paperId};
        $paperSeenNoSnippet{$paperId} = 1;
        
        # Explain why there is no snippet
        my $excuse;
        my $short;
        if (!defined $paper->{access}) {
          ;
        } elsif ($paper->{access} eq "full") {
          $short = "no snippet";
          $excuse = "This term was not found in the full text, sorry.";
        } elsif ($paper->{isOpen} == 1) {
          if ($paper->{access} eq "abstract") {
            $short = "no snippet";
            $excuse = "This paper is open access but PaperBLAST only searched the the abstract.";
          } else {
            $short = "no snippet";
            $excuse = "This paper is open access but PaperBLAST did not search either the full text or the abstract.";
          }
        } elsif ($paper->{isOpen} eq "") {
          # this happens if the link is from GeneRIF
          $short = "no snippet";
          $excuse = "PaperBLAST did not search either the full text or the abstract.";
        } elsif ($paper->{journal} eq "") {
          $short = "secret";
          $excuse = "PaperBLAST does not have access to this paper, sorry";
        } else {
          $short = "secret";
          $excuse = "$paper->{journal} is not open access, sorry";
        }
        if ($excuse) {
          
          my $href = a({-title => $excuse}, $short);
          $paper_header .= " " . small("(" . $href . ")"); 
        }
      }
      my $pieces = join($li_with_style, @pieces);
      $pieces = join("", $ul_with_style, $li_with_style, $pieces, "</UL>")
        if $pieces;
      push @content, $paper_header . $pieces;
    }
  }
  my $content = join($li_with_style, @content);
  $content = join("", $ul_with_style, $li_with_style, $content, "</UL>")
    if $content;
  push @headers, $coverage_html;
  return p({-style => "margin-top: 1em; margin-bottom: 0em;"},
           join("<BR>", @headers) . $content) . "\n";
}

# type, proteinid => onclick javascript code for a link related to the specified protein
sub loggerjs($$) {
  my ($type, $prot) = @_;
  my $string = $type . "::" . $prot;
  return qq{logger(this, '$string')};
}

sub GetMotd {
  my $motd = "";
  if (open(MOTD, "<", "../motd")) {
    $motd = join("\n", <MOTD>);
    close(MOTD);
    $motd =~ s/\r//g;
    $motd =~ s/\s+$//;
  }
  $motd = p($motd) if $motd ne "";
  return $motd;
}

# Given an accession, from either uniq.faa *or* the duplicate_id field of SeqToDuplicate,
# return its sequence
sub FetchFasta($$$) {
    my ($dbh, $db, $acc) = @_;
    my $tmpDir = "../tmp";
    my $procId = $$;
    my $timestamp = int (gettimeofday() * 1000);
    my $prefix = "$tmpDir/$procId$timestamp";

    # First, check if the sequence is a duplicate
    my $uniqIds = $dbh->selectcol_arrayref("SELECT sequence_id FROM SeqToDuplicate WHERE duplicate_id = ?",
                                           {}, $acc);
    my $acc2 = $acc;
    $acc2 = $uniqIds->[0] if @$uniqIds > 0;

    die "Invalid def2 argument: $acc2" if $acc2 eq "" || $acc2 =~ m/\s/ || $acc2 =~ m/,/;
    my $fastacmd = "../bin/blast/fastacmd";
    die "No such executable: $fastacmd" unless -x $fastacmd;
    system("$fastacmd","-s",$acc2,"-d",$db,"-o", "$prefix.fetch");
    open(SEQ, "<", "$prefix.fetch") || die "Cannot read $prefix.fetch -- fastacmd failed?";
    my @lines = <SEQ>;
    close(SEQ) || die "Error reading $prefix.fetch";
    unlink("$prefix.fetch");
    (@lines > 0 && $lines[0] =~ m/^>/) || die "Unknown accession: $acc";
    shift @lines;
    @lines = map { chomp; $_; } @lines;
    return join("", @lines);
}

sub HmmToFile($) {
  my ($hmmId) = @_;
  return undef unless $hmmId;
  if ($hmmId && $hmmId =~ m/^[a-zA-Z0-9_.-]+$/) {
    my @hmmdir = ("../static/pfam", "../static/tigrfam", "../tmp");
    foreach my $hmmdir (@hmmdir) {
      return "$hmmdir/$hmmId.hmm" if -e "$hmmdir/$hmmId.hmm";
    }
  }
  # if reached
  my @glob = glob("../static/pfam/$hmmId.*.hmm");
  return @glob > 0 ? $glob[0] : undef;
}

sub TopDivHtml($$) {
  my ($banner, $URL) = @_;
  return <<END
<div style="background-color: #40C0CB; display: block; position: absolute; top: 0px; left: -1px;
  width: 100%; padding: 0.25em; z-index: 400;">
<H2 style="margin: 0em;">
<A HREF="$URL" style="color: gold; font-family: 'Montserrat', sans-serif; font-style:italic;
  text-shadow: 1px 1px 1px #000000; text-decoration: none;">
$banner
</A></H2></div>
<P style="margin: 0em;">&nbsp;</P>
<SCRIPT src="../static/pb.js"></SCRIPT>
END
;
}

# Run with parameters like 'title' => 'Curated BLAST for Genomes'
# Optional parameters:
# 'banner' (defaults to the PaperBLAST banner; otherwise specify arbitrary html)
# 'bannerURL' (the URL that the banner links to)
sub start_page {
  my %param = @_;
  my $title = $param{title} || die "Must specify title";
  my $banner = $param{banner} || "PaperBLAST &ndash; <small>Find papers about a protein or its homologs</small>";
  my $bannerURL = $param{bannerURL} || "litSearch.cgi";
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

  # utf-8 because that is the encoding used by EuropePMC
  print
    header(-charset => 'utf-8'),
    start_html(-head => Link({-rel => "shortcut icon", -href => "../static/favicon.ico"}),
               -style => { -code => $style },
             -script => [{ -type => "text/javascript", -src => "../static/autocomplete_uniprot.js" }],
             -title => $title),
    &TopDivHtml($banner, $bannerURL),
    h2($title),
    "\n";
}

# This should be called only after start_page()
sub finish_page() {
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

# Utilities for PaperBLAST's web site
package pbweb;
use strict;
use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use Time::HiRes qw{gettimeofday};
use LWP::Simple qw{get};
use JSON;
use DBI;
use IO::String;
use Bio::SeqIO;
use URI::Escape;

our (@ISA,@EXPORT);
@ISA = qw(Exporter);
@EXPORT = qw(UniqToGenes SubjectToGene AddCuratedInfo
             GeneToHtmlLine GenesToHtml
             GetMotd FetchFasta HmmToFile loggerjs start_page finish_page
             VIMSSToFasta RefSeqToFasta UniProtToFasta FBrowseToFasta DBToFasta pdbToFasta
             commify
             LinkifyComment FormatStepPart DataForStepParts HMMToURL
             parseSequenceQuery sequenceToHeader
             warning fail
);

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
    $gene->{comment} = b("mutant phenotype:") . " " . $gene->{comment};
    $gene->{priority} = 5;
    my ($orgId, $locusId) = split /:/, $protId;
    die "Invalid protId $protId" unless $locusId;
    $gene->{URL} = "http://fit.genomics.lbl.gov/cgi-bin/domains.cgi?orgId=$orgId&locusId=$locusId";
  } elsif ($db eq "REBASE") {
    $gene->{priority} = 4;
    $gene->{URL} = "http://rebase.neb.com/rebase/enz/$protId.html";
  } elsif ($db eq "BRENDA") {
    $gene->{priority} = 2.5; # just behind Swiss-Prot
    $gene->{source} = "BRENDA";
    $gene->{URL} = "http://www.brenda-enzymes.org/sequences.php?AC=" . $protId;
  } elsif ($db eq "UniProt") { # used in GapMind, i.e. see curatedClusters.cgi
    $gene->{priority} = 6;
    $gene->{source} = "UniProt";
    $gene->{URL} = "http://www.uniprot.org/uniprot/$protId";
  } elsif ($db eq "TCDB") {
    $gene->{priority} = 2.8; # just behind Brenda
    $gene->{source} = "TCDB";
    if ($gene->{id2} =~ m/,/) { # in more than one system
      $gene->{URL} = "http://www.tcdb.org/search/result.php?acc=$protId";
    } else {
      $gene->{URL} = "http://www.tcdb.org/search/result.php?tc=" . $gene->{id2}
        . "&acc=" . $protId;
    }
    my @comments = split /_:::_/, $gene->{comment};
    @comments = map { s/[;. ]+$//; $_; } @comments;
    my @out = ();
    foreach my $comment (@comments) {
     my @words = split / /, $comment;
     $words[0] = "TCDB comment:" if $words[0] eq "COMMENT:";
     $words[0] = b(lc($words[0]));
     push @out, join(" ", @words);
   }
    @out = sort @out if $db eq "TCDB"; # substrates before comments
    $gene->{comment} = join("<BR>\n", @out);
  } elsif ($db eq "biolip") {
    $gene->{priority} = 5; # behind most, but ahead of UniProt entries
    my $entry = $protId; # i.e., 101mA
    $entry =~ s/[A-Z]$//;
    $gene->{URL} = "https://www.rcsb.org/structure/" . uc($entry);
    # to show which ligand goes with which PDB structure
    $gene->{comment} .= " " . "(" . a({ -href => $gene->{URL} }, $protId) . ")";
  } else {
    die "Unexpected curated database $db";
  }
  my @ids = ( $gene->{name}, $gene->{id2} );
  push @ids, $protId if $db eq "SwissProt";
  unshift @ids, $protId if $db eq "biolip"; # always show pdb entry/chain first
  if ($db eq "TCDB") {
    my @tcids = split /,/, $gene->{id2};
    @ids = map { "TC $_" } @tcids;
    push @ids, $protId;
  }
  @ids = grep { defined $_ && $_ ne "" } @ids;
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
  }
  # else look in Gene table or HasSites
  my $gene = $dbh->selectrow_hashref("SELECT * FROM Gene WHERE geneId = ?", {}, $subjectId);
  if ($gene) {
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
      $gene->{URL} = "https://www.uniprot.org/uniprot/$subjectId";
      $gene->{source} = "SwissProt/TReMBL";
    } else {
      die "Cannot build a URL for subject $subjectId";
    }
  } elsif ($subjectId =~ m/^[a-zA-Z]+:/) {
    # Look in HasSites. This only works for ids like SwissProt:P90754 or PDB:1t0b:A
    my @parts = split /:/, $subjectId;
    die "Invalid subject $subjectId is not in Gene or HasSite"
      unless @parts >= 2 && @parts <= 3;
    my ($db, $id, $chain) = @parts;
    $chain = "" if !defined $chain;
    $gene = $dbh->selectrow_hashref("SELECT * FROM HasSites WHERE db = ? AND id = ? AND chain = ?",
                                    {}, $db, $id, $chain);
    die "Cannot find $subjectId in Gene or HasSites table" unless $gene;
    $gene->{subjectId} = $subjectId;
    $gene->{source} = $db;
    $gene->{protein_length} = $gene->{length};
    if ($db eq "SwissProt") {
      $gene->{URL} = "https://www.uniprot.org/uniprot/$id";
    } elsif ($db eq "PDB") {
      $gene->{URL} = "https://rcsb.org/structure/$id";
    } else {
      die "Cannot build a URL for subject $subjectId -- unhandled db $db";
    }
    $gene->{priority} = 4.5; # above literature mined
  } else {
    die "Cannot handle subject $subjectId not in Gene table";
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
  if (!defined $gene->{showName}) {
    $gene->{showName} = join(", ", @terms);
    if ($gene->{showName} eq "") {
      # for HasSites entries, there won't be any papers
      $gene->{showName} = $subjectId;
      $gene->{showName} =~ s/^[^:]+://;
      $gene->{showName} =~ s/://g;
    }
  }
  return $gene;
}

my $li_with_style = qq{<LI style="list-style-type: none;" margin-left: 6em; >};
my $ul_with_style = qq{<UL style="margin-top: 0em; margin-bottom: 0em;">};

sub GeneToHtmlLine($) {
  my ($gene) = @_;

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
      a({ -href => "http://pubmed.ncbi.nlm.nih.gov/pubmed/?term=" . join(",",@pmIds). "&sort=pubdate",,
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
  return join(" ", @pieces);
}

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
    my $line = GeneToHtmlLine($gene);
    my (@pieces) = $line;

    # Skip the header if this is a UniProt entry that is redundant with a curated
    # (Swiss-Prot) entry
    push @headers, $line
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
    # also handle file names like PF11902.8.hmm from hmmId PF11902
    my @glob = glob("../static/pfam/$hmmId.*.hmm");
    return @glob > 0 ? $glob[0] : undef;
  }
  return undef;
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
/* for showing protein sequences */
/* base style, and for unusual items */
.aa  { background-color: white; }
.aa1 { background-color: lightblue; } /* for AILMFWV */
.aa2 { background-color: red; } /* for KR */
.aa3 { background-color: magenta; } /* for ED */
.aa4 { background-color: lightgreen; } /* for NQST */
.aa5 { background-color: pink; } /* for C */
.aa6 { background-color: orange; } /* for G */
.aa7 { background-color: yellow; } /* for P */
.aa8 { background-color: cyan; } /* for HY */

/* for showing the alignment line */
.alnS, .alnS1, .alnS0 { font-size: smaller; }
.alnS1 { font-weight: bold; color: darkgreen; } /* site, matching */
.alnS0 { font-weight: bold; color: darkred; } /* site, not matching */

/* for an alignment column */
.alnCol, .alnLabel {
  vertical-align:top; float:left; display:inline-block;
  margin-top: 1em;
  font-family: "Courier New", monospace;
}
.alnCol { background-color: #EEEEEE; border: solid 1px #EEEEEE; }
/* In many browsers, the label prevents the columns from wrapping all the way to the left,
   so, make it subtly smaller */
.alnLabel { margin-right: 0.5em; padding-top: 0.04em; border: none; font-size: 95% }

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


# Given a locus tag or VIMSSnnnn query, get it in FASTA format
sub VIMSSToFasta($) {
  my ($short) = @_;
  die unless defined $short;
  my $mo_dbh = DBI->connect('DBI:mysql:genomics:pub.microbesonline.org', "guest", "guest")
    || die $DBI::errstr;
  my $locusId;
  if ($short =~ m/^VIMSS(\d+)$/i) {
    $locusId = $1;
  } else {
    # try to find the locus tag
    ($locusId) = $mo_dbh->selectrow_array( qq{SELECT locusId FROM Synonym JOIN Locus USING (locusId,version)
						WHERE name = ? AND priority = 1 },
                                              {}, $short );
  }
  return undef unless $locusId;
  my ($aaseq) = $mo_dbh->selectrow_array( qq{SELECT sequence FROM Locus JOIN AASeq USING (locusId,version)
                                             WHERE locusId = ? AND priority=1 },
                                          {}, $locusId);

  unless (defined $aaseq) {
    print p("$short matches",
            a({-href => "http://www.microbesonline.org/cgi-bin/fetchLocus.cgi?locus=$locusId"},
              "VIMS$locusId"),
            "from MicrobesOnline, which is not a protein-coding gene.");
    return undef;
  }
  my ($desc) = $mo_dbh->selectrow_array( qq{SELECT description FROM Locus JOIN Description USING (locusId,version)
                                             WHERE locusId = ? AND priority=1 },
                                          {}, $locusId);
  return ">$short $desc\n$aaseq\n" if $desc;
  return ">$short\n$aaseq\n" if $desc;
}

sub RefSeqToFasta($) {
  my ($short) = @_;
  die unless defined $short;
  if ($short =~ m/^[A-Z]P_\d+[.]?\d*$/ || $short =~ m/^[A-Z][A-Z][A-Z]\d+[.]?\d*$/) {
    # Potentially an NCBI protein identifier like WP_093840703.1 or AAC76544.1
    # (the version number like .1 is optional)
    my $url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.cgi?db=Protein&rettype=fasta&id=$short";
    my $results = get($url);
    return $results if defined $results && $results =~ m/^>/;
  }
  return undef unless $short =~ m/^[A-Za-z][A-Za-z0-9]+_[A-Za-z0-9]+[.]?\d?$/;

  my $url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.cgi?db=Nucleotide&retmode=json&term=$short";
  my $json = from_json(get($url));
  return undef unless defined $json;
  my $id = $json->{esearchresult}{idlist}[0];
  return undef unless $id;

  # Fetch genbank format for entry $id
  print "<P>Looking for $short in genbank entry $id\n";
  $url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=Nucleotide&rettype=gbwithparts&retmode=text&id=$id";
  my $gb = get($url);
  my $seqio = Bio::SeqIO->new(-fh => IO::String->new($gb), -format => "genbank");
  while (my $seq = $seqio->next_seq) {
    foreach my $ft ($seq->get_SeqFeatures) {
      next unless $ft->primary_tag eq "CDS"
        && $ft->has_tag("translation")
        && (($ft->has_tag("locus_tag") && ($ft->get_tag_values("locus_tag"))[0] eq $short)
            || ($ft->has_tag("old_locus_tag") && ($ft->get_tag_values("old_locus_tag"))[0] eq $short));
      my ($aaseq) = $ft->get_tag_values("translation");
      my $defline = $short;
      $defline .= " " . ($ft->get_tag_values("protein_id"))[0] if $ft->has_tag("protein_id");
      $defline .= " " . ($ft->get_tag_values("product"))[0] if $ft->has_tag("product");
      $defline .= " [" . $seq->desc . "]";
      $defline =~ s/[\t\r\n]+/ /g; # not sure if this can ever occur
      return ">$defline\n$aaseq";
    }
  }
  return undef;
}

sub UniProtToFasta($) {
  my ($short) = @_;
  die unless defined $short;
  # include=no -- no isoforms
  my $url = "http://www.uniprot.org/uniprot/?query=${short}&format=fasta&sort=score&include=no&limit=2";
  my $results = get($url);
  if (defined $results && $results =~ m/^>/) {
    # select the first hit only
    my @lines = split /\n/, $results;
    my @out = ();
    my $nHeader = 0;
    foreach my $line (@lines) {
      $nHeader++ if substr($line, 0, 1) eq ">";
      push @out, $line if $nHeader <= 1;
    }
    return join("\n", @out)."\n";
  }
  # else
  return undef;
}

# The first argument is the Fitness Browser data directory
sub FBrowseToFasta($$) {
  my ($fbdata, $short) = @_;
  return undef unless -d $fbdata;
  my $fastacmd = "../bin/blast/fastacmd";
  die "No such executable: $fastacmd" unless -x $fastacmd;

  my $fbdbh = DBI->connect("dbi:SQLite:dbname=$fbdata/feba.db","","",{ RaiseError => 1 }) || die $DBI::errstr;
  my $gene = $fbdbh->selectrow_hashref("SELECT * FROM Gene WHERE locusId = ? OR sysName = ? OR gene = ? LIMIT 1",
                                       {}, $short, $short, $short);
  my $fasta;
  if ($gene) {
    my $seqId = $gene->{orgId} . ":" . $gene->{locusId};
    die "Missing aaseqs file for fitness browser: $fbdata/aaseqs\n"
      unless -e "$fbdata/aaseqs";

    my $procId = $$;
    my $timestamp = int (gettimeofday() * 1000);
    my $filename = $procId . $timestamp;
    my $seqFile = "/tmp/$filename.fasta";

    if (system($fastacmd, "-s", $seqId, "-o", $seqFile, "-d", "$fbdata/aaseqs") == 0) {
      open(SEQ, "<", $seqFile) || die "Cannot read $seqFile";
      my $seq = "";
      while (my $line = <SEQ>) {
        next if $line =~ m/^>/;
        chomp $line;
        die "Invalid output from fastacmd" unless $line =~ m/^[A-Z*]+$/;
        $seq .= $line;
      }
      close(SEQ) || die "Error reading $seqFile";
      unlink($seqFile);

      my $showId = $gene->{sysName} || $gene->{locusId};
      my $org = $fbdbh->selectrow_hashref("SELECT * FROM Organism WHERE orgId = ?",
                                          {}, $gene->{orgId})
        || die "No such organism in $fbdata/feba.db: $gene->{orgId}";
      my $orgdesc = join(" ", $org->{genus}, $org->{species}, $org->{strain});
      $fasta = ">$showId $gene->{gene} $gene->{desc} ($orgdesc)\n$seq\n";
    }
  }
  $fbdbh->disconnect();
  return $fasta;
}

# Given an identifier that might be from pdb, returns a string in fasta format, or undef
# Supported identifiers are of the form 3osd, 2eq7C, 2eq7_C, or 2eq7_2 (case insensitive)
sub pdbToFasta($) {
  my ($query) = @_;
  $query = uc($query);
  return undef unless $query =~ m/^([0-9][0-9A-Z][0-9A-Z][0-9A-Z])([A-Z]?_?[A-Z]?\d*)$/i;
  my ($entry, $chainSpec) = ($1,$2);
  return undef unless $chainSpec eq ""
    || $chainSpec =~ m/_?[A-Z]/i
    || $chainSpec =~ m/_\d+/;
  $chainSpec =~ s/_//g; # now it should be a number or a chain letter
  my $URL = "https://www.rcsb.org/fasta/entry/$entry/display"; # either case works
  my $pdbFasta = get($URL);
  return undef unless $pdbFasta =~ m/>/;
  my @lines = split /\n/, $pdbFasta;
  my %seqs = ();
  my $id;
  foreach my $line (@lines) {
    if ($line =~ m/^>(.*)$/) {
      $id = $1;
      fail("Duplicate fasta for $id from " . a({-href => $URL}, "RCSB"). " -- $line")
        if exists $seqs{$id};
      $seqs{$id} = "";
    } else {
      fail("Missing header from " . a({-href => $URL}, "RCSB"))
        unless defined $id;
      $seqs{$id} .= $line;
    }
  }
  my @ids = keys %seqs;
  if (@ids == 1 && $chainSpec eq "") {
    $id = $ids[0];
  } else {
    my %chain = (); # chain specifier such as A, B, or 1, to id
    foreach my $id (@ids) {
      # Example ids:
      # 2EQ7_1|Chains A, B|2-oxoglutarate dehydrogenase E3 component|Thermus thermophilus (300852)
      # 3OSD_1|Chain A|putative glycosyl hydrolase|Bacteroides thetaiotaomicron (226186)
      my ($spec, $chains) = split /[|]/, $id;
      $spec =~ s/^.*_//;
      $chain{$spec} = $id;
      $chains =~ s/^chains? //i;
      $chains =~ s/, //g;
      foreach my $chain (split //, $chains) {
        $chain{uc($chain)} = $id;
      }
    }
    $chainSpec = "A" if $chainSpec eq "";
    $id = $chain{$chainSpec};
    fail("Fasta entry " . a({-href => $URL}, $entry) . " does not have chain $chainSpec")
      unless defined $id;
  }
  return ">$id\n$seqs{$id}\n";
}

sub DBToFasta($$$) {
  my ($dbh, $blastdb, $query) = @_;
  my ($gene, $geneId);
  $gene = $dbh->selectrow_hashref("SELECT * from Gene WHERE geneId = ?", {}, $query);
  $geneId = $gene->{geneId} if $gene;
  if (! $gene) {
    $gene = $dbh->selectrow_hashref("SELECT * from CuratedGene WHERE db = 'SwissProt' AND protId = ?",
                                    {}, $query);
    $geneId = "SwissProt::".$query if $gene;
  }
  if (! $gene && $query =~ m/^(.*)::(.*)$/) {
    my ($db,$protId) = ($1,$2);
    $gene = $dbh->selectrow_hashref("SELECT * from CuratedGene WHERE db = ? AND protId = ?",
                                    {}, $db, $protId);
    $geneId = "${db}::${protId}" if $gene;
  }

  if (defined $geneId) {
    # look for the duplicate
    my $desc = $gene->{desc};
    my $org = $gene->{organism};
    my $dup = $dbh->selectrow_hashref("SELECT * from SeqToDuplicate WHERE duplicate_id = ?", {}, $geneId);
    my $seqId = $dup ? $dup->{sequence_id} : $geneId;
    my $fastacmd = "../bin/blast/fastacmd";
    die "No such executable: $fastacmd" unless -x $fastacmd;
    # Note some genes are dropped from the database during construction so it
    # may fail to find it

    my $procId = $$;
    my $timestamp = int (gettimeofday() * 1000);
    my $filename = $procId . $timestamp;
    my $seqFile = "/tmp/$filename.fasta";

    if (system($fastacmd, "-s", $seqId, "-o", $seqFile, "-d", $blastdb) == 0) {
      open(SEQ, "<", $seqFile) || die "Cannot read $seqFile";
      my $seq = "";
      while (my $line = <SEQ>) {
        next if $line =~ m/^>/;
        chomp $line;
        die "Invalid output from fastacmd" unless $line =~ m/^[A-Z*]+$/;
        $seq .= $line;
      }
      close(SEQ) || die "Error reading $seqFile";
      return ">$geneId $desc ($org)\n$seq\n";
    }
  }
  return undef;
}

sub commify($) {
    local $_  = shift;
    1 while s/^(-?\d+)(\d{3})/$1,$2/;
    return $_;
}

# Given text string, such as the comment in a GapMind steps file, replace terms
# with links. Examples of terms that are replaced:
# pfam:PF14574
# uniprot:Q8AAD7_BACTN
# PMID:7523119
# PMC179677
# metacyc:PHENOLDEG-PWY
# URL:https://mspace.lib.umanitoba.ca/handle/1993/30289
# EC:2.8.3.6
sub LinkifyComment($) {
  my ($comment) = @_;
  my @words = split /\s/, $comment;
  my @out = ();
  foreach my $word (@words) {
    # Pull of leading brackets or parentheses
    my $pre = "";
    $pre = $1 if $word =~ m/^([\[({]+)/;
    $word =~ s/^([\[({]+)//;
    if ($word =~ m/^pfam:(PF\d+)/i) {
      my $pfam = $1;
      $word =~ s/^pfam:PF\d+//i;
      push @out, $pre . a({ -href => "hmmSearch.cgi?hmmId=$pfam" }, $pfam) . $word;
    } elsif ($word =~ m/^uniprot:([A-Z][A-Z0-9_]+)/i) {
      my $uniprotId = $1;
      $word =~ s/^uniprot:[A-Z][A-Z0-9_]+//i;
      push @out, $pre . a({ -href => "https://www.uniprot.org/uniprot/$uniprotId" }, $uniprotId) . $word;
    } elsif ($word =~ m/^pmid:(\d+)/i) {
      my $pmId = $1;
      $word =~ s/^pmid:\d+//i;
      push @out, $pre . a({ -href => "https://pubmed.ncbi.nlm.nih.gov/$pmId/" }, "PMID:$pmId") . $word;
    } elsif ($word =~ m/^(PMC\d+)/) {
      my $pmcId = $1;
      $word =~ s/^PMC\d+//;
      push @out, $pre . a({ -href => "http://www.ncbi.nlm.nih.gov/pmc/articles/" . lc($pmcId) . "/" },
                   $pmcId) . $word;
    } elsif ($word =~ m/^metacyc:([0-9A-Z][A-Z0-9-]+)/i) {
      my $metacycId = $1;
      $word =~ s/^metacyc:([0-9A-Z][A-Z0-9-]+)//i;
      push @out, $pre. a({ -href => "https://metacyc.org/META/NEW-IMAGE?object=$metacycId" },
                         "link") . $word;
    } elsif ($word =~ m!^URL:(http[A-Za-z0-9_,:./?&-]+)!i) {
      my $URL = $1;
      $word =~ s!^URL:(http[A-Za-z0-9_,:./?&-]+)!!;
      push @out, $pre. a({ -href => $URL }, "link") . $word;
    } elsif ($word =~ m/^EC:([0-9][.][0-9.]+[0-9])/i) {
      my $ec = $1;
      $word =~ s/^EC:([0-9][.][0-9.]+[0-9])//i;
      push @out, $pre . "EC " . a({ -href => "https://enzyme.expasy.org/EC/$ec" }, $ec) . $word;
    } else {
      push @out, $pre . $word;
    }
  }
  return join(" ", @out);
}

# Fetch the sequence queries (uniprotQuery and curatedQuery) for a step
sub DataForStepParts($$$) {
  my ($dbhS, $pathwayId, $stepId) = @_;
  my $stepParts = $dbhS->selectall_arrayref("SELECT * from StepPart WHERE pathwayId = ? AND stepId = ?",
                                            { Slice => {} }, $pathwayId, $stepId);
  my $stepQueries = $dbhS->selectall_arrayref("SELECT * from StepQuery WHERE pathwayId = ? AND stepId = ?",
                                              { Slice => {} }, $pathwayId, $stepId);
  my %curatedQuery = (); # curatedId (not ids -- a single compoennt) to stepquery row
  # (Includes entries of type curated or ignore)
  my %uniprotQuery = (); # uniprotId to stepquery row
  foreach my $sq (@$stepQueries) {
    if ($sq->{queryType} eq "curated" || $sq->{queryType} eq "ignore") {
      $sq->{desc} =~ s/;;/. /g;
      foreach my $id (split /,/, $sq->{curatedIds}) {
        $curatedQuery{$id} = $sq;
      }
    } elsif ($sq->{queryType} eq "uniprot") {
      $uniprotQuery{$sq->{uniprotId}} = $sq;
    }
  }
  return { 'uniprotQuery' => \%uniprotQuery,
           'curatedQuery' => \%curatedQuery };
}

# Given the data, the step part row (from the database), the set, and the orgId (or ""),
# and the hash of organisms, return HTML for the step part
sub FormatStepPart($$$$$) {
  my ($data, $stepPart, $set, $orgId, $orgs) = @_;
  die unless defined $data->{curatedQuery} && defined $data->{uniprotQuery};
  my $type = $stepPart->{partType};
  my $value = $stepPart->{value};

  if ($type eq "EC") {
    # Use local URLs for Curated BLAST links, instead of using the official papers.genomics.lbl.gov
    # site, because the genome may not exist at the public site
    my $out = "Curated proteins or TIGRFams with EC "
      . a({-href => "https://enzyme.expasy.org/EC/$value" }, $value);
    $out .= " (" . a({ -href => "genomeSearch.cgi?gdb="
                       . $orgs->{$orgId}{gdb}
                       . "&gid=" . $orgs->{$orgId}{gid}
                       . "&query=$value",
                       -title => "Run Curated BLAST" },
                     "search")
      . ")" if $orgId ne "";
    return $out;
  } elsif ($type eq "hmm") {
    return "HMM " . a({-href => HMMToURL($value) }, $value);
  } elsif ($type eq "term") {
    # nead to uri_escape the query because it may contain %
    my $URL = "curatedClusters.cgi?set=$set&word=1&query=" . uri_escape($value);
    $URL = "genomeSearch.cgi?gdb=" . $orgs->{$orgId}{gdb}
          . "&gid=" . $orgs->{$orgId}{gid}
          . "&word=1"
          . "&query=" . uri_escape($value)
            if $orgId ne "";
    return "Curated proteins matching "
      . a({ -href => $URL,
            -title => $orgId eq "" ? "" : "Run Curated BLAST" }, $value);
  } elsif ($type eq "curated") {
    my $URL = "http://papers.genomics.lbl.gov/cgi-bin/litSearch.cgi?query=".$value;
    my $show_id = $value; $show_id =~ s/^.*://;
    # Find the relevant step query
    return "Curated sequence " . a({-href => $URL, -title => "View in PaperBLAST"}, $show_id)
      . ": " . $data->{curatedQuery}{$value}{desc};
  } elsif ($type eq "uniprot") {
    my $URL = "https://www.uniprot.org/uniprot/".$value;
    return "UniProt sequence " . a({-href => $URL, -title => "View in UniProt"}, $value)
      . ": " . $data->{uniprotQuery}{$value}{desc};
  } elsif ($type eq "ignore_other") {
    my $URL = "http://papers.genomics.lbl.gov/cgi-bin/curatedSearch.cgi?word=1"
      . "&query=" . uri_escape($value); # value could contain %
    return "Ignore hits to items matching "
      . a({-href => $URL}, $value)
      . " when looking for 'other' hits";
  } elsif ($type eq "ignore") {
    my $URL = "http://papers.genomics.lbl.gov/cgi-bin/litSearch.cgi?query=$value";
    my $showId = $value; $showId =~ s/^.*://;
    my $desc = $data->{curatedQuery}{$value}{desc} || "entry not in GapMind's database";
    return  join(" ",
                 "Ignore hits to ",
                 a({-href => $URL, -title => "View in PaperBLAST"}, $showId),
                 "when looking for 'other' hits ($desc)");
  } elsif ($type eq "ignore_hmm") {
    return "Do not include HMM " . a({ -href => HMMToURL($value) }, $value)
      . " (when considering EC numbers)";
  }
  die "Unknown StepPart type $type";
}

sub HMMToURL($) {
  my ($hmmId) = @_;
  if ($hmmId =~ m/^TIGR/) {
    return "https://www.ncbi.nlm.nih.gov/Structure/cdd/".$hmmId;
  } elsif ($hmmId =~ m/^PF/) {
    my $hmmIdShort = $hmmId; $hmmIdShort =~ s/[.]\d+$//;
    return "http://pfam.xfam.org/family/".$hmmIdShort;
  }
  return "";
}

sub warning {
  my @warnings = @_;
  my $URL = CGI::url(-relative => 1);
  if ($URL ne "") {
    print p({ -style => "color: red;" }, @warnings), "\n";
  } else {
    print join(" ", @warnings)."\n";
  }
}

sub fail($) {
  my ($notice) = @_;
  my $URL = CGI::url(-relative => 1);
  if ($URL ne "") {
    print
      p($notice),
      p(a({-href => $URL}, "New query")),
      end_html;
  } else {
    print "Failed: $notice\n";
  }
  exit(0);
}

# Takes as arguments as hash that includes
# -query -- the input
# -dbh -- the PaperBLAST database handle
# -blastdb -- the filename for PaperBLAST's BLAST database (usually uniq.faa)
# -fbdata -- the filename for the fitness browser database (optional)
#
# The input should be in fasta format, raw sequence (no definition line, linebreaks allowed),
# uniprot format, or be an identifier from various databases.
#
# Any * characters in the input sequence are replaced with X.
#
# On success, returns (definition line, sequence)
#   where the definition line will be empty if the
#   query was sequence only, with no definition line.
# If the input is empty (or all whitespace), returns ().
# If it does not recognize the identifier, or it's not valid input, it uses fail().
sub parseSequenceQuery {
  my (%param) = @_;
  foreach my $field (qw{-query -dbh -blastdb}) {
    die "parameter $field must be provided to parseSequenceQuery()"
      unless defined defined $param{$field};
  }
  my $query = $param{-query};
  my $dbh = $param{-dbh};
  my $blastdb = $param{-blastdb};
  my $fbdata = $param{-fbdata};

  # remove leading and trailing whitespace
  $query =~ s/^\s+//;
  $query =~ s/\s+$//;
  return if $query eq "";

  # a single word query is assumed to be a gene id if it contains any non-sequence character
  # But, putting a protein sequence on a line is allowed (if all uppercase)
  if ($query ne "" && $query !~ m/\n/ && $query !~ m/ / && $query =~ m/[^A-Z*]/) {
    my $short = $query;
    $query = undef;
    fail("Sorry, query has a FASTA header but no sequence") if $short =~ m/^>/;

    # Is it a VIMSS id?
    $query = &VIMSSToFasta($short) if $short =~ m/^VIMSS\d+$/i;

    # Is it in the database?
    if (!defined $query) {
      $query = &DBToFasta($dbh, $blastdb, $short);
    }

    # Is it a PDB identifier
    if (!defined $query) {
      $query = &pdbToFasta($short);
    }

    # is it a fitness browser locus tag?
    if (!defined $query && $fbdata && $short =~ m/^[0-9a-zA-Z_]+$/) {
      $query = &FBrowseToFasta($fbdata, $short);
    }

    # is it a UniProt id or gene name or protein name?
    if (!defined $query) {
      $query = &UniProtToFasta($short);
    }

    # is it in VIMSS as a locus tag or other synonym?
    if (!defined $query) {
      $query = &VIMSSToFasta($short);
    }

    # is it in Nucleotide/RefSeq? (Locus tags not in refseq may not be indexed)
    if (!defined $query) {
      $query = &RefSeqToFasta($short);
    }

    my $shortSafe = HTML::Entities::encode($short);
    &fail("Sorry -- we were not able to find a protein sequence for the identifier <b>$shortSafe</b>. We checked it against our database of proteins that are linked to papers, against UniProt (including their ID mapping service), against MicrobesOnline, against the NCBI protein database (RefSeq and Genbank), and against PDB/RCSB. Please use the sequence as a query instead.")
      if !defined $query;
  }

  my $seq = "";
  my $def = "";
  my @lines = split /[\r\n]+/, $query;
  if (@lines > 0 && $lines[0] =~ m/^>/) {
    $def = shift @lines;
    $def =~ s/^>//;
  }
  foreach (@lines) {
    s/[ \t]//g;
    s/^[0-9]+//; # leading digit/whitespace occurs in UniProt format
    next if $_ eq "//";
    &fail("Error: more than one sequence was entered.") if m/^>/;
    &fail("Unrecognized characters in sequence")
      unless m/^[a-zA-Z*]*$/;
    s/[*]/X/g;
    $seq .= uc($_);
  }

  my $seqlen = length($seq);
  fail("Sequence is too short") unless length($seq) >= 10;
  my @nt = $seq =~ m/[ACGTUN]/g;
  my $fACGTUN = scalar(@nt) / $seqlen;
  if ($fACGTUN >= 0.9) {
    warning(sprintf("Warning: sequence is %.1f%% nucleotide characters -- are you sure this is a protein query?",
                    100 * $fACGTUN));
  }
  return ($def, $seq);
}

sub sequenceToHeader($) {
  my ($seq) = @_;
  my $seqlen = length($seq);
  my $initial = substr($seq, 0, 10);
  $initial .= "..." if $seqlen > 10;
  return length($seq) . " a.a. (" . $initial . ")";
}

1

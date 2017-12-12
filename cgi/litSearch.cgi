#!/usr/bin/perl -w
#######################################################
## litSearch.cgi
##
## Copyright (c) 2017 University of California
##
## Authors: Morgan Price
#######################################################
#
# Optional CGI garameters:
# query -- this should be the protein sequence in FASTA or UniProt or plain format,
#	or a VIMSS id, or a gene identifier in the database,
#	or a locus tag in MicrobesOnline,
#	or a UniProt id or other gene name that UniProt recognizes
#	or a Genbank protein (or RefSeq protein) identifier.
# more -- which geneId (in the database) to show the full list of papers for
#	(not to be used with query)
#
# If none of these is specified, shows the query box, an example link, and some documentation

use strict;
use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use Time::HiRes qw{gettimeofday};
use DBI;
use LWP::Simple qw{get};
use IO::Handle; # for autoflush

sub fail($);
sub simstring($$$$$$$$$$$$$);
sub SubjectToGene($);
sub commify($);
sub loggerjs($$); # type, proteinid => onclick javascript code for a link related to the specified protein

my $tmpDir = "../tmp";
my $blastall = "../bin/blast/blastall";
my $nCPU = 4;
my $base = "../data";
my $blastdb = "$base/uniq.faa";
my $sqldb = "$base/litsearch.db";
die "No such executable: $blastall" unless -x $blastall;
die "No such file: $blastdb" unless -e $blastdb;
die "No such file: $sqldb" unless -e $sqldb;

# If a gene has more papers than this, there is a "more" lnik to show all the information
# for that gene on a separate page

my $cgi=CGI->new;
my $query = $cgi->param('query') || "";
my $more_subjectId = $cgi->param('more') || "";
my $maxPapers = $more_subjectId ? 10000 : 8;

my $dbh = DBI->connect("dbi:SQLite:dbname=$sqldb","","",{ RaiseError => 1 }) || die $DBI::errstr;

my %stats = ();
open(STATS, "<", "$base/stats") || die "Cannot read $base/stats";
while(my $line = <STATS>) {
  chomp $line;
  my ($key,$value) = split /\t/, $line, -1;
  $stats{$key} = $value;
}
close(STATS) || die "Error reading $base/stats";
foreach my $key (qw{nSeq nPaper}) {
  $stats{$key} = commify( $stats{$key} );
}

my $documentation = <<END

<H3><A NAME="stats">Statistics</A></H3>

The PaperBLAST database links $stats{nSeq} different protein sequences to $stats{nPaper} scientific articles. Searches against EuropePMC were last performed on $stats{date}.

<H3><A NAME="works">How It Works</A></H3>

<P>PaperBLAST builds a database of protein sequences that are linked
to scientific articles. These links come from automated text searches
against the articles in <A HREF="http://europepmc.org/">EuropePMC</A>
and from manually-curated information from <A
HREF="https://www.ncbi.nlm.nih.gov/gene/about-generif" title="Gene
Reference into Function (NCBI)">GeneRIF</A>, <A
HREF="http://www.uniprot.org/">Swiss-Prot</A>,
<A HREF="http://www.cazy.org/" title="Carbohydrate-Active enZYmes Database">CAZy</A> (as made available by <A HREF="http://csbl.bmb.uga.edu/dbCAN/download.php">dbCAN</A>),
<A HREF="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3245046/" title="Database of experimentally characterized protein annotations">CharProtDB</A>,
<A HREF="http://metacyc.org/">MetaCyc</A>,
<A HREF="http://ecocyc.org">EcoCyc</A>,
<A HREF="http://rebase.neb.com/rebase/rebase.html" title="The Restriction Enzyme Database">REBASE</A>,
and the <A HREF="http://fit.genomics.lbl.gov/" title="Reannotations from genome-wide fitness data">Fitness Browser</A>.
Given this database and a protein sequence query,
PaperBLAST uses <A
HREF="https://en.wikipedia.org/wiki/BLAST">protein-protein BLAST</A>
to find similar sequences with E &lt; 0.001.

<P>To build the database, we query EuropePMC with locus tags, with <A
HREF="https://www.ncbi.nlm.nih.gov/refseq/">RefSeq</A> protein
identifiers, and with <A HREF="http://www.uniprot.org/">UniProt</A>
accessions.  We obtain the locus tags from RefSeq or from <A
HREF="http://www.microbesonline.org/">MicrobesOnline</A>.  We use
queries of the form "locus_tag AND genus_name" to try to ensure that
the paper is actually discussing that gene.  Because EuropePMC indexes
most recent biomedical papers, even if they are not open access, some
of the links may be to papers that you cannot read or that our
computers cannot read.  We query each of these identifiers that
appears in the open access part of EuropePMC, as well as every locus
tag that appears in the 500 most-referenced genomes, so that a gene
may appear in the PaperBLAST results even though none of the papers
that mention it are open access. We also incorporate text mined links
from EuropePMC that link open access articles to UniProt or RefSeq
identifiers.  (This yields some additional links because EuropePMC
uses different heuristics for their text mining than we do.)

<P>For every article that mentions a locus tag, a RefSeq protein
identifier, or a UniProt accession, we try to select one or two
snippets of text that refer to the protein. If we cannot get access to
the full text, we try to select a snippet from the abstract, but
unfortunately, unique identifiers such as locus tags are rarely
provided in abstracts.

<P>PaperBLAST also incorporates manually-curated links between protein sequences and
articles: <UL> <LI>Proteins from NCBI's RefSeq are included if a
<A HREF="https://www.ncbi.nlm.nih.gov/gene/about-generif">GeneRIF</A>
entry links the gene to an article in <A
HREF="http://www.ncbi.nlm.nih.gov/pubmed/">PubMed</A><sup>&reg;</sup>.
GeneRIF also provides a short summary of the article's claim about the
protein, which is shown instead of a snippet.  <LI>Proteins
from Swiss-Prot (the curated part of <A HREF="http://uniprot.org">UniProt</A>) are included if the curators
identified experimental evidence for the protein's function (evidence
code ECO:0000269). For these proteins, the fields of the Swiss-Prot entry that
describe the protein's function are shown (with bold headings).
  <LI>Every protein from <A HREF="http://ecocyc.org">EcoCyc</A>, a curated
database of the proteins in <i> Escherichia coli</i> K-12, is included, regardless
of whether they are characterized or not.
<LI>Proteins from the <A HREF="http://metacyc.org">MetaCyc</A> metabolic pathway database are included if they are linked to a paper in PubMed.
<LI>Every protein from <A HREF="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3245046/">CharProtDB</A>, a database of experimentally characterized protein annotations, is included.
<LI>Proteins from the <A HREF="http://www.cazy.org/">CAZy</A> database of carbohydrate-active enzymes are included if they are associated with an Enzyme Classification number. Even though CAZy does not provide links from individual protein sequences to papers, these should all be experimentally-characterized proteins.
<LI>Proteins from the <A HREF="http://rebase.neb.com/rebase/rebase.html">REBASE</A> database of restriction enzymes if they have known specificity.
<LI>Every protein with an evidence-based reannotation (based on mutant phenotypes) in the <A HREF="http://fit.genomics.lbl.gov/">Fitness Browser</A> is included.
  </UL>Except for GeneRIF,
the curated entries include a short curated
description of the protein's function.  Many of these entries also
link to articles in <A
HREF="http://www.ncbi.nlm.nih.gov/pubmed/">PubMed</A>.

<P>For more information see the <A title="PaperBLAST: Text Mining Papers for Information about Homologs" HREF="http://msystems.asm.org/content/2/4/e00039-17">PaperBLAST paper</A> (<i>mSystems</i> 2017) or the
<A
HREF="https://github.com/morgannprice/PaperBLAST">code</A>.

Also note some changes since the paper was written:

<UL>
<LI>December 2017: incorporated MetaCyc, CharProtDB, CAZy, REBASE, and the reannotations from the Fitness Browser.
<LI>September 2017: EuropePMC no longer returns some table entries in their search results. This has shrunk PaperBLAST's database, but has also reduced the number of low-relevance hits.
</UL>

<H3><A NAME="secret">Secrets</A></H3>

<P>PaperBLAST cannot provide snippets for many of the papers that are
published in non-open-access journals. This limitation applies even if
the paper is marked as "free" on the publisher's web site and is
available in PubmedCentral or EuropePMC. If a journal that you publish
in is marked as "secret," please consider publishing elsewhere.

<H3><A NAME="omission">Omissions from the PaperBLAST Database</A></H3>

<P>Some important articles are missing from PaperBLAST, either because
the article's full text is not in EuropePMC (as for many older
articles) or because of PaperBLAST's heuristics. If you notice an
article that characterizes a protein's function but is missing from
PaperBLAST, please notify the curators at <A
HREF="http://www.uniprot.org/update">UniProt</A> or add an entry to <A
HREF="https://www.ncbi.nlm.nih.gov/gene/submit-generif">GeneRIF</A>.
Entries in either of these databases will eventually be incorporated
into PaperBLAST.  Note that to add an entry to UniProt, you will need
find the UniProt identifier for the protein.  If the protein is not
already in UniProt, you can ask them to create an entry.  To add an
entry to GeneRIF, you will need an NCBI Gene identifier, but
unfortunately many prokaryotic proteins in RefSeq do not have
corresponding Gene identifers.

<center>by <A HREF="http://morgannprice.org/">Morgan Price</A>,
<A HREF="http://genomics.lbl.gov/">Arkin group</A><BR>
Lawrence Berkeley National Laboratory
</center>
END
    ;

my $title = "PaperBLAST";
# utf-8 because that is the encoding used by EuropePMC
print
    header(-charset => 'utf-8'),
    start_html(-head => Link({-rel => "shortcut icon", -href => "../static/favicon.ico"}),
               -title => $title),
    qq{<div style="background-color: #40C0CB; display: block; position: absolute; top: 0px; left: -1px; width: 100%; padding: 0.25em; z-index: 400;"><H2 style="margin: 0em;"><A HREF="litSearch.cgi" style="color: gold; font-family: 'Montserrat', sans-serif; font-style:italic; text-shadow: 1px 1px 1px #000000; text-decoration: none;">PaperBLAST &ndash; <small>Find papers about a protein or its homologs</small></A></H2></div><P style="margin: 0em;">&nbsp;</P>\n},
    qq{<SCRIPT src="../static/pb.js"></SCRIPT>\n};

my $procId = $$;
my $timestamp = int (gettimeofday() * 1000);
my $filename = $procId . $timestamp;
my $seqFile = "$tmpDir/$filename.fasta";

# remove leading and trailing whitespace
$query =~ s/^\s+//;
$query =~ s/\s+$//;
# a single word query is assumed to be a gene id if it contains any non-sequence character
# But, putting a protein sequence on a line is allowed (if all uppercase)

if ($query ne "" && $query !~ m/\n/ && $query !~ m/ / && $query =~ m/[^A-Z*]/) {
  my $short = $query;
  $query = undef;
  fail("Sorry, query has a FASTA header but no sequence") if $short =~ m/^>/;

  # Is it a VIMSS id?
  $query = &VIMSSToQuery($short) if $short =~ m/^VIMSS\d+$/i;

  # Is it in the database?
  if (!defined $query) {
    my $gene = $dbh->selectrow_hashref("SELECT * from Gene WHERE geneId = ?", {}, $short);
    my $geneId;
    if ($gene) {
      $geneId = $gene->{geneId};
    } else {
      $gene = $dbh->selectrow_hashref("SELECT * from CuratedGene WHERE db = 'SwissProt' AND protId = ?",
                                      {}, $short);
      $geneId = "SwissProt::".$short;
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
      if (system($fastacmd, "-s", $seqId, "-o", $seqFile, "-d", $blastdb) == 0) {
        open(SEQ, "<", $seqFile) || die "Cannot read $seqFile";
        my $seq = "";
        while (my $line = <SEQ>) {
          next if $line =~ m/^>/;
          chomp $line;
          die "Invalid output from fastacmd" unless $line =~ m/^[A-Z*]+$/;
          $seq .= $line;
        }
        $query = ">$geneId $desc ($org)\n$seq\n";
      }
    }
  }

  # is it a UniProt id or gene name or protein name?
  if (!defined $query) {
    $query = &UniProtToQuery($short);
  }

  # is it in VIMSS as a locus tag or other synonym?
  if (!defined $query) {
    $query = &VIMSSToQuery($short);
  }

  # is it in RefSeq?
  if (!defined $query) {
    $query = &RefSeqToQuery($short);
  }

  &fail("Sorry -- we were not able to find a protein sequence for the identifier <b>$short</b>. We checked it against our database of proteins that are linked to papers, against UniProt (including their ID mapping service), against MicrobesOnline, and against the NCBI protein database (RefSeq and Genbank). Please use the sequence as a query instead.")
    if !defined $query;
}

my $seq;
my $def = "";
my $hasDef = 0;
if ($query =~ m/[A-Za-z]/) {
    $seq = "";
    my @lines = split /[\r\n]+/, $query;
    $def = shift @lines if @lines > 0 && $lines[0] =~ m/^>/;
    $def =~ s/^>//;
    foreach (@lines) {
        s/[ \t]//g;
        s/^[0-9]+//; # leading digit/whitespace occurs in UniProt format
        next if $_ eq "//";
        &fail("Error: more than one sequence was entered.") if m/^>/;
        &fail("Unrecognized characters in $_") unless m/^[a-zA-Z*]*$/;
        s/[*]/X/g;
        $seq .= uc($_);
    }
    $hasDef = 1 if $def ne "";
    $def = length($seq) . " a.a." if $def eq "";
}

my $motd = "";
if (open(MOTD, "<", "../motd")) {
  $motd = join("\n", <MOTD>);
  close(MOTD);
  $motd =~ s/\r//g;
  $motd =~ s/\s+$//;
}
$motd = p($motd) if $motd ne "";

if (!defined $seq && ! $more_subjectId) {
    my $exampleId = "3615187";
    my $refseqId = "WP_012018426.1";
    print
        $motd,
        start_form( -name => 'input', -method => 'GET', -action => 'litSearch.cgi'),
        p(br(),
          b("Enter a protein sequence in FASTA or Uniprot format,<BR>or an identifier from UniProt, RefSeq, or MicrobesOnline: "),
          br(),
          textarea( -name  => 'query', -value => '', -cols  => 70, -rows  => 10 )),
        p(submit('Search'), reset()),
        end_form,
        p("Or see",
          a({ -href => "litSearch.cgi?query=VIMSS$exampleId" }, "example results"),
            "for the putative alcohol dehydrogenase",
          a({ -href => "https://www.ncbi.nlm.nih.gov/protein/$refseqId",
                      -style => "color: black;" },
                    small("$refseqId,") ),
            "which is actually the regulator",
            a({ -href => "litSearch.cgi?query=VIMSS$exampleId",
                -title => "Show PaperBLAST hits" },
              i("ercA")) . "."),
          $documentation;
} else {
    if ($more_subjectId) {
      print h3("Full List of Papers Linked to $more_subjectId");
    } else {
      die "No sequence to search" unless $seq;
      my $initial = substr($seq, 0, 10);
      my $seqlen = length($seq);
      $initial .= "..." if $seqlen > 10;
      $initial = "$seqlen a.a., $initial" if $hasDef;
      print
        $motd,
        h3("PaperBLAST Hits for $def ($initial)");

      my @nt = $seq =~ m/[ACGTUN]/g;
      my $fACGTUN = scalar(@nt) / $seqlen;
      if ($fACGTUN >= 0.9) {
        printf("<P><font color='red'>Warning: sequence is %.1f%% nucleotide characters -- are you sure this is a protein query?</font>",
               100 * $fACGTUN);
      }
    }

    autoflush STDOUT 1; # show preliminary results
    print "\n";

    my @hits = ();
    if ($more_subjectId) {
      push @hits, [ $more_subjectId, $more_subjectId, 100.0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ];
    } else {
      open(SEQ, ">", $seqFile) || die "Cannot write to $seqFile";
      print SEQ ">$def\n$seq\n";
      close(SEQ) || die "Error writing to $seqFile";
      my $hitsFile = "$tmpDir/$filename.hits";
      system($blastall, "-p", "blastp", "-d", $blastdb, "-i", $seqFile, "-o", $hitsFile,
             "-e", 0.001, "-m", 8, "-a", $nCPU, "-F", "m S") == 0 || die "Error running blastall: $!";
      open(HITS, "<", $hitsFile) || die "Cannot read $hitsFile";
      while(<HITS>) {
        chomp;
        my @F = split /\t/, $_;
        push @hits, \@F;
      }
      close(HITS) || die "Error reading $hitsFile";
      unlink($seqFile);
      unlink($hitsFile);
    }
    my $nHits = scalar(@hits);
    if ($nHits == 0) {
        print p("Sorry, no hits to proteins in the literature.");
    } else {
        print p("Found $nHits similar proteins in the literature:"), "\n"
          unless $more_subjectId;

        my %seen_subject = ();
        my $li_with_style = qq{<LI style="list-style-type: none;" margin-left: 6em; >};
        my $ul_with_style = qq{<UL style="margin-top: 0em; margin-bottom: 0em;">};
        foreach my $row (@hits) {
            my ($queryId,$subjectId,$percIdentity,$alnLength,$mmCnt,$gapCnt,$queryStart,$queryEnd,$subjectStart,$subjectEnd,$eVal,$bitscore) = @$row;
            next if exists $seen_subject{$subjectId};
            $seen_subject{$subjectId} = 1;

            my $dups = $dbh->selectcol_arrayref("SELECT duplicate_id FROM SeqToDuplicate WHERE sequence_id = ?",
                                                {}, $subjectId);
            my @subject_ids = ($subjectId);
            push @subject_ids, @$dups;

            my @genes = map { &SubjectToGene($_) } @subject_ids;
            @genes = sort { $a->{priority} <=> $b->{priority} } @genes;
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
            foreach my $gene (@genes) {
                die "No subjectId" unless $gene->{subjectId};
                $gene->{desc} = "No description" unless $gene->{desc}; # could be missing in MicrobesOnline or EcoCyc
                foreach my $field (qw{showName priority subjectId desc protein_length source}) {
                    die "No $field for $subjectId" unless $gene->{$field};
                }
                die "URL not set for $subjectId" unless exists $gene->{URL};
                my $fromText = $gene->{organism} ? " from " . i($gene->{organism}) : "";
                my @pieces = ( a({ -href => $gene->{URL}, -title => $gene->{source},
                                   -onmousedown => loggerjs("curated", $gene->{showName}) },
                                 $gene->{showName}),
                               ($gene->{curated} ? b($gene->{desc}) : $gene->{desc}) . $fromText);
                # The alignment to show is always the one reported, not necessarily the one for this gene
                # (They are all identical, but only $subjectId is guaranteed to be in the blast database
                # and to be a valid argument for showAlign.cgi)
                push @pieces, &simstring(length($seq), $gene->{protein_length},
                                         $queryStart,$queryEnd,$subjectStart,$subjectEnd,
                                         $percIdentity,$eVal,$bitscore,
                                         $def, $gene->{showName}, $seq, $subjectId)
                    if $gene->{subjectId} eq $genes[0]{subjectId} && ! $more_subjectId;
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
                      push @content, a({-href => "litSearch.cgi?more=".$subjectId},
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
            print p({-style => "margin-top: 1em; margin-bottom: 0em;"},
                    join("<BR>", @headers) . $content) . "\n";
        }
    }

    print qq{<script src="http://fit.genomics.lbl.gov/d3js/d3.min.js"></script>
             <script src="http://fit.genomics.lbl.gov/images/fitblast.js"></script>
             <H3><A title="Fitness BLAST searches for similarity to bacterial proteins that have mutant phenotypes" HREF="http://fit.genomics.lbl.gov/" NAME="#fitness">Fitness Blast Results</A></H3>
             <P><DIV ID="fitblast_short"></DIV></P>
             <script>
             var server_root = "http://fit.genomics.lbl.gov/";
             var seq = "$seq";
             fitblast_load_short("fitblast_short", server_root, seq);
             </script>
    } unless $more_subjectId;

    if (! $more_subjectId) {
      my @pieces = $seq =~ /.{1,60}/g;
      print h3("Query Sequence"),
        p({-style => "font-family: monospace;"}, small(join(br(), ">$def", @pieces)));
    }
    print h3(a({-href => "litSearch.cgi"}, "New Search")),
      $documentation,
        end_html;
}

sub fail($) {
    my ($notice) = @_;
    print
        p($notice),
        p(a({-href => "litSearch.cgi"}, "New search")),
        end_html;
    exit(0);
}

sub simstring($$$$$$$$$$$$$) {
    my ($qLen, $sLen, $queryStart,$queryEnd,$subjectStart,$subjectEnd,$percIdentity,$eVal,$bitscore,
        $def1, $def2, $seq1, $acc2) = @_;
    $percIdentity = sprintf("%.0f", $percIdentity);
    # the minimum of coverage either way
    my $cov = ($queryEnd-$queryStart+1) / ($qLen > $sLen ? $qLen : $sLen);
    my $percentCov = sprintf("%.0f", 100 * $cov);
    my $title ="$queryStart:$queryEnd/$qLen of query is similar to $subjectStart:$subjectEnd/$sLen of hit (E = $eVal, $bitscore bits)";
    return "(" .
        a({ -title => $title,
            -href => "showAlign.cgi?def1=$def1&def2=$def2&seq1=$seq1&acc2=$acc2" },
          "$percIdentity% identity, $percentCov% coverage")
        . ")";
}

# The entry will include:
# showName, URL, priority (for choosing what to show first), subjectId, desc, organism, protein_length, source,
# and other entries that depend on the type -- either papers for a list of GenePaper/PaperAccess items,
# or pmIds (a list of pubmed identifiers)
sub SubjectToGene($) {
  my ($subjectId) = @_;
  if ($subjectId =~ m/::/) { # curated gene
    my ($db, $protId) = split /::/, $subjectId;
    my $gene = $dbh->selectrow_hashref("SELECT * FROM CuratedGene WHERE db = ? AND protId = ?", {}, $db, $protId);
    die "Unrecognized subject $subjectId" unless defined $gene;
    $gene->{subjectId} = $subjectId;
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
        my $cofactor = $words[0] eq "COFACTOR:";
        $words[0] = b(lc($words[0]));
        $words[1] = b(lc($words[1])) if @words > 1 && $words[1] =~ m/^[A-Z]+:$/;
        my $out = join(" ", @words);
        if ($cofactor) {
          # Remove Evidence= and Xref= fields., as often found in the cofactor entry
          $out =~ s/ Evidence=[^ ]*;?//g;
          $out =~ s/ Xref=[^ ]+;?//g;
          # Transform Name=x; to x;
          $out =~ s/ Name=([^;]+);/ $1;/g;
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
    } else {
      die "Unexpeced database $db";
    }

    my @ids = ( $gene->{name}, $gene->{id2} );
    push @ids, $protId if $db eq "SwissProt";
    @ids = grep { $_ ne "" } @ids;
    $gene->{showName} = join(" / ", @ids) || $protId;
    $gene->{showName} = $protId if $db eq "REBASE";
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

# Given a locus tag or VIMSSnnnn query, get it in FASTA format
sub VIMSSToQuery($) {
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

  &fail("Sorry, VIMSS$locusId is not a protein in MicrobesOnline") unless defined $aaseq;
  my ($desc) = $mo_dbh->selectrow_array( qq{SELECT description FROM Locus JOIN Description USING (locusId,version)
                                             WHERE locusId = ? AND priority=1 },
                                          {}, $locusId);
  return ">$short $desc\n$aaseq\n" if $desc;
  return ">$short\n$aaseq\n" if $desc;
}

sub RefSeqToQuery($) {
  my ($short) = @_;
  die unless defined $short;
  return undef unless $short =~ m/_/;
  my $url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.cgi?db=Protein&rettype=fasta&id=$short";
  my $results = get($url);
  return $results if $results =~ m/^>/;
  return undef;
}

sub UniProtToQuery($) {
  my ($short) = @_;
  die unless defined $short;
  # include=no -- no isoforms
  my $url = "http://www.uniprot.org/uniprot/?query=${short}&format=fasta&sort=score&include=no&limit=2";
  my $results = get($url);
  if ($results =~ m/^>/) {
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

sub commify($) {
    local $_  = shift;
    1 while s/^(-?\d+)(\d{3})/$1,$2/;
    return $_;
}

sub loggerjs($$) {
  my ($type, $prot) = @_;
  my $string = $type . "::" . $prot;
  return qq{logger(this, '$string')};
}


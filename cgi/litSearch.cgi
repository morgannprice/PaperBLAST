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
use CGI qw(:standard Vars start_ul);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use Time::HiRes qw{gettimeofday};
use DBI;
use LWP::Simple qw{get};
use HTML::Entities;
use IO::Handle; # for autoflush
use lib "../lib";
use pbweb; # for SubjectToGene()

sub fail($);
sub simstring($$$$$$$$$$$$$);
sub commify($);

my $tmpDir = "../tmp";
my $blastall = "../bin/blast/blastall";
my $nCPU = 4;
my $base = "../data";
my $blastdb = "$base/uniq.faa";
my $sqldb = "$base/litsearch.db";
my $fastacmd = "../bin/blast/fastacmd";
die "No such executable: $blastall" unless -x $blastall;
die "No such executable: $fastacmd" unless -x $fastacmd;
die "No such file: $blastdb" unless -e $blastdb;
die "No such file: $sqldb" unless -e $sqldb;

# A symbolic link to the Fitness Browser data directory is used (if it exists)
# to allow quick access to proteins from the fitness browser.
# That directory must include feba.db (sqlite3 database) and aaseqs (in fasta format)
my $fbdata = "../fbrowse_data"; # path relative to the cgi directory

# If a gene has more papers than this, there is a "more" link to show all the information
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
against the articles in <A HREF="http://europepmc.org/" title="Europe PMC (PubMedCentral)">EuropePMC</A>
and from manually-curated information from <A
HREF="https://www.ncbi.nlm.nih.gov/gene/about-generif" title="Gene
Reference into Function (NCBI)">GeneRIF</A>, <A
HREF="http://www.uniprot.org/" title="The manually annotated and reviewed section of UniProt">UniProtKB/Swiss-Prot</A>,
<A HREF="http://www.brenda-enzymes.org/index.php" title="The Comprehensive Enzyme Information System">BRENDA</A>,
<A HREF="http://www.cazy.org/" title="Carbohydrate-Active enZYmes Database">CAZy</A> (as made available by <A HREF="http://csbl.bmb.uga.edu/dbCAN/download.php">dbCAN</A>),
<A HREF="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3245046/" title="Database of experimentally characterized protein annotations">CharProtDB</A>,
<A HREF="http://metacyc.org/" title="MetaCyc Metabolic Pathway Database">MetaCyc</A>,
<A HREF="http://ecocyc.org" title="EcoCyc: Encyclopedia of E. coli Genes and Metabolic Pathways">EcoCyc</A>,
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
that mention it are open access. We also incorporate text-mined links
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
<LI>Proteins from <A HREF="http://www.brenda-enzymes.org/index.php">BRENDA</A>, a curated database of enzymes, are included if they are linked to a paper in PubMed and their full sequence is known.
  <LI>Every protein from <A HREF="http://ecocyc.org">EcoCyc</A>, a curated
database of the proteins in <i> Escherichia coli</i> K-12, is included, regardless
of whether they are characterized or not.
<LI>Proteins from the <A HREF="http://metacyc.org">MetaCyc</A> metabolic pathway database are included if they are linked to a paper in PubMed and their full sequence is known.
<LI>Every protein from <A HREF="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3245046/">CharProtDB</A>, a database of experimentally characterized protein annotations, is included.
<LI>Proteins from the <A HREF="http://www.cazy.org/">CAZy</A> database of carbohydrate-active enzymes are included if they are associated with an Enzyme Classification number. Even though CAZy does not provide links from individual protein sequences to papers, these should all be experimentally-characterized proteins.
<LI>Proteins from the <A HREF="http://rebase.neb.com/rebase/rebase.html">REBASE</A> database of restriction enzymes are included if they have known specificity.
<LI>Every protein with an evidence-based reannotation (based on mutant phenotypes) in the <A HREF="http://fit.genomics.lbl.gov/">Fitness Browser</A> is included.
  </UL>Except for GeneRIF,
the curated entries include a short curated
description of the protein's function.  Many of these entries also
link to articles in <A
HREF="http://www.ncbi.nlm.nih.gov/pubmed/">PubMed</A>.

<P>For more information see the
<A title="PaperBLAST: Text Mining Papers for Information about Homologs" HREF="http://msystems.asm.org/content/2/4/e00039-17">PaperBLAST paper</A> (<i>mSystems</i> 2017)
or the <A HREF="https://github.com/morgannprice/PaperBLAST">code</A>.
You can download PaperBLAST's database <A HREF="https://github.com/morgannprice/PaperBLAST#download">here</A>.

<P>Changes since the paper was written:

<UL>
<LI>April 2019: EuropePMC now returns table entries in their search results. This has expanded PaperBLAST's database, but most of the new entries are of low relevance, and the resulting snippets are often just lists of locus tags with annotations.
<LI>February 2018: the alignment page reports the conservation of the hit's functional sites (if available from from Swiss-Prot or UniProt)
<LI>January 2018: incorporated BRENDA.
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

<P>Many important articles are missing from PaperBLAST, either because
the article's full text is not in EuropePMC (as for many older
articles), or because the paper does not mention a protein identifier such as a locus tag, or because of PaperBLAST's heuristics. If you notice an
article that characterizes a protein's function but is missing from
PaperBLAST, please notify the curators at <A
HREF="http://www.uniprot.org/update">UniProt</A>
or add an entry to <A
HREF="https://www.ncbi.nlm.nih.gov/gene/submit-generif">GeneRIF</A>.
Entries in either of these databases will eventually be incorporated
into PaperBLAST.  Note that to add an entry to UniProt, you will need
to find the UniProt identifier for the protein.  If the protein is not
already in UniProt, you can ask them to create an entry.  To add an
entry to GeneRIF, you will need an NCBI Gene identifier, but
unfortunately many prokaryotic proteins in RefSeq do not have
corresponding Gene identifers.

<H3>References</H3>

<P><A HREF="http://msystems.asm.org/content/2/4/e00039-17">PaperBLAST: Text-mining papers for information about homologs.</A><BR>
 <small>M. N. Price and A. P. Arkin (2017). mSystems, 10.1128/mSystems.00039-17.</small>

<P><A HREF="https://doi.org/10.1093/nar/gkx1005">Europe PMC in 2017.</A><BR>
<small>M. Levchenko et al (2017). Nucleic Acids Research, 10.1093/nar/gkx1005.</small>

<P><A HREF="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1480312/">Gene indexing: characterization and analysis of NLM's GeneRIFs.</A><BR>
<small>J. A. Mitchell et al (2003). AMIA Annu Symp Proc 2003:460-464.</small>

<P><A HREF="https://doi.org/10.1093/nar/gkw1099">UniProt: the universal protein knowledgebase.</A><BR>
<small>The UniProt Consortium (2016). Nucleic Acids Research, 10.1093/nar/gkw1099.</small>

<P><A HREF="https://doi.org/10.1093/nar/gkw952">BRENDA in 2017: new perspectives and new tools in BRENDA.</A><BR>
<small>S. Placzek et al (2017). Nucleic Acids Research, 10.1093/nar/gkw952.</small>

<P><A HREF="https://doi.org/10.1093/nar/gkw1003">The EcoCyc database: reflecting new knowledge about Escherichia coli K-12.</A><BR>
<small>I. M. Keeseler et al (2016). Nucleic Acids Research, 10.1093/nar/gkw1003.</small>

<P><A HREF="https://doi.org/10.1093/nar/gkx935">The MetaCyc database of metabolic pathways and enzymes.</A><BR>
<small>R. Caspi et al (2018). Nucleic Acids Research, 10.1093/nar/gkx935.</small>

<P><A HREF="https://academic.oup.com/nar/article/40/D1/D237/2903195">CharProtDB: a database of experimentally characterized protein annotations.</A><BR>
<small>R. Madupu et al (2012). Nucleic Acids Research, 10.1093/nar/gkr1133.</small>

<P><A HREF="https://doi.org/10.1093/nar/gkt1178">The carbohydrate-active enzymes database (CAZy) in 2013.</A><BR>
<small>V. Lombard et al (2014). Nucleic Acids Research, 10.1093/nar/gkt1178.</small>

<P><A HREF="https://doi.org/10.1093/nar/gku1046">REBASE - a database for DNA restriction and modification: enzymes, genes and genomes.</A><BR>
<small>R. J. Roberts et al (2015). Nucleic Acids Research, 10.1093/nar/gku1046.</small>

<P><A HREF="http://dx.doi.org/10.1101/072470">Deep annotation of protein function across diverse bacteria from mutant phenotypes.</A><BR>
<small>M. N. Price et al (2016). bioRxiv, 10.1101/072470.</small>

END
    ;

my $title = "PaperBLAST";
start_page('title' => $title);
print <<END
<SCRIPT src="../static/pb.js"></SCRIPT>
<SCRIPT src="http://fit.genomics.lbl.gov/d3js/d3.min.js"></SCRIPT>
<SCRIPT src="http://fit.genomics.lbl.gov/images/fitblast.js"></SCRIPT>
END
  ;

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
    my ($gene, $geneId);
    $gene = $dbh->selectrow_hashref("SELECT * from Gene WHERE geneId = ?", {}, $short);
    $geneId = $gene->{geneId} if $gene;
    if (! $gene) {
      $gene = $dbh->selectrow_hashref("SELECT * from CuratedGene WHERE db = 'SwissProt' AND protId = ?",
                                               {}, $short);
      $geneId = "SwissProt::".$short if $gene;
    }
    if (! $gene && $short =~ m/^(.*)::(.*)$/) {
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
        close(SEQ) || die "Error reading $seqFile";
        $query = ">$geneId $desc ($org)\n$seq\n";
      }
    }
  }

  # is it a fitness browser locus tag?
  if (!defined $query && $short =~ m/^[0-9a-zA-Z_]+$/ && -e $fbdata ) {
    my $fbdbh = DBI->connect("dbi:SQLite:dbname=$fbdata/feba.db","","",{ RaiseError => 1 }) || die $DBI::errstr;
    my $gene = $fbdbh->selectrow_hashref("SELECT * FROM Gene WHERE locusId = ? OR sysName = ? OR gene = ? LIMIT 1",
                                         {}, $short, $short, $short);
    if ($gene) {
      my $seqId = $gene->{orgId} . ":" . $gene->{locusId};
      die "Missing aaseqs file for fitness browser: $fbdata/aaseqs\n"
        unless -e "$fbdata/aaseqs";
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
        my $showId = $gene->{sysName} || $gene->{locusId};
        my $org = $fbdbh->selectrow_hashref("SELECT * FROM Organism WHERE orgId = ?",
                                            {}, $gene->{orgId})
          || die "No such organism in $fbdata/feba.db: $gene->{orgId}";
        my $orgdesc = join(" ", $org->{genus}, $org->{species}, $org->{strain});
        $query = ">$showId $gene->{gene} $gene->{desc} ($orgdesc)\n$seq\n";
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
  my $shortSafe = HTML::Entities::encode($short);
  &fail("Sorry -- we were not able to find a protein sequence for the identifier <b>$shortSafe</b>. We checked it against our database of proteins that are linked to papers, against UniProt (including their ID mapping service), against MicrobesOnline, and against the NCBI protein database (RefSeq and Genbank). Please use the sequence as a query instead.")
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

if (!defined $seq && ! $more_subjectId) {
    my $exampleId = "3615187";
    my $refseqId = "WP_012018426.1";
    print
        GetMotd(),
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
              -title => "Show PaperBLAST hits" }, i("ercA")) . "."),
        h3("Related Tools"),
        start_ul(),
        li(a({ -href => "genomeSearch.cgi",
               -title => "Search a genome for proteins that are related to a query term" },
           "Curated BLAST for Genomes")),
        li(a({ -href => 'gapView.cgi', -title => "Annotate the amino acid biosynthesis pathways in a genome" },
             "GapMind for amino acid biosynthesis")),
        li(a({ -href => "hmmSearch.cgi", -title => "Literature for all proteins that match an HMM" },
             "Family Search vs. Papers")),
        li(a({ -href => "http://papers.genomics.lbl.gov/vspdb",
               -title => "Papers vs. PDB: Well-studied proteins that lack structural information" },
            "Papers vs. PDB")),
        li("Papers vs. PFam",
           start_ul(),
           li(a({-href => "http://papers.genomics.lbl.gov/vspfam/dufs.html"}, "DUFs with characterized representatives")),
           li(a({-href => "http://papers.genomics.lbl.gov/vspfam/nopfam.html"}, "Characterized proteins not in PFam")),
           end_ul()),
        end_ul(),
        $documentation;
    finish_page();
} else {
    if ($more_subjectId) {
      print h3("Full List of Papers Linked to", HTML::Entities::encode($more_subjectId));
    } else {
      die "No sequence to search" unless $seq;
      my $initial = substr($seq, 0, 10);
      my $seqlen = length($seq);
      $initial .= "..." if $seqlen > 10;
      $initial = "$seqlen a.a., $initial" if $hasDef;
      print
        GetMotd(),
        h3("PaperBLAST Hits for", HTML::Entities::encode($def), "($initial)");

      my @nt = $seq =~ m/[ACGTUN]/g;
      my $fACGTUN = scalar(@nt) / $seqlen;
      if ($fACGTUN >= 0.9) {
        printf("<P><font color='red'>Warning: sequence is %.1f%% nucleotide characters -- are you sure this is a protein query?</font>",
               100 * $fACGTUN);
      }

      my $newline = "%0A";
      my $cdd_url = "http://www.ncbi.nlm.nih.gov/Structure/cdd/wrpsb.cgi?seqinput=>"
        . HTML::Entities::encode($def) . $newline . $seq;
      print qq{ <DIV style="float:right; padding-left: 10px; width:25%; background-color: #EEEEEE;">
                <P>Also see
                <A HREF="$cdd_url" TITLE="Compare your sequence to NCBI's Conserved Domains Database">Conserved Domains</A>
                or
                <A TITLE="Fitness BLAST compares your sequence to bacterial proteins that have mutant phenotypes"
                   NAME="#fitness">Fitness BLAST:</A>
                <SPAN ID="fitblast_short"><SMALL>loading...</SMALL></SPAN>
                <SCRIPT>
                var server_root = "http://fit.genomics.lbl.gov/";
                var seq = "$seq";
                fitblast_load_short("fitblast_short", server_root, seq);
                </SCRIPT>
                </DIV>
             };
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
      foreach my $row (@hits) {
        my ($queryId,$subjectId,$percIdentity,$alnLength,$mmCnt,$gapCnt,
            $queryStart,$queryEnd,$subjectStart,$subjectEnd,
            $eVal,$bitscore) = @$row;
        next if exists $seen_subject{$subjectId};
        $seen_subject{$subjectId} = 1;
        my @genes = &UniqToGenes($dbh, $subjectId);
        # Ok if not using the uniq id in the link -- litSearch.cgi can disambiguate
        my $coverage_html = "";
        $coverage_html = &simstring(length($seq), $genes[0]{protein_length},
                                    $queryStart, $queryEnd, $subjectStart, $subjectEnd,
                                    $percIdentity,$eVal,$bitscore,
                                    $def, $genes[0]{showName}, $seq, $subjectId)
          unless $more_subjectId;
        print GenesToHtml($dbh, $subjectId, \@genes, $coverage_html, $maxPapers);
        print "\n";
      }
    }

    if (! $more_subjectId) {
      my @pieces = $seq =~ /.{1,60}/g;
      print h3("Query Sequence"),
        p({-style => "font-family: monospace;"}, small(join(br(), ">" . HTML::Entities::encode($def), @pieces)));
    }
    print h3(a({-href => "litSearch.cgi"}, "New Search")),
      $documentation;
    finish_page();
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
    return a({ -title => $title,
               -href => "showAlign.cgi?def1=$def1&def2=$def2&seq1=$seq1&acc2=$acc2",
               -style => "font-family: sans-serif; font-size: smaller;" },
             "$percIdentity% identity, $percentCov% coverage");
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
  return $results if defined $results && $results =~ m/^>/;
  return undef;
}

sub UniProtToQuery($) {
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

sub commify($) {
    local $_  = shift;
    1 while s/^(-?\d+)(\d{3})/$1,$2/;
    return $_;
}


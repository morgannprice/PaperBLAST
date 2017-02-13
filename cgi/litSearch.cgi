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
# query (protein sequence in FASTA or UniProt or plain format)
# vimss -- specify the vimss id to search for (mostly to make testing easier)

use strict;
use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use Time::HiRes qw{gettimeofday};
use DBI;

sub fail($);
sub simstring($$$$$$$$$$$);

my $tmpDir = "../tmp";
my $blastall = "../bin/blast/blastall";
my $nCPU = 4;
my $base = "../data";
my $blastdb = "$base/litsearch.faa";
my $sqldb = "$base/litsearch.db";
die "No such executable: $blastall" unless -x $blastall;
die "No such file: $blastdb" unless -e $blastdb;
die "No such file: $sqldb" unless -e $sqldb;

my $cgi=CGI->new;
my $query = $cgi->param('query') || "";

my $dbh = DBI->connect("dbi:SQLite:dbname=$sqldb","","",{ RaiseError => 1 }) || die $DBI::errstr;
my $mo_dbh = DBI->connect('DBI:mysql:genomics:pub.microbesonline.org', "guest", "guest")
    || die $DBI::errstr;

my $documentation = <<END
<H3><A NAME="works">How It Works</A></H3>

<P>PaperBLAST relies on a database of protein sequences that are linked to scientific publications. These links come from automated text searches against <A HREF="http://europepmc.org/">EuropePMC</A> and from the manually curated information in <A HREF="http://www.uniprot.org/">Swiss-Prot</A> (the manually reviewed part of UniProt). Given this database and a query sequence, we use <A HREF="https://en.wikipedia.org/wiki/BLAST">protein-protein BLAST</A> to find similar sequences with E &lt; 0.001. As of January 2017, PaperBLAST contains over 200,000 proteins and links to over 50,000 papers.

<P>To link proteins in sequenced genomes to papers in EuropePMC, we query every locus tag or <A HREF="https://www.ncbi.nlm.nih.gov/refseq/">RefSeq</A> protein id that appears in the open access part of EuropePMC. We obtain the protein sequences and identifiers from <A HREF="http://www.microbesonline.org/">MicrobesOnline</A> or from RefSeq. We use queries of the form "locus_tag AND genus_name" to try to ensure that the paper is actually discussing the gene as opposed to something else whose identifier happens to match a locus tag. Note that EuropePMC indexes secret papers as well as open-access papers, so some of the links may be to papers that you cannot read (and that it would be illegal for our computers to read).</P>
<P>We also query EuropePMC for every locus tag that appears in the 300 most-referenced genomes. So, a gene may appear in the PaperBLAST results even though all of the papers that mention it are secret.</P>

<P>Finally, we index proteins from Swiss-Prot if their curators identified experimental evidence for the protein\'s function, as indicated by evidence code ECO:0000269. Most of these entries include links to articles in <A HREF="http://www.ncbi.nlm.nih.gov/pubmed/">PubMed</A><sup>&reg;</sup>. (We also use PubMed<sup>&reg;</sup> abstracts to help identify snippets of text that contain the locus tag from the secret papers.)

<P>The code for PaperBLAST is available <A HREF="https://github.com/morgannprice/PaperBLAST">here</A>.

<H3><A NAME="secret">Secrets</A></H3>

<P>PaperBLAST cannot provide snippets for most of the "secret" papers that are published in non-open-access journals. This limitation applies even if the paper is marked as "free" on the publisher\'s web site and is available in PubmedCentral or EuropePMC. Because these papers are not open access, PaperBLAST cannot download these papers en masse and scan them for snippets. If a journal that you publish in is marked as secret, please consider publishing elsewhere. Overall, PaperBLAST has access to the full text for about half of the papers.

<center><A HREF="http://morgannprice.org/">Morgan Price</A><BR>
<A HREF="http://genomics.lbl.gov/">Arkin group</A><BR>
Lawrence Berkeley National Laboratory<BR>
January 2017</center>
END
    ;

my $title = "PaperBLAST";
# utf-8 because that is the encoding used by EuropePMC
print
    header(-charset => 'utf-8'),
    start_html($title);

my $seq;
my $def = "";

if ($cgi->param('vimss')) {
    my $locusId = $cgi->param('vimss');
    die "Invalid vimss argument" unless $locusId =~ m/^\d+$/;
    my ($aaseq) = $mo_dbh->selectrow_array(
        qq{SELECT sequence FROM Locus JOIN AASeq USING (locusId,version)
            WHERE locusId = ? AND priority=1 },
        {}, $locusId);
    die "Unknown VIMSS id $locusId" unless defined $aaseq;
    $query = ">VIMSS$locusId\n$aaseq\n";
}
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

if (!defined $seq) {
    my $exampleId = "3615187";
    my $refseqId = "WP_012018426.1";
    print
        h2($title),
        p(b("Find papers about a protein or its homologs")),
        start_form( -name => 'input', -method => 'POST', -action => 'litSearch.cgi'),
        p("Enter a sequence in FASTA or Uniprot format: ",
          br(),
          textarea( -name  => 'query', -value => '', -cols  => 70, -rows  => 10 )),
        p(submit('Search'), reset()),
        end_form,
        p("Example: an 'alcohol dehydrogenase'",
          small("(" .
                a({ -href => "https://www.ncbi.nlm.nih.gov/protein/$refseqId",
                    -style => "color: black;" },
                  "$refseqId" )
                . ")"),
          "is actually the regulator",
          a({ -href => "litSearch.cgi?vimss=$exampleId",
              -title => "Show PaperBLAST hits" },
            i("ercA")) . "."),
        $documentation;
} else {
    my $procId = $$;
    my $timestamp = int (gettimeofday() * 1000);
    my $filename = $procId . $timestamp;
    my $seqFile = "$tmpDir/$filename.fasta";
    my $seqlen = length($seq);

    my $initial = substr($seq, 0, 10);
    $initial .= "..." if $seqlen > 10;
    $initial = "$seqlen a.a., $initial" if $hasDef;

    print
        h2("PaperBLAST Hits for $def ($initial)"),
        "\n";
    open(SEQ, ">", $seqFile) || die "Cannot write to $seqFile";
    print SEQ ">$def\n$seq\n";
    close(SEQ) || die "Error writing to $seqFile";
    my $hitsFile = "$tmpDir/$filename.hits";
    system($blastall, "-p", "blastp", "-d", $blastdb, "-i", $seqFile, "-o", $hitsFile,
           "-e", 0.001, "-m", 8, "-a", $nCPU, "-F", "m S") == 0 || die "Error running blastall: $!";
    my @hits = ();
    open(HITS, "<", $hitsFile) || die "Cannot read $hitsFile";
    while(<HITS>) {
        chomp;
        my @F = split /\t/, $_;
        push @hits, \@F;
    }
    close(HITS) || die "Error reading $hitsFile";
    unlink($seqFile);
    unlink($hitsFile);
    my $nHits = scalar(@hits);
    if ($nHits == 0) {
        print p("Sorry, no hits to proteins in the literature.");
    } else {
        print p("Found $nHits similar proteins in the literature:"), "\n";
        print $cgi->start_ul, "\n";
        foreach my $row (@hits) {
            my ($queryId,$subjectId,$percIdentity,$alnLength,$mmCnt,$gapCnt,$queryStart,$queryEnd,$subjectStart,$subjectEnd,$eVal,$bitscore) = @$row;
            my $gene = $dbh->selectrow_hashref("SELECT * FROM Gene WHERE geneId = ?", {}, $subjectId);
            if (defined $gene) {
                my $papers = $dbh->selectall_arrayref(qq{ SELECT DISTINCT * FROM GenePaper
                                                          LEFT JOIN PaperAccess USING (pmcId,pmId)
                                                          WHERE geneId = ? ORDER BY year DESC },
                                                      { Slice => {} }, $subjectId);
                my $n = scalar(@$papers);
                next unless $n > 0; # not sure this should be possible
                my @terms = map { $_->{queryTerm} } @$papers;
                my %terms = map { $_ => 1 } @terms;
                @terms = sort keys %terms;
                my $subjectShow = join(", ", @terms);
                my $desc = $gene->{desc};
                my $URL = "";

                if ($subjectId =~ m/^VIMSS(\d+)$/) {
                    my $locusId = $1;
                    $URL = "http://www.microbesonline.org/cgi-bin/fetchLocus.cgi?locus=$locusId";
                    if ($desc eq "") {
                        # Fetch the gene description from MicrobesOnline if not in the litsearch db
                        ($desc) = $mo_dbh->selectrow_array(
                            qq{SELECT description FROM Locus JOIN Description USING (locusId,version)
                           WHERE locusId = ? AND priority=1 },
                            {}, $locusId);
                    }
                } elsif ($subjectId =~ m/^[A-Z]+_[0-9]+[.]\d+$/) { # refseq
                    $URL = "http://www.ncbi.nlm.nih.gov/protein/$subjectId";
                } elsif ($subjectId =~ m/^[A-Z][A-Z0-9]+$/) { # SwissProt/TREMBL
                    $URL = "http://www.uniprot.org/uniprot/$subjectId";
                }
                die "Cannot build URL for subject $subjectId" if $URL eq "";
                $subjectShow = a({ -href => $URL,
                                   -title => $desc || "no description" },
                                 $subjectShow);
                print
                    li($subjectShow, "from", i($gene->{organism}),
                       &simstring(length($seq), $gene->{protein_length},
                                  $queryStart,$queryEnd,$subjectStart,$subjectEnd,$percIdentity,$eVal,$bitscore,
                                  $def, join(", ", @terms), $seq, $subjectId )),
                    $cgi->start_ul(),
                    "\n";
                foreach my $paper (@$papers) {
                    my $url = undef;
                    my $pubmed_url = "http://www.ncbi.nlm.nih.gov/pubmed/" . $paper->{pmId};
                    if ($paper->{pmcId}) {
                        $url = "http://www.ncbi.nlm.nih.gov/pmc/articles/" . $paper->{pmcId};
                    } elsif ($paper->{pmid}) {
                        $url = $pubmed_url;
                    } elsif ($paper->{doi}) {
                        $url = "http://doi.org/" . $paper->{doi};
                    }
                    my $title = $paper->{title};
                    $title = a({-href => $url}, $title) if defined $url;
                    my $authorShort = $paper->{authors};
                    $authorShort =~ s/ .*//;
                    my $extra = "";
                    $extra = small(a({-href => $pubmed_url}, "(PubMed)"))
                        if !$paper->{pmcId} && $paper->{pmId};
                    print li("$title,",
                             a({-title => $paper->{authors}}, "$authorShort,"),
                             $paper->{journal}, $paper->{year}, $extra);
                    my $snippets = [];
                    # SELECT DISTINCT because results can be repeated; i.e. some pubmed abstracts
                    # show up multiple times
                    $snippets = $dbh->selectall_arrayref(
                        "SELECT DISTINCT * from Snippet WHERE geneId = ? AND pmcId = ? AND pmId = ?",
                        { Slice => {} },
                        $subjectId, $paper->{pmcId}, $paper->{pmId})
                        if $paper->{pmcId} || $paper->{pmId};
                    if (@$snippets > 0) {
                        print $cgi->start_ul(), "\n";
                        foreach my $snippet (@$snippets) {
                            my $text = $snippet->{snippet};
                            my $term = $snippet->{queryTerm};
                            $text =~ s!$term!<span style="color: red;">$term</span>!g;
                            print li(small(qq{"...$text..."}));
                        }
                        print $cgi->end_ul(), "\n";
                    } else {
                        # Explain why there is no snippet
                        my $excuse;
                        if ($paper->{access} eq "full") {
                            $excuse = "This term was not found in the full text, sorry.";
                        } elsif ($paper->{isOpen} == 1) {
                            if ($paper->{access} eq "abstract") {
                                $excuse = "No snippet, sorry (this paper is open access but only the abstract was searched)";
                            } else {
                                $excuse = "No snippet, sorry (this paper is open access but neither the abstract nor the full text was searched)";
                            }
                        } elsif ($paper->{journal} eq "") {
                            $excuse = qq{This paper is <A HREF="#secret">secret</A>, sorry};
                        } else {
                            $excuse = qq{$paper->{journal} is <A HREF="#secret">secret</A>, sorry};
                        }
                        print $cgi->start_ul() . li(small($excuse)) . $cgi->end_ul()
                            if defined $excuse;
                    }
                }
                print $cgi->end_ul(), "\n";
            } elsif ($subjectId =~ m/^gnl\|ECOLI\|(.*)$/) { # EcoCyc gene
                my $protein_id = $1;
                my $URL = "http://ecocyc.org/gene?orgid=ECOLI&id=$protein_id";
                my $gene = $dbh->selectrow_hashref("SELECT * FROM EcoCyc WHERE protein_id = ?", {}, $protein_id);
                die "Unrecognized subject $protein_id" unless defined $gene;
                my $organism = "Escherichia coli K-12";
                my @ids = ( $gene->{protein_name}, $gene->{bnumber} );
                @ids = grep { $_ ne "" } @ids;
                my $show = join(" or ", @ids) || $protein_id;
                my $papers = $dbh->selectcol_arrayref(
                    "SELECT pmId FROM EcoCycToPubMed WHERE protein_id = ?",
                    {}, $protein_id);
                my $show_papers = "";
                my $pubmed_url = "http://www.ncbi.nlm.nih.gov/pubmed/" . join(",", @$papers);
                my $nPapers = scalar(@$papers);
                if ($nPapers > 0) {
                    $show_papers = " (see " . a({-href => $pubmed_url},
                                                $nPapers > 1 ? "papers" : "paper")
                        . ")";
                }
                print li(a({-href => $URL, -title => "see details in EcoCyc"}, "$show (EcoCyc)"),
                          " from <I>$organism</I>: ",
                         $gene->{desc}, $show_papers);

            } else { # UniProt gene
                $gene = $dbh->selectrow_hashref("SELECT * FROM UniProt WHERE acc = ?", {}, $subjectId);
                die "Unrecognized subject $subjectId" unless defined $gene;
                my @comments = split /_:::_/, $gene->{comment};
                @comments = map { s/[;. ]+$//; $_; } @comments;
                @comments = grep m/^FUNCTION|COFACTOR|CATALYTIC|ENZYME|DISRUPTION/, @comments;
                my $comment = "<LI>" . join("<LI>\n", @comments);
                my @pmIds = $comment =~ m!{ECO:0000269[|]PubMed:(\d+)}!;
                $comment =~ s!{ECO:[A-Za-z0-9|:,.| -]+}!!g;
                if (@pmIds > 0) {
                    my %seen = ();
                    @pmIds = grep { my $keep = !exists $seen{$_}; $seen{$_} = 1; $keep; } @pmIds;
                    my $note = @pmIds > 1 ? scalar(@pmIds) . " papers" : "paper";
                    $comment .= li("See ",
                                   a({-href => "http://www.ncbi.nlm.nih.gov/pubmed/" . join(",",@pmIds) },
                                     $note));
                }
                print li(a({-href => "http://www.uniprot.org/uniprot/$gene->{acc}"}, $gene->{acc}),
                         b($gene->{desc}),
                         "from",
                         i($gene->{organism}),
                         &simstring(length($seq), $gene->{protein_length},
                                    $queryStart,$queryEnd,$subjectStart,$subjectEnd,$percIdentity,$eVal,$bitscore,
                                    $def, $subjectId, $seq, $subjectId ),
                         $cgi->start_ul,
                         small($comment),
                         $cgi->end_ul );
            }
        }
        print $cgi->end_ul, "\n";
    }

    print qq{<script src="http://fit.genomics.lbl.gov/d3js/d3.min.js"></script>
             <script src="http://fit.genomics.lbl.gov/images/fitblast.js"></script>
             <H3><A title="Fitness BLAST searches for similarity to bacterial proteins that have mutant phenotypes" HREF="http://fit.genomics.lbl.gov/" NAME="#fitness">Fitness Blast</A></H3>
             <P><DIV ID="fitblast_short"></DIV></P>
             <script>
             var server_root = "http://fit.genomics.lbl.gov/";
             var seq = "$seq";
             fitblast_load_short("fitblast_short", server_root, seq);
             </script>
    };

    my @pieces = $seq =~ /.{1,60}/g;
    print
        h3("Query Sequence"),
        p({-style => "font-family: monospace;"}, small(join(br(), ">$def", @pieces))),
        h3(a({-href => "litSearch.cgi"}, "New Search")),
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

sub simstring($$$$$$$$$) {
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

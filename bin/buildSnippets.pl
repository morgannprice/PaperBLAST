#!/usr/bin/perl -w
use strict;
use XML::LibXML;
use Getopt::Long;

my $snippetBefore = 15;
my $snippetAfter = 30;

my $usage = <<END
buildSnippets.pl -in xml -out snippets oa.papers ... refseq.papers
	The -in file should be a file in xml format (or gzipped, if its name
	ends in .gz).

	The papers files are from parseEuropePMCHits.pl and include
	queryId, queryTerm, pmcId.

	Output will be tab-delimited with the fields queryId,
	queryTerm, pmcId, snippet

	Optional arguments control the size of each snippet:
	-before $snippetBefore
        -after $snippetAfter

	There is a limit of 2 snippets per gene per paper
END
    ;

sub ProcessArticle($$);
sub NodeText($);

{
    my ($infile, $outfile);
    my @paperfiles;
    die $usage
        unless GetOptions( 'in=s' => \$infile,
                           'out=s' => \$outfile,
                           'before=i' => \$snippetBefore,
                           'after=i' => \$snippetAfter)
        && defined $infile
        && defined $outfile
        && scalar(@ARGV) > 0;
    @paperfiles = @ARGV;
    die "-before cannot be negative" if $snippetBefore < 0;
    die "-after must be positive" if $snippetAfter < 1;

    die "No such file: $infile\n" unless -e $infile;
    my $gunzip = $infile =~ m/[.]gz$/;
    if ($gunzip) {
        open(IN, "zcat $infile |") || die "Cannot run zcat on $infile";
    } else {
        open(IN, "<", $infile) || die "Cannot read $infile";
    }

    my %papers = (); # pmcId => queryTerm => list of queryId
    foreach my $paperfile (@paperfiles) {
        die "No such file: $paperfile\n" unless -e $paperfile;
    }
    foreach my $paperfile (@paperfiles) {
        open(PAPERS, "<", $paperfile) || die "Error reading $paperfile";
        while(my $line = <PAPERS>) {
            chomp $line;
            my ($queryId, $queryTerm, $pmcId) = split /\t/, $line;
            die "Not enough fields in line\n$line\nin $paperfile" unless defined $pmcId;
            next if $pmcId eq "";
            push @{ $papers{$pmcId}{$queryTerm} }, $queryId;
        }
        close(PAPERS) || die "Error reading $paperfile";
    }
    print STDERR "Read " . scalar(keys %papers) . " pmcIds to search\n";

    my $firstline = <IN>;
    chomp $firstline;
    die "First line should be: <articles>\n" unless $firstline eq "<articles>";
    open(OUT, ">", $outfile) || die "Error writing to $outfile";
    my $n = 0;
    my $nSeen = 0;
    while(my $line = <IN>) {
        chomp $line;
        next if $line eq "";
        if ($line eq "</articles>") {
            next;
            last;
        }
        $line =~ m/^<article[> ]/ || die "Does not begin with an article: " . substr($line, 0, 100);
        my @lines = ( $line );
        while ($lines[-1] !~ m!</article>!) {
            $line = <IN>;
            if (defined $line) {
                chomp $line;
                push @lines, $line;
            } else {
                last;
            }
        }
        if ($lines[-1] !~ m!</article>!) {
            print STDERR "Last article is truncated? Starts with " . substr($lines[0], 0, 100);
        }
        # without recover, get errors like
        # :1: parser error : Premature end of data in tag p line 1
        # l 2004 news article &#x0201c;Reaching across the Border with the SBRP&#x0201d; [

        my $dom = XML::LibXML->load_xml(string => join("\n",@lines), recover => 1);
        $n++;
        print STDERR "Parsed $n articles\n" if ($n % 1000) == 0;
        $nSeen += &ProcessArticle($dom, \%papers);
    }
    print STDERR "Processed $n articles with $nSeen relevant pmcIds\n";
}

sub ProcessArticle($$) {
    my ($dom,$papers) = @_;

    # Find the pmcid, i.e. <article-id pub-id-type="pmcid">1240566</article-id>
    # First, use an XPath query to find article-id objects with pub-id-type attributes,
    my @ids = $dom->findnodes(qq{//article-id[attribute::pub-id-type]});
    my $pmcid = "";
    # then get the pmcid values, if any, from each
    foreach my $id (@ids) {
        foreach my $attr ($id->attributes) {
            if ($attr->nodeName eq "pub-id-type" && $attr->getValue eq "pmcid") {
                $pmcid = $id->textContent;
            }
        }
    }
    if (@ids == 0) {
        print STDERR "Warning: no article-id elements\n";
    } elsif ($pmcid eq "") {
        print STDERR "Warning: no pmcid\n";
    }
    return 0 if $pmcid eq "";
    $pmcid = "PMC" . $pmcid;
    return 0 if !exists $papers->{$pmcid};

    my $text = NodeText($dom);
    $text =~ s/\P{IsASCII}//g; # remove wide characters
    my @words = split /\s+/, $text;
    my @subwords = @words;
    # removing leading and trailing punctuation from each word
    @subwords = map { s/[.,;()-]+$//; s/^[.,;()-]+//; uc($_); } @subwords;

    while (my ($queryTerm, $queryIds) = each %{ $papers->{$pmcid} }) {
        my @hits = grep { $subwords[$_] eq uc($queryTerm) } (0..(scalar(@words)-1));
        my @snippets = ();
        for (my $i = 0; $i < 2 && $i < scalar(@hits); $i++) {
            my $j = $hits[$i]; # which word matches
            my $j1 = $j - $snippetBefore;
            my $j2 = $j + $snippetAfter;
            $j1 = 0 if $j1 < 0;
            $j2 = scalar(@words)-1 if $j2 >= scalar(@words);
            push @snippets, join(" ", @words[$j1..$j2]);
        }

        my %queryIds = map { $_ => 1 } @$queryIds;
        my @queryIds = keys(%queryIds);
        foreach my $queryId (@queryIds) {
            foreach my $snippet (@snippets) {
                print OUT join("\t", $pmcid, $queryTerm, $queryId, $snippet)."\n";
            }
        }
        print STDERR "Warning: no snippets for $queryTerm in $pmcid\n" if scalar(@snippets) == 0;
    }
    return 1;
}

sub NodeText($) {
    my ($begin_node) = @_;
    my @work = ( $begin_node );
    # Depth first search -- work remembers the work left to do
    my @pieces = ();
    while(@work > 0) {
        my $node = shift @work;
        if ($node->nodeType == XML_TEXT_NODE) {
            push @pieces, $node->data;
        }
        unshift @work, $node->childNodes();
    }
    return join(" ", grep { $_ ne "" } @pieces);
}

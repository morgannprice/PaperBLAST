#!/usr/bin/perl -w
# Given an XML articles file from PMC Europe,
# find all the locus tag like words
use strict;
use XML::LibXML;

sub ProcessArticle($);
sub NodeText($);

{
    die "Run as a filter\n" unless @ARGV==0;
    my $firstline = <STDIN>;
    chomp $firstline;
    die "First line should be: <articles>\n" unless $firstline eq "<articles>";
    my $n = 0;
    while(my $line = <STDIN>) {
        chomp $line;
        next if $line eq "";
        if ($line eq "</articles>") {
            next;
            last;
        }
        $line =~ m/^<article[> ]/ || die "Does not begin with an article: " . substr($line, 0, 100);
        my @lines = ( $line );
        while ($lines[-1] !~ m!</article>!) {
            $line = <STDIN>;
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
        ProcessArticle($dom);
    }
    print STDERR "Processed $n articles\n";
}

sub ProcessArticle($) {
    my ($dom) = @_;
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

    # For starters, just print the text
    my $text = NodeText($dom);
    $text =~ s/\P{IsASCII}//g; # remove wide characters
    my @words = split /\s+/, $text;
    my %seen = ();
    foreach my $word (@words) {
        # Looking for potential locus tags with a letter, ending with at least 3 numbers in a row,
        # and possibly trailing punctuation
        $word =~ s/[.,;-]+$//;
        next unless $word =~ m/^[a-zA-Z][a-zA-Z90-9_]+\d\d\d$/;
        print "$pmcid\t$word\n" unless exists $seen{$word};
        $seen{$word} = 1;
    }
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

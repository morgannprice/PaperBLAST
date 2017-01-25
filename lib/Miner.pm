package Miner;
require Exporter;
use strict;
use XML::LibXML;
use LWP::Simple;

our (@ISA,@EXPORT);
@ISA = qw(Exporter);
@EXPORT = qw( RemoveWideCharacters TextSnippets NodeText XMLFileToText ElsevierFetch );

sub RemoveWideCharacters($) {
    my ($text) = @_;
    $text =~ s/\P{IsASCII}//g;
    return($text);
}

# text (with wide characters removed) to short list of snippets
sub TextSnippets($$$$$) {
    my ($text, $queryTerm, $snippetBeforeWords, $snippetAfterWords, $maxCharacters) = @_;
    my @words = split /\s+/, $text;
    my @subwords = @words;
    # removing leading and trailing punctuation from each word
    @subwords = map { s/[.,;()-]+$//; s/^[.,;()-]+//; uc($_); } @subwords;

    my @hits = grep { $subwords[$_] eq uc($queryTerm) } (0..(scalar(@words)-1));
    my @snippets = ();
    for (my $i = 0; $i < 2 && $i < scalar(@hits); $i++) {
        my $j = $hits[$i]; # which word matches
        my $j1 = $j - $snippetBeforeWords;
        my $j2 = $j + $snippetAfterWords;
        $j1 = 0 if $j1 < 0;
        $j2 = scalar(@words)-1 if $j2 >= scalar(@words);
            my $snippet = join(" ", @words[$j1..$j2]);
            while (length($snippet) > $maxCharacters) {
                $j1++;
                $j2--;
                if ($j2 <= $j1) {
                    last;
                } else {
                    $snippet = join(" ", @words[$j1..$j2]);
                }
            }
            push @snippets, $snippet if $snippet ne "";
    }
    return @snippets;
}

# XML to a string
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

sub ElsevierFetch($$$$) {
    my ($pubmedId, $key, $cachefile, $jour)= @_;
    my $results = LWP::Simple::get("http://api.elsevier.com/content/article/pubmed_id/$pubmedId?APIKey=$key&httpAccept=text/xml");
    $results = "" if !defined $results;
    if ($results =~ m/<coredata>/) {
        # success
        open(CACHE, ">:encoding(UTF-8)", "$cachefile.tmp") || die "Error writing to $cachefile.tmp";
        print CACHE $results;
        close(CACHE) || die "Error writing to $cachefile.tmp";
        rename("$cachefile.tmp",$cachefile) || die "Error renaming $cachefile.tmp to $cachefile";
        return 1;
    }
    # else
    if ($results eq "") {
        print STDERR "Failed: empty results from $jour for $pubmedId\n";
    } else {
        $results =~ s/\n//g;
        print STDERR "Failed: $results from $jour for $pubmedId\n";
    }
    return 0;
    
}

sub XMLFileToText($) {
    my ($file) = @_;
    return undef if ! -e $file;
    open(CACHE, "<:encoding(UTF-8)", $file) || die "Error reading cache file $file";
    my @lines = <CACHE>;
    close(CACHE) || die "Error reading cache file $file";
    my $dom = XML::LibXML->load_xml(string => join("\n",@lines), recover => 1);
    return undef if ! $dom;
    return &RemoveWideCharacters( &NodeText($dom) );
}

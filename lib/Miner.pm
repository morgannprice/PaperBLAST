package Miner;
require Exporter;
use strict;
use XML::LibXML;
use LWP::Simple;
use LWP::UserAgent;
use HTTP::Request;
use JSON;


our (@ISA,@EXPORT);
@ISA = qw(Exporter);
@EXPORT = qw( RemoveWideCharacters TextSnippets NodeText
              XMLFileToText TextFileToText PDFFileToText
              ElsevierFetch CrossrefToURL CrossrefFetch );

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

# Note use of UTF-8 encoding
sub SaveResultTo($$) {
    my ($results, $file) = @_;
    open(CACHE, ">:encoding(UTF-8)", "$file.tmp") || die "Error writing to $file.tmp";
    print CACHE $results;
    close(CACHE) || die "Error writing to $file.tmp";
    rename("$file.tmp",$file) || die "Error renaming $file.tmp to $file";
}

sub ElsevierFetch($$$$) {
    my ($pubmedId, $key, $cachefile, $jour)= @_;
    my $results = LWP::Simple::get("http://api.elsevier.com/content/article/pubmed_id/$pubmedId?APIKey=$key&httpAccept=text/xml");
    $results = "" if !defined $results;
    if ($results =~ m/<coredata>/) {
        # success
        &SaveResultTo($results, $cachefile);
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

# Returns empty string or file type (xml or pdf) and URL
sub CrossrefToURL {
    my ($doi, $jour) = @_;
    my $results = LWP::Simple::get("http://api.crossref.org/works/$doi");
    if (! $results) {
        print STDERR "Failed: no results from CrossRef for $jour, $doi\n";
        return "";
    }
    my $json = from_json($results);
    if (!defined $json) {
        print STDERR "Invalid JSON from CrossRef for $jour, $doi\n";
        return "";
    }
    my $msg = $json->{message};
    unless ($json->{status} eq "ok" && defined $msg) {
        print STDERR "Failed: status from CrossRef not ok for $jour, $doi\n";
        return "";
    }
    if (!exists $msg->{link}) {
        print STDERR "No link in CrossRef for $jour, $doi\n";
        return "";
    }
    my $xml = undef;
    my $pdf = undef;
    my $other = undef;
    foreach my $link (@{ $msg->{link} }) {
        next unless $link->{URL};
        if ($link->{"content-type"} eq "application/xml") {
            $xml = $link->{URL};
        } elsif ($link->{"content-type"} eq "application/pdf") {
            $pdf = $link->{URL};
        } else {
            $other = $link->{URL};
        }
    }
    return ("xml", $xml) if defined $xml;
    return ("pdf", $pdf) if defined $pdf;
    return ("other", $other) if defined $other;
    print STDERR "No link in CrossRef for $jour, $doi\n";
    return "";
}

# returns empty on failure and returns the reported Content-Type on success
sub CrossrefFetch($$$$) {
    my ($URL, $type, $file, $jour) = @_;
    print STDERR "Fetching $jour $URL $file\n"; # XXX
    my $request = HTTP::Request->new("GET" => $URL);
    $request->header("CR-Clickthrough-Client-Token" => "e5a0777f-ac991935-bc17d839-edb4695c");
    my $ua = LWP::UserAgent->new;
    my $res = $ua->request($request);
    if (! $res->is_success) {
        print STDERR "Failure: " . $res->status_line . "\n";
        return 0;
    }
    # If it is actually a PDF, then do not use UTF-8
    my $mode = $res->header('Content-Type') eq "application/pdf" ?
        ">" : ">:encoding(UTF-8)";
    open(FILE, $mode, "$file.tmp") || die "Error writing to $file.tmp";
    print FILE $res->content();
    close(FILE) || die "Error writing to $file.tmp";
    rename("$file.tmp",$file) || die "Error renaming $file.tmp to $file";
    return $res->header('Content-Type');
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

sub TextFileToText($) {
    my ($file) = @_;
    open(FILE, "<", $file) || die "Cannot read $file";
    my @lines = <FILE>;
    close(FILE) || die "Error reading $file";
    return join("", @lines);
}

sub PDFFileToText($) {
    my ($file) = @_;
    my $txtfile = "$file.$$.txt";
    unless (system("pdftotext $file $txtfile >& /dev/null") == 0) {
        print STDERR "Error: cannot run pdftotext on $file: $!";
        unlink($txtfile);
        return "";
    }
    my $lines  = &RemoveWideCharacters( &TextFileToText($txtfile) );
    unlink($txtfile);
    return $lines;
}

1;

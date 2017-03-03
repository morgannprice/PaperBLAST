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
              ElsevierFetch CrossrefToURL CrossrefFetch
              DomToPMCId FetchPubMed );

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
    # remove trailing text like .1 or .2 if this is a refseq protein (in case the text includes the version)
    @subwords = map { s/[.]\d$//; $_; } @subwords
      if ($queryTerm =~ m/^[A-Z]P_\d+$/);

    my @hits = grep { $subwords[$_] eq uc($queryTerm) } (0..(scalar(@words)-1));
    my @snippets = ();
    my $end1 = undef; # last word in 1st snippet
    for (my $i = 0; $i < scalar(@hits) && scalar(@snippets) < 2; $i++) {
        my $j = $hits[$i]; # which word matches
        my $j1 = $j - $snippetBeforeWords;
        my $j2 = $j + $snippetAfterWords;
        $j1 = 0 if $j1 < 0;
        $j2 = scalar(@words)-1 if $j2 >= scalar(@words);
        my $snippet = join(" ", @words[$j1..$j2]);
        my $nWordsTrim = 0;
        my $fBefore = $snippetBeforeWords/($snippetBeforeWords + $snippetAfterWords);
        # final boundaries after trimming
        my ($j1f, $j2f) = ($j1, $j2);
        while (length($snippet) > $maxCharacters) {
            $nWordsTrim++;
            my $trimLeft = int(0.5 + $nWordsTrim * $fBefore);
            my $trimRight = $nWordsTrim - $trimLeft;
            $j1f = $j1 + $trimLeft;
            $j2f = $j2 - $trimRight;
            if ($j2f < $j || $j1f > $j || $j2f <= $j1f) {
                last;
            } else {
                $snippet = join(" ", @words[$j1f..$j2f]);
            }
        }
        # Ignore overlapping snippets
        next if defined $end1 && $j1f <= $end1;
        if ($snippet ne "") {
            $end1 = $j2f if !defined $end1;
            push @snippets, $snippet;
        }
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
    # use ASCII7 to convert ligatures to standard text
    unless (system("pdftotext -enc ASCII7 $file $txtfile >& /dev/null") == 0) {
        print STDERR "Error: cannot run pdftotext on $file: $!";
        unlink($txtfile);
        return "";
    }
    my $lines  = &RemoveWideCharacters( &TextFileToText($txtfile) );
    unlink($txtfile);
    return $lines;
}

# Find the pmcid, i.e. <article-id pub-id-type="pmcid">1240566</article-id>
# or <article-id pub-id-type="pmc">2268633</article-id>
sub DomToPMCId($) {
    my ($dom) = @_;

    # Find the pmcid, i.e. <article-id pub-id-type="pmcid">1240566</article-id>
    # First, use an XPath query to find article-id objects with pub-id-type attributes,
    my @ids = $dom->findnodes(qq{//article-id[attribute::pub-id-type]});
    if (@ids == 0) {
        @ids = $dom->findnodes(qq{/article/front/article-meta/article-id[attribute::pub-id-type]});
    }
    if (@ids == 0) {
        print STDERR "Warning: no article-id elements in dom\n";
        return undef;
    }

    # then get the pmcid value, if any
    foreach my $id (@ids) {
        foreach my $attr ($id->attributes) {
            if ($attr->nodeName eq "pub-id-type"
                && ($attr->getValue eq "pmcid" || $attr->getValue eq "pmc")) {
                return $id->textContent;
            }
        }
    }

    # or, return the PPR id for preprints
    foreach my $id (@ids) {
        foreach my $attr ($id->attributes) {
            if ($attr->nodeName eq "pub-id-type"
                && $attr->getValue eq "archive") {
              return $id->textContent if $id->textContent =~ m/^PPR/;
            }
        }
    }

    return undef;
}

# returns a list of XML objects, one per pm id
sub FetchPubMed {
  my @pmids = @_;
  return () if @pmids == 0;
  my %known = ();
  @pmids = grep { my $keep = !exists $known{$_}; $known{$_} = 1; $keep } @pmids;
  my @out = ();
  while (@pmids > 0) {
    my $maxTries = 2;
    my $size = 50;
    $size = @pmids if $size > @pmids;
    my @set = ();
    foreach my $i (1..$size) {
      push @set, (shift @pmids);
    }
    my $iTry;
    for ($iTry = 0; $iTry < $maxTries; $iTry++) {
      sleep(1) if $iTry > 0; # wait for NCBI to recover
      my $URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.cgi?db=PubMed&rettype=xml&id="
        . join(",",@set);
      my $results = LWP::Simple::get($URL);
      unless (defined $results && $results =~ m/PubmedArticleSet/) {
        print STDERR "Warning: Could not fetch results for ids $set[0] to $set[-1], retrying\n";
        next;
      }
      my $dom = XML::LibXML->load_xml(string => $results, recover => 1);
      if (! $dom) {
        print STDERR "Warning: Could not parse results for ids $set[0] to $set[-1], retrying\n";
        next;
      }
      next unless $dom;
      my @children = $dom->findnodes("//PubmedArticle");
      if (scalar(@children) == 0) {
        print STDERR "Warning: no PubmedArticle nodes for ids $set[0] to $set[-1], retrying\n";
        next;
      }
      print STDERR "Warning: expected " . scalar(@set) . " entries for ids $set[0] to $set[-1] but got " . scalar(@children) . "\n"
        if scalar(@children) < scalar(@set);
      push @out, @children;
      last;
    }
    print STDERR "Warning: failed to fetch pubmd ids $set[0] to $set[-1]\n"
      if $iTry == $maxTries;
  }
  return @out;
}

1;

#!/usr/bin/perl -w
use strict;
use XML::LibXML;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use Miner;

my $snippetBefore = 15;
my $snippetAfter = 30;
my $maxChar = 1000;

my $usage = <<END
buildSnippets.pl -in xml -out snippets oa.papers ... refseq.papers
  or
buildSnippets.pl -list xml_file_list -out snippets oa.papers ... refseq.papers

	The -in file should be a file in xml format (or gzipped, if its name
	ends in .gz).

	The -list file should have a list of xml files (or xml.gz files), one per line.

	The papers files are from parseEuropePMCHits.pl and include
	queryId, queryTerm, pmcId. If queryTerm is quoted, the quotes are removed.

	Output will be tab-delimited with the fields
	pmcId, pmId, queryTerm, queryId, snippet

	Also writes to output_file.access, with tab-delimited records of the form
	pmcId, pmId, "full"
	if it was able to read the paper

	Optional arguments control the size of each snippet in words or characters:
	-before $snippetBefore
        -after $snippetAfter
        -maxChar $maxChar

	There is a limit of 2 snippets per gene per paper
END
    ;

my %papers = (); # pmcId => queryTerm => list of queryId
my %pmc2pm = (); # pmcId => pmId
my %doi2pmc = (); # doi => pmcId
my %pm2pmc = (); # pm => pmcId

sub ProcessArticle($);

{
    my ($infile, $listfile, $outfile);
    my @paperfiles;
    die $usage
        unless GetOptions( 'in=s' => \$infile,
                           'list=s' => \$listfile,
                           'out=s' => \$outfile,
                           'before=i' => \$snippetBefore,
                           'after=i' => \$snippetAfter,
                           'maxChar=i' => \$maxChar )
        && (defined $infile xor defined $listfile)
        && defined $outfile
        && scalar(@ARGV) > 0;
    @paperfiles = @ARGV;
    die "-before cannot be negative" if $snippetBefore < 0;
    die "-after must be positive" if $snippetAfter < 1;
    die "-maxChar must be at least 50" if $maxChar < 50;

    my $gunzip;
    my @files = ();
    my $in;
    if (defined $infile) {
      die "No such file: $infile\n" unless -e $infile;
      @files = ( $infile );
    } else {
      die "No such file: $listfile" unless -e $listfile;
      open(LIST, "<", $listfile) || die "Error reading $listfile";
      while(my $line = <LIST>) {
        chomp $line;
        die "No such file: $line" unless -e $line;
        push @files, $line;
      }
      close(LIST) || die "Error reading $listfile";
      print STDERR "Read list of " . scalar(@files) . " input files\n";
    }

    foreach my $paperfile (@paperfiles) {
        die "No such file: $paperfile\n" unless -e $paperfile;
    }
    foreach my $paperfile (@paperfiles) {
        open(PAPERS, "<", $paperfile) || die "Error reading $paperfile";
        while(my $line = <PAPERS>) {
            chomp $line;
            my ($queryId, $queryTerm, $pmcId, $pmId, $doi) = split /\t/, $line;
            die "Not enough fields in line\n$line\nin $paperfile" unless defined $pmId;
            next if $pmcId eq "";
            $queryTerm =~ s/^"//;
            $queryTerm =~ s/"$//;
            push @{ $papers{$pmcId}{$queryTerm} }, $queryId;
            $pmc2pm{$pmcId} = $pmId;
            $doi2pmc{$doi} = $pmcId if $doi ne "";
            $pm2pmc{$pmId} = $pmcId if $pmId ne "";
        }
        close(PAPERS) || die "Error reading $paperfile";
    }
    print STDERR "Read " . scalar(keys %papers) . " pmcIds to search\n";

    my $accfile = "$outfile.access";
    open(OUT, ">", $outfile) || die "Error writing to $outfile";
    open(ACC, ">", $accfile) || die "Error writing to $accfile";
    my $n = 0;
    my $nSeen = 0;
    foreach my $file (@files) {
      $gunzip = $file =~ m/[.]gz$/;
      my $fhIn;
      if ($gunzip) {
        open($fhIn, "zcat $file |") || die "Cannot run zcat on $file";
      } else {
        open($fhIn, "<", $file) || die "Cannot read $file";
      }
      my $state = {};
      while (my $dom = ReadXMLArticle($fhIn, $state)) {
        $n++;
        print STDERR "Parsed $n articles\n" if ($n % 1000) == 0;
        $nSeen += &ProcessArticle($dom);
      }
      close($fhIn) || die "Error reading $file" . ($gunzip ? " . with zcat" : "")
    }
    print STDERR "Processed $n articles with $nSeen relevant pmcIds\n";
    close(OUT) || die "Error writing to $outfile";
    close(ACC) || die "Error writing to $accfile";
}

sub ProcessArticle($) {
    my ($dom) = @_;

    my @ids = &DomToIds($dom);
    my $pmcId = IdsToPMCId(@ids);
    my $pmId = IdsToPubmedId(@ids);
    my $doi = IdsToDOI(@ids);
    $pmcId = $pm2pmc{$pmId} if !defined $pmcId && defined $pmId && exists $pm2pmc{$pmId};
    $pmcId = $doi2pmc{$doi} if !defined $pmcId && defined $doi && exists $doi2pmc{$doi};
    return 0 if !defined $pmcId;
    $pmcId = "PMC" . $pmcId if $pmcId =~ m/^\d/; # do not use this prefix for PPR ids

    my $pmcLookup = $pmcId; # it may show up in a different form in the search results
    if (!exists $papers{$pmcLookup} && $pmcId =~ m/^PPR/) {
      $pmcLookup =~ s/^PPR/PMC-/;
    }
    return 0 if !exists $papers{$pmcLookup};

    my $text = &RemoveWideCharacters( &NodeText($dom) );
    print ACC join("\t", $pmcId, $pmc2pm{$pmcId} || "", "full")."\n";

    while (my ($queryTerm, $queryIds) = each %{ $papers{$pmcLookup} }) {
        my @snippets = TextSnippets($text, $queryTerm, $snippetBefore, $snippetAfter, $maxChar);
        my %queryIds = map { $_ => 1 } @$queryIds;
        my @queryIds = keys(%queryIds);
        foreach my $queryId (@queryIds) {
            foreach my $snippet (@snippets) {
                print OUT join("\t", $pmcId, $pmc2pm{$pmcId} || "", $queryTerm, $queryId, $snippet)."\n";
            }
        }
        print STDERR "Warning: no snippets for $queryTerm in $pmcId\n" if scalar(@snippets) == 0;
    }
    return 1;
}

#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use JSON;

my $usage = <<END
parseEuropePMCHits.pl -in queryprot -hits epmc.hits -out out

queryprot must contain the organism, the query_term (or locus_tag),
the queryId, the sequence, and (optionally) the description

epmc.hits must contain the queryId and the hits in JSON format.

Writes:
out.queries with queryId, organism, protein_length, description
out.papers with queryId, query_term, pmcId, pmId, doi, title, authors, journal, year, isOpen
out.faa in FASTA format, with queryId as the item name
END
    ;

{
    my ($infile, $hitsfile, $outpre) = undef;
    die $usage
        unless GetOptions('in=s' => \$infile,
                          'hits=s' => \$hitsfile,
                          'out=s' => \$outpre)
        && @ARGV == 0
        && defined $infile && defined $hitsfile && defined $outpre;
    die "No such file: $infile\n" unless -e $infile;
    die "No such file: $hitsfile\n" unless -e $hitsfile;

    my %seen = (); # queryId => query => 1 if it is in the hits file already
    my %hasHit = (); # queryId => 1 if it has a paper about it
    open(HITS, "<:encoding(UTF-8)", $hitsfile) || die "Cannot read $hitsfile";
    open(PAPERS, ">:encoding(UTF-8)", "$outpre.papers") || die "Cannot write to $outpre.papers";
    while(my $line = <HITS>) {
        chomp $line;
        my ($queryId, $jsonstring) = split /\t/, $line;
        die "Invalid hits line\n$line\n" if !defined $jsonstring || $queryId eq "";
        if ($jsonstring eq "") {
            print STDERR "No hits entry for $queryId\n";
        }
        my $json;
        eval{ $json = from_json($jsonstring) };
        if (!defined $json) {
            print STDERR "Skipping bad hits line for $queryId\n";
            next;
        }
        # Sometimes this information is in "queryString" instead of in "query"
        my $queryText = $json->{request}{query} || $json->{request}{queryString};
        die "No query text for $queryId" unless defined $queryText;
        if (exists $seen{$queryId}{$queryText}) {
            print STDERR "Skipping duplicate hits line for $queryId $queryText\n";
            next;
        }
        $seen{$queryId}{$queryText} = 1;
        if ($json->{hitCount} > 0) {
            $hasHit{$queryId} = 1;
            $queryText =~ m/^(.*) AND/ || die "Invalid query term $queryText";
            my $queryTerm = $1;
            foreach my $result (@{ $json->{resultList}{result} }) {
                my $pmcId = $result->{pmcid} || "";
                my $id = $result->{id} || "";
                # I.e. use the PPR identifier not PMC-7028
                $pmcId = $id if $pmcId =~ m/-/ && $id;
                my $pmId = $result->{pmid} || "";
                my $doi = $result->{doi} || "";
                my $title = $result->{title} || "";
                $title =~ s/[.]$//;
                my $authors = $result->{authorString} || "";
                my $journal = $result->{journalTitle} || $result->{journalInfo}{journal}{title} || "";
                my $year = $result->{pubYear} || "";
                # Occasional tab characters at beginning or end of these fields, or \n at end, so
                my @out = ($queryId, $queryTerm, $pmcId, $pmId, $doi,
                                     $title, $authors, $journal, $year,
                                     $result->{isOpenAccess} =~ m/y/i ? 1 : 0);
                foreach my $field (@out) {
                  $field =~ s/[\r\n\t]/ /g;
                  $field =~ s/^ +//;
                  $field =~ s/ +$//;
                }
                print PAPERS join("\t", @out)."\n";
            }
        }
    }
    close(HITS) || die "Error reading $hitsfile";
    close(PAPERS) || die "Error writing to $outpre.papers";
    print STDERR "Parsed hits for " . scalar(keys %seen) . " queries, "
        . scalar(keys %hasHit) . " with hits\n";
    print STDERR "Wrote $outpre.papers\n";

    my %hasSeq = (); # queryId => 1 if has sequence
    my $nSeqSaved = 0;
    my $nDup = 0;
    open(IN, "<", $infile) || die "Cannot read from $infile";
    open(FAA, ">", "$outpre.faa") || die "Cannot write to $outpre.faa";
    open(QUERIES, ">", "$outpre.queries") || die "Cannot write to $outpre.queries";
    foreach my $line (<IN>) {
        chomp $line;
        my ($organism, undef, $queryId, $aaseq, $desc) = split /\t/, $line;
        die "Not enough columns in $infile" unless defined $aaseq;
        $desc = "" if !defined $desc; # allow it to be missing
        # "*" can be in the translations from MicrobesOnline because of errors
        # (like a programmed frame shift that is not recorded)
        # Also there may be a trailing translated stop codon
        $aaseq =~ s/[*]//g;
        die "Invalid input in $infile" if $organism eq "" || $queryId eq "" || $aaseq eq "";
        $aaseq =~ s/-//g; # a rare occurence in MicrobesOnline
        die "Invalid sequence $aaseq for $queryId" unless $aaseq =~ m/^[A-Z]+$/;
        if (exists $hasHit{$queryId}) {
            if (exists $hasSeq{$queryId}) {
                $nDup++;
            } else {
                print FAA ">$queryId\n$aaseq\n";
                $hasSeq{$queryId} = 1;
                $nSeqSaved++;
            }
            print QUERIES join("\t", $queryId, $organism, length($aaseq), $desc)."\n";
        }
    }
    close(IN) || die "Error reading $infile";
    close(FAA) || die "Error writing to $outpre.faa";
    close(QUERIES) || die "Error writing to $outpre.queries";
    my $nNoSeq = scalar(keys %hasHit) - $nSeqSaved;
    print STDERR "Warning: No sequence for $nNoSeq queries that have hits\n"
        if $nNoSeq > 0;
    print STDERR "Queries with hits and more than one locus tag entry: $nDup\n"
        if $nDup > 0;
    print STDERR "Wrote $nSeqSaved entries to\n$outpre.faa\n$outpre.queries\n";
}

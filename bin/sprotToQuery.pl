#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my $usage = <<END
sprotToQuery.pl -words wordsfile < uniprot_fasta > queryprot
	Given a fasta file from trembl or swissprot, with headers like
	>sp|accession|entry_name description OS=organism name KEY=...
	or
	>tr|accession|entry_name description OS=organism name KEY=...

	and a list of words (1 per line), outputs a set of queries.
END
    ;

sub ProcessFasta($$$); # header, sequence, ref to hash of words

{
    my ($wordsfile);
    die $usage
        unless GetOptions('words=s' => \$wordsfile) 
        && @ARGV == 0;
    die "No such file: $wordsfile\n" unless -e $wordsfile;
    my %words = ();
    open(WORDS, "<", $wordsfile) || die "Error reading $wordsfile\n";
    while (my $line = <WORDS>) {
        chomp $line;
        $words{$line} = 1;
    }
    close(WORDS) || die "Error reading $wordsfile";
    print STDERR "Read " . scalar(keys %words) . " words from $wordsfile\n";

    my $nSeq = 0;
    my $header = undef;
    my $seq = "";
    while(my $line = <STDIN>) {
        chomp $line;
        if ($line =~ m/^>(.*)$/) {
            $nSeq++;
            my $newheader = $1;
            &ProcessFasta($header, $seq, \%words) if defined $header;
            $header = $newheader;
            $seq = "";
        } else {
            die "Invalid sequence line $line" unless $line =~ m/^[A-Z]+$/;
            $seq .= $line;
        }
    }
    &ProcessFasta($header, $seq, \%words) if defined $header;
    print STDERR "Read $nSeq sequences from STDIN\n";
}

sub ProcessFasta($$$) {
    my ($header, $seq, $words) = @_;
    die "Empty sequence for: $header" if $seq eq "";
    $header =~ m/^(sp|tr)[|]([0-9A-Z]+)[|]([0-9A-Z_]+) (.*) OS=(.*)$/
        || die "Cannot parse header $header";
    my ($source, $acc,$entry_name,$desc,$org) = ($1,$2,$3,$4,$5);
    $org =~ s/^uncultured //;
    $org =~ s/ [A-Z]+=.*//;
    print join("\t", $org, $acc, $acc, $seq, $desc)."\n"
        if exists $words->{$acc};
    print join("\t", $org, $entry_name, $acc, $seq, $desc)."\n"
        if exists $words->{$entry_name};
}


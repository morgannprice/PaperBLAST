#!/usr/bin/perl -w
# Given a set of words that appear in a corpus,
# find the relevant tags in RefSeq
use strict;

# hash of CDS attributes, organism, and hash of locus_tags of interest
sub WriteCDS($$$);

die "Usage: findRefSeqQueries.pl words < refseq.gbff > output\n"
    unless @ARGV == 1;

my $nCDS = 0;
my $nCDSKeep = 0;

{
    my ($wordsfile) = @ARGV;
    open(WORDS, "<", $wordsfile) || die "Cannot read $wordsfile\n";
    my %words = ();
    foreach my $line (<WORDS>) {
        chomp $line;
        die "Invalid input line\n$line\n -- must have just one word per line" unless $line =~ m/^\S+$/;
        $words{$line} = 1;
    }
    close(WORDS) || die "Error reading $wordsfile";
    print STDERR "Read " . scalar(keys %words) . " words from $wordsfile\n";

    my $inCDS = 0;
    my $cds = {};
    my $contline = 0;
    my $org = "";
    my ($field,$value);
    while(my $line = <STDIN>) {
        chomp $line;
        if ($line =~ m/^     CDS/) {
            &WriteCDS($cds,$org,\%words) if $inCDS;
            $inCDS = 1;
            $cds = {};
        } elsif ($line =~ m/^     [a-zA-Z]/ || $line eq "//") {
            &WriteCDS($cds,$org,\%words) if $inCDS;
            $cds = {};
            $inCDS = 0;
            $org = "" if $line eq "//";
        } elsif ($inCDS) {
            if ($line =~ m!^                     /(.*)="(.*)"$!) {
                $field = $1;
                $value = $2;
                $cds->{$field} = $value;
            } elsif ($line =~ m!^                     /(.*)=["](.*)$!) {
                # unterminated --multi-line value
                $field = $1;
                $value = $2;
                $contline = 1;
            } elsif ($contline && $line =~ m/^                     (.*)$/) {
                my $add = $1;
                print STDERR "Unexpected in continuation context:\t$line\n"
                     if $add =~ m!^/[A-Za-z_]+=!;
                my $is_last = $add =~ m/["]$/;
                $add =~ s/["]$// if $is_last;
                $value .= " " unless $field eq "translation";
                $value .= $add;
                if ($is_last) {
                    $contline = 0;
                    $cds->{$field} = $value;
                }
            } elsif ($contline) {
                die "Unexpected line in continuation context\n$line\n";
            }
        } elsif ($line =~ m/^  ORGANISM +(.*)$/) {
            $org = $1;
        }
    }
    &WriteCDS($cds,$org,\%words) if $inCDS;
    print STDERR "Read $nCDS CDS, kept $nCDSKeep\n";
}

sub WriteCDS($$$) {
    my ($hash,$org,$words) = @_;
    my @words = (); # the words that appear in the corpus and may refer to this CDS
    foreach my $field (qw{old_locus_tag locus_tag gene gene_synonym protein_id}) {
        my $value = $hash->{$field};
        next unless defined $value;
        # ignore trailing version number (like .1 or .2) on protein ids
        $value =~ s/[.]\d+$// if $field eq "protein_id";
        push @words, $value if exists $words->{$value};
    }
    $nCDS++;
    if (scalar(@words) > 0) {
        $nCDSKeep++;
        foreach my $word (@words) {
            print join("\t",
                       $word, 
                       "organism=$org",
                       map { $_ . "=" . $hash->{$_} } sort keys %$hash)
                ."\n";
        }
    }
}

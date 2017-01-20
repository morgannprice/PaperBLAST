#!/usr/bin/perl -w
#######################################################
## showAlign.cgi
##
## Copyright (c) 2017 University of California
##
## Authors: Morgan Price
#######################################################
#
# Required CGI parameters:
# seq1 -- 1st sequence
# Either seq2 -- 2nd sequence or acc2 -- accession for 2nd sequence, to look it up in the database
#
# Optional CGI garameters:
# def1 -- name of 1st sequence
# def2 -- name of 2nd sequence

use strict;
use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use Time::HiRes qw{gettimeofday};

sub fetch_fasta($$);

my $cgi=CGI->new;
my $def1 = $cgi->param('def1') || "query";
my $def2 = $cgi->param('def2') || "subject";
my $title = "Align $def1 vs. $def2";
print
    header(-charset => 'utf-8'),
    start_html($title);

print h2("Align $def1 vs. $def2"), "\n";

my $tmpDir = "../tmp";
my $bl2seq = "../bin/blast/bl2seq";
die "No such executable: $bl2seq" unless -x $bl2seq;
my $fastacmd = "../bin/blast/fastacmd";
die "No such executable: $fastacmd" unless -x $fastacmd;

my $procId = $$;
my $timestamp = int (gettimeofday() * 1000);
my $prefix = "$tmpDir/$procId$timestamp";

my $seq1 = $cgi->param('seq1') || die "No seq1 argument";
my $seq2 = $cgi->param('seq2');
if (!defined $seq2) {
    my $acc2 = $cgi->param('acc2') || die "Either seq2 or acc2 must be specified";
    # Find it in the database using fastacmd
    my $db = "../data/litsearch.faa";
    die "No such file: $db" unless -e $db;
    $seq2 = fetch_fasta($acc2, $db) || die "Cannot fetch sequence for $def2";
}
$seq1 =~ m/^[A-Z*]+$/ || die "Invalid seq1";
$seq1 =~ s/[*]//g;
$seq2 =~ m/^[A-Z*]+$/ || die "Invalid seq2";
$seq2 =~ s/[*]//g;
die "Empty seq1" if $seq1 eq "";
die "Empty seq2" if $seq2 eq "";

my $len1 = length($seq1);
my $len2 = length($seq2);
my $cdd_base = "http://www.ncbi.nlm.nih.gov/Structure/cdd/wrpsb.cgi?seqinput=>";
my $newline = "%0A";
print
    p("$def1: $len1 amino acids", br(),
      "$def2: $len2 amino acids"),
    p("Or search the Conserved Domains Database with",
      a({ -href => "$cdd_base$def1$newline$seq1" }, $def1),
      "or",
      a({ -href => "$cdd_base$def2$newline$seq2" }, $def2));

open(FAA1, ">", "$prefix.faa1");
print FAA1 ">$def1\n$seq1\n";
close(FAA1) || die "Error writing to $prefix.faa1";

open(FAA2, ">", "$prefix.faa2");
print FAA2 ">$def2\n$seq2\n";
close(FAA2) || die "Error writing to $prefix.faa2";

system($bl2seq,
       "-p", "blastp",
       "-e", 0.01,
       "-F", "m S",
       "-i", "$prefix.faa1",
       "-j", "$prefix.faa2",
       "-o", "$prefix.bl2seq") == 0 || die "bl2seq failed";
my @lines = ();
open(IN, "<", "$prefix.bl2seq") || die "Cannot read $prefix.bl2seq";
@lines = <IN>;
@lines = map { chomp; $_; } @lines;
close(IN) || die "Error reading $prefix.bl2seq";

unlink("$prefix.faa1");
unlink("$prefix.faa2");
unlink("$prefix.bl2seq");

# remove lines up to the first alignment ("Score")
do {
    shift @lines
} until $lines[0] =~ m/Score/i || $lines[0] =~ m/No hits/i || @lines == 1;
my @out = ();
foreach my $line (@lines) {
    if ($line =~ m/^Lambda\s+K\s+H$/) {
        last;
    } else {
        push @out, $line;
    }
}

print
    h3("Protein alignments from bl2seq"),
    pre(join("\n", @out));

sub fetch_fasta($$) {
    my ($def, $db) = @_;
    unlink("$prefix.fetch");
    die "Invalid def2 argumnet: $def" if $def eq "" || $def =~ m/\s/ || $def =~ m/,/;
    system("$fastacmd","-s",$def,"-d",$db,"-o", "$prefix.fetch");
    open(SEQ, "<", "$prefix.fetch") || die "Cannot read $prefix.fetch -- fastacmd failed?";
    my @lines = <SEQ>;
    close(SEQ) || die "Error reading $prefix.fetch";
    unlink("$prefix.fetch");
    (@lines > 0 && $lines[0] =~ m/^>/) || die "Unknown def: $def";
    shift @lines;
    @lines = map { chomp; $_; } @lines;
    return join("", @lines);
}

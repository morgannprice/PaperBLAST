#!/usr/bin/perl -w

use strict;
use XML::LibXML;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use Miner; # for TextSnippets()

my $snippetBefore = 15;
my $snippetAfter = 30;
my $maxChar = 1000;

my $usage = <<END
abstractSnippets.pl [ -in abstracts | -list files ] -out pubmed.snippets oa.papers ... refseq.papers

	The input file(s) should be tab-delimited with the 1st column
	having the pubmedId and the second column having the abstract
	(as from pubmedparse.pl)

	The papers files are from parseEuropePMCHits.pl

	Output will be tab-delimited with the fields
	pmcId, pmId, queryTerm, queryId, snippet

	Also writes to output_file.access, with tab-delimited records of the form
	pmcId, pmId, "abstract"
	if it was able to read the abstract for the paper

	Optional arguments control the size of each snippet in words or characters:
	-before $snippetBefore
        -after $snippetAfter
        -maxChar $maxChar

	There is a limit of 2 snippets per gene per paper
END
    ;

sub ProcessLine($$$);

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

    foreach my $paperfile (@paperfiles) {
        die "No such file: $paperfile\n" unless -e $paperfile;
    }
    my %papers = (); # pmId => queryTerm => list of queryId
    my %pm2pmc = (); # pmId => pmcId
    foreach my $paperfile (@paperfiles) {
        open(PAPERS, "<", $paperfile) || die "Error reading $paperfile";
        while(my $line = <PAPERS>) {
            chomp $line;
            my ($queryId, $queryTerm, $pmcId, $pmId) = split /\t/, $line;
            die "Not enough fields in line\n$line\nin $paperfile" unless defined $pmId;
            next if $pmId eq "";
            $queryTerm =~ s/^"//;
            $queryTerm =~ s/"$//;
            push @{ $papers{$pmId}{$queryTerm} }, $queryId;
            $pm2pmc{$pmId} = $pmcId;
        }
        close(PAPERS) || die "Error reading $paperfile";
    }
    print STDERR "Read " . scalar(keys %papers) . " pubmedIds to search\n";

    my $nRelevant = 0;
    open(OUT, ">", $outfile) || die "Error writing to $outfile";
    my $accfile = "$outfile.access";
    open(ACC, ">", $accfile) || die "Error writing to $accfile";
    if (defined $infile) {
      open(IN, "<", $infile) || die "Cannot read $infile";
      while(my $line = <IN>) {
        $nRelevant++ if ProcessLine($line, \%papers, \%pm2pmc);
      }
      close(IN) || die "Error reading $infile";
    } elsif (defined $listfile) {
      open(LIST , "<", $listfile) || die "Cannot read $listfile";
      while(my $file = <LIST>) {
        chomp $file;
        open(IN, "<", $file) || die "Cannot read $file";
        while(my $line = <IN>) {
          $nRelevant++ if ProcessLine($line, \%papers, \%pm2pmc);
        }
        close(IN) || die "Error reading $file";
      }
      close(LIST) || die "Error reading $listfile";
    }
    close(OUT) || die "Error writing to $outfile";
    close(ACC) || die "Error writing to $accfile";
    print STDERR "Processed $nRelevant relevant abstracts\n";
}

sub ProcessLine($$$) {
  my ($line,$papers,$pm2pmc) = @_;
  chomp $line;
  my ($pmId, $abstract) = split /\t/, $line;
  return 0 unless $abstract && exists $papers->{$pmId};
  my $hash = $papers->{$pmId};
  my $pmcId = $pm2pmc->{$pmId};
  print ACC join("\t", $pmcId, $pmId, "abstract")."\n";
  while (my ($queryTerm, $queryIds) = each %$hash) {
    if ($abstract =~ m/$queryTerm/) { # require case-sensitive match
      my @snippets = TextSnippets($abstract, $queryTerm, $snippetBefore, $snippetAfter, $maxChar);
      my %queryIds = map { $_ => 1 } @$queryIds;
      my @queryIds = keys(%queryIds);
      foreach my $queryId (@queryIds) {
        foreach my $snippet (@snippets) {
          print OUT join("\t", $pmcId, $pmId, $queryTerm, $queryId, $snippet)."\n";
        }
      }
    }
  }
  return 1;
}

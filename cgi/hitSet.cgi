#!/usr/bin/perl -w
# Required CGI parameter: set -- set.fasta must be a file in ../tmp/sets/
# Optional parameters:
# minIdentity -- default 40
# minCoverage -- default 75
# minYear -- defaults to last year

use strict;
use CGI qw(:standard Vars start_ul);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use Time::HiRes qw{gettimeofday};
use DBI;
use IO::Handle; # for autoflush
use lib "../lib";
use pbweb;
use pbutils qw{ReadFastaDesc NewerThan};
use POSIX qw(strftime);
use URI::Escape;


my $blastall = "../bin/blast/blastall";
my $nCPU = 20;
my $base = "../data";
my $blastdb = "$base/uniq.faa";
my $sqldb = "$base/litsearch.db";
my $dbh = DBI->connect("dbi:SQLite:dbname=$sqldb","","",{ RaiseError => 1 }) || die $DBI::errstr;

die "No such executable: $blastall" unless -x $blastall;
die "No such file: $blastdb" unless -e $blastdb;
die "No such file: $sqldb" unless -e $sqldb;

my $cgi=CGI->new;
my $set = $cgi->param('set') || die "Must specify the set argument";
die "Invalid set" unless $set =~ m/^[a-zA-Z0-9_.-]+$/;
my $setPre = "../tmp/sets/$set";
my $faaFile = "$setPre.fasta";
die "No such file: $faaFile" unless -e $faaFile;

start_page('title' => "Recent papers about $set");
print <<END
<SCRIPT src="../static/pb.js"></SCRIPT>
END
;

autoflush STDOUT 1; # show preliminary results

my %h = ReadFastaDesc($faaFile);
die "Error reading $faaFile" if exists $h{error};
my $seqs = $h{seq} || die;
my $descs = $h{desc} || die;
my $lens = $h{len} || die;

my $hitsFile = "$setPre.hits";
unless (NewerThan($hitsFile, $blastdb)) {
  print p("Running BLAST; this might take a few minutes"), "\n";
  my $tmpFile = "$hitsFile.$$.tmp";
  system("$blastall", "-p", "blastp",  "-d", $blastdb, "-i", $faaFile, "-o", $tmpFile,
         "-e", "1e-5", "-m 8", "-a", $nCPU, -F, "m S") == 0
           || die "Error running $blastall: $!";
  print p("BLAST succeeded"), "\n";
  rename($hitsFile, "$hitsFile.old"); # save old hits
  rename($tmpFile, $hitsFile) || die "Rename failed";
}

my $minIdentity = $cgi->param('minIdentity');
$minIdentity = 40 if !defined $minIdentity || $minIdentity eq "";
die "Invalid minIdentity" unless $minIdentity =~ m/^\d+$/;

my $minCoverage = $cgi->param('minCoverage');
$minCoverage = 75 if !defined $minCoverage || $minCoverage eq "";
die "Invalid minCoverage" unless $minCoverage =~ m/^\d+$/;

my %hits = (); # id to list of rows
open(my $fhHits, "<", $hitsFile) || die "Cannot read $hitsFile";
while (my $line = <$fhHits>) {
  chomp $line;
  my @row = split /\t/, $line;
  my ($queryId,$subjectId,$percIdentity,$alnLength,$mmCnt,$gapCnt,
      $queryStart,$queryEnd,$subjectStart,$subjectEnd,
      $eVal,$bitscore) = @row;
  die "Unknown query $queryId" unless exists $lens->{$queryId};
  next unless $percIdentity >= $minIdentity
    || ($queryEnd-$queryStart+1) >= $lens->{queryId} & $minCoverage/100.0;
  push @{ $hits{$queryId} }, \@row;
}

print p("Proteins with a hit at ${minIdentity}% identity and ${minCoverage}% coverage:",
        scalar(keys %hits), "of", scalar(keys %$lens)), "\n";

# Now, show which ones have hits, and link to the full results at litSearch.cgi

my $thisYear = strftime("%Y", localtime);
my $minYear = $cgi->param('minYear');
$minYear = $thisYear - 1 if !defined $minYear || $minYear eq "";
die "Invalid minYear" unless $minYear =~ m/^\d+$/;

my %hasRecent = (); # names with a paper

foreach my $queryId (sort keys %hits) {
  foreach my $row (@{ $hits{$queryId} }) {
    my $hasPaper = 0;
    my $uniqId = $row->[1];
    my $dups = $dbh->selectcol_arrayref("SELECT duplicate_id FROM SeqToDuplicate WHERE sequence_id = ?",
                                      {}, $uniqId);
    my @subjectIds = $uniqId;
    push @subjectIds, @$dups;
    my $keep = 0;
    foreach my $subjectId (@subjectIds) {
      my $years = $dbh->selectcol_arrayref(qq{ SELECT DISTINCT year FROM GenePaper WHERE geneId = ? },
                                           {}, $subjectId);
      foreach my $year (@$years) {
        if ($year >= $minYear) {
          $keep = 1;
          last;
        }
      }
      last if $keep;
    }
    $hasRecent{$queryId} = 1 if $keep;
  }
}

print p("Proteins with papers from $minYear or later:", scalar(keys %hasRecent)), "\n";

if (keys(%hasRecent) > 0) {
  print start_ul();
  foreach my $name (sort keys %hasRecent) {
    my $encodedQuery = uri_escape(">$name $descs->{$name}\n$seqs->{$name}\n");
    print li(a({ -href => "litSearch.cgi?query=$encodedQuery", -title => "see PaperBLAST hits" },
               $name), $descs->{$name});
  }
  print end_ul();
}

finish_page();

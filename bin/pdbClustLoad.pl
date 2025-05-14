#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw{$RealBin};
use lib "$RealBin/../lib";
use pbutils qw{ReadFastaEntry};
sub csv_quote($);

my $usage = <<END
Usage: pdbClustLoad.pl -in pdbClust.faa -db litsearch.db

Empties and refills the PdbClustInfo table in the sqlite3 database, using the headers from pdbClust.faa
END
  ;

my ($inFile, $dbFile);
die $usage
  unless GetOptions('inFile=s' => \$inFile, 'db=s' => \$dbFile)
  && @ARGV == 0
  && defined $inFile && defined $dbFile;

die "No such file: $dbFile\n" unless -e $dbFile;
open(my $fh, "<", $inFile) || die "Cannot read $inFile\n";
my $infoFile = "$inFile.info";
open(my $fhInfo, ">", $infoFile) || die "Cannot write to $infoFile";
my $state = {};
my $nRows = 0;
while(my ($header, undef) = ReadFastaEntry($fh, $state)) {
  $header =~ m/^([0-9a-zA-Z]+)_([0-9a-zA-Z]+) (.*)$/ || die "Cannot parse header: $header";
  my ($id, $chain, $desc) = ($1,$2,$3);
  print $fhInfo join("\t", $id, $chain, csv_quote($desc))."\n";
  $nRows++;
}
close($fhInfo) || die "Error writing to $infoFile";
close($fh) || die "Error reading from $inFile";

open(SQLITE, "|-", "sqlite3", $dbFile) || die "Cannot run sqlite3 on $dbFile";
print SQLITE <<END
.mode tabs
DELETE FROM PdbClustInfo;
.import $inFile.info PdbClustInfo
.quit
END
  ;
close(SQLITE) || die "Load into PdbClustInfo failed";

my $nActual = `sqlite3 $dbFile 'SELECT COUNT(*) FROM PdbClustInfo;'`;
chomp $nActual;
die "Wrong number of entries in PdbClustInfo: $nActual vs. expected $nRows\n"
  unless $nActual == $nRows;
print STDERR "The PdbClustInfo table now has $nActual rows\n";

sub csv_quote($) {
  my ($in) = @_;
  return $in unless $in =~ m/"/;
  $in =~ s/"/""/g;
  return '"' . $in . '"';
}

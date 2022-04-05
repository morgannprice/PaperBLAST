#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Digest::MD5 qw{md5_hex};
use FindBin qw{$RealBin};
use lib "$RealBin/../lib";
use pbutils qw{ReadFastaEntry};

my $indir = ".";
my $usage = <<END
Usage: loadSiteTables.pl -db sites.db

Optional arguments:
  -indir $indir -- the directory with the input
    tab-delimited files
    (Sites.tab, HasSites.tab, PDBLigands.tab)

The schema should already have been set up (using litsearch.sql)
END
;

my $dbfile;
die $usage
  unless GetOptions('db=s' => \$dbfile,
                    'indir=s' => \$indir);
die "Not a directory: $indir\n" unless -d $indir;
die "Must specify -db:\n$usage" unless defined $dbfile && $dbfile ne "";
die "database file $dbfile must already exist (with schema loaded)\n"
  unless -s $dbfile;

my %siteSeq = (); # db => id => chain => sequence
open(my $fhSeq, "<", "$indir/hassites.faa")
  || die "Cannot read $indir/hassites.faa";
my $state = {};
while (my ($header,$seq) = ReadFastaEntry($fhSeq,$state)) {
  my ($db, $id, $chain) = split /:/, $header;
  die "Invalid header $header" unless defined $id;
  $chain = "" if !defined $chain;
  $siteSeq{$db}{$id}{$chain} = $seq;
}
close($fhSeq) || die "Error reading $indir/hassites.faa";

my %hasSite = (); # db => id => chain => 1
open(my $fhSites, "<", "$indir/Sites.tab")
  || die "Cannot read $indir/Sites.tab";
while(my $line = <$fhSites>) {
  chomp $line;
  my ($db, $id, $chain) = split /\t/, $line;
  die "$db:$id:$chain has sites but no sequence"
    unless exists $siteSeq{$db}{$id}{$chain};
  $hasSite{$db}{$id}{$chain} =1 ;
}
close($fhSites) || die "Error reading $indir/Sites.tab";

my $tmpDir = $ENV{TMP} || "/tmp";
my $tmpFile = "$tmpDir/SeqHasSite.$$.tab";

open(my $fhOut, ">", $tmpFile)
  || die "Cannot write to $tmpFile";
foreach my $db (sort keys %hasSite) {
  my $hash = $hasSite{$db};
  foreach my $id (sort keys %$hash) {
    my @chains = sort keys %{ $hasSite{$db}{$id} };
    foreach my $chain (@chains) {
      my $seq = $siteSeq{$db}{$id}{$chain} || die;
      my $hash = md5_hex($seq);
      print $fhOut join("\t", $hash, length($seq), $db, $id, $chain)."\n";
    }
  }
}
close($fhOut) || die "Error writing $tmpFile";

open(my $fhSQL, "|-", "sqlite3", $dbfile)
  || die "Cannot run sqlite3 on $dbfile\n";
print $fhSQL <<END
.mode tabs
.import $indir/Sites.tab Site
.import $indir/HasSites.tab HasSites
.import $indir/PDBLigands.tab PDBLigand
.import $tmpFile SeqHasSite
END
;
close($fhSQL) || die "Error running loading commands for sqlite3 on $dbfile\n";
unlink($tmpFile);
print STDERR "Loaded site tables into $dbfile from ${indir}/\n";

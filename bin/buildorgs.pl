#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw{$Bin};
use lib "$Bin/../lib";
use Steps qw{WriteAssemblyAsOrgProteins WriteSixFrameAsOrgProteins};
use FetchAssembly qw{CacheAssembly SetFitnessBrowserPath SetPrivDir};
use pbutils qw{ReadFastaEntry};

my $cachedir = "PaperBLAST/tmp";
my $privdir = "PaperBLAST/private"; # for JGI access key
my $febadata = "feba/cgi_data";
my $usearch = "./usearch";

my $usage = <<END
./buildorgs.pl -out prefix -orgs gdb1:gid1 ... gdbN:gidN

where gdb is one of the genome databases supported by FetchAssembly
(NCBI, JGI, UniProt, MicrobesOnline, FitnessBrowser)
and gid is a genome identifier.

Writes to prefix.faa and prefix.org

Optional arguments:

-sixframe -- ignore the protein annotation and use the 6-frame
 translation instead
-cache cachedir -- the directory for caching the fetched genomes
-febadata febadata -- the feba/cgi_data directory that
  contains feba.db and aaseqs
-privdir privdir -- the directory with the JGI key
-usearch $usearch -- the location of the usearch executable
END
;

my $prefix;
my @orgspec;
my $sixframe;

die $usage
  unless GetOptions('orgs=s{1,}' => \@orgspec,
                    'out=s' => \$prefix,
                    'cache=s' => \$cachedir,
                    'febadata=s' => \$febadata,
                    'privdir=s' => \$privdir,
                    'usearch=s' => \$usearch,
                    'sixframe' => \$sixframe)
  && @ARGV == 0;
die "Must specify -out and -orgs\n"
  unless $prefix && @orgspec > 0;
foreach my $dir ($cachedir, $febadata, $privdir) {
  die "Not a directory: $dir\n"
    unless -d $dir;
}
SetFitnessBrowserPath($febadata);
SetPrivDir($privdir);
my $blastdir = "$Bin/blast";
die "No such directory: $blastdir\n" unless -d $blastdir;
my $formatdb = "$blastdir/formatdb";
die "No such executable: $formatdb\n" unless -x $formatdb;

# Fetch the assemblies
my %gid = (); # gid => assembly
my @assemblies = ();
foreach my $orgspec (@orgspec) {
  die "Whitespace not allowed in organism specifier: $orgspec\n"
    if $orgspec =~ m/\s/;
  my @parts = split /:/, $orgspec;
  die "organism specifier $orgspec should have two parts"
    unless @parts == 2;
  my ($gdb, $gid) = @parts;
  die "Duplicate genome identifier $gid -- each genome can be included only once\n"
    if exists $gid{$gid};
  my $assembly = CacheAssembly($gdb, $gid, $cachedir)
    || die "Cannot fetch $orgspec\n";
  print STDERR "Loaded $orgspec\n";
  push @assemblies, $assembly;
  $gid{$gid} = $assembly;
}

open(my $fhOut, ">", "$prefix.faa")
  || die "Cannot write to $prefix.faa\n";
foreach my $assembly (@assemblies) {
  if (defined $sixframe) {
    $assembly->{orgId} = WriteSixFrameAsOrgProteins($assembly, $fhOut, $usearch);
  } else {
    $assembly->{orgId} = WriteAssemblyAsOrgProteins($assembly, $fhOut);
  }
}
close($fhOut) || die "Error writing to $prefix.faa\n";
print STDERR "Wrote $prefix.faa\n";

# Index $prefix.faa
system($formatdb, "-p", "T", "-i", "$prefix.faa", "-o", "T") == 0
  || die "$formatdb on $prefix.faa failed: $!\n";

open(my $fhOrg, ">", "$prefix.org")
  || die "Cannot write to $prefix.org\n";
print $fhOrg join("\t", qw{orgId gdb gid genomeName nProteins})."\n";
foreach my $assembly (@assemblies) {
  my $faaFile = $assembly->{faafile} || die "No faa file";
  open(my $fhIn, "<", $faaFile) || die "Cannot read $faaFile";
  my $nSeq = 0;
  my $state = {};
  while (my ($header, $seq) = ReadFastaEntry($fhIn, $state)) {
    $nSeq++;
  }
  close($fhIn) || die "Error reading $faaFile\n";
  # assembly->{orgId} is from WriteAssemblyAsOrgProteins()
  print $fhOrg join("\t", $assembly->{orgId}, $assembly->{gdb}, $assembly->{gid},
                    $assembly->{genomeName}, $nSeq)."\n";
}
close($fhOrg) || die "Error writing to $prefix.faa\n";
print STDERR "Wrote $prefix.org\n";

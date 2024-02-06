#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw{$Bin};
use lib "$Bin/../lib";
use Steps qw{WriteAssemblyAsOrgProteins WriteSixFrameAsOrgProteins};
use FetchAssembly qw{CacheAssembly SetFitnessBrowserPath SetPrivDir AASeqToAssembly};
use pbutils qw{ReadFastaEntry ReadFastaDesc};

my $cachedir = "$Bin/../tmp";
my $privdir = "$Bin/../private"; # for JGI access key
my $febadata = "$Bin/../fbrowse_data";
my $usearch = "$Bin/usearch";

my $usage = <<END
./buildorgs.pl -out prefix -orgs gdb1:gid1 ... gdbN:gidN
./buildorgs.pl -out prefix -orgfile file

gdb is one of the genome databases supported by FetchAssembly
(NCBI, JGI, UniProt, MicrobesOnline, FitnessBrowser)
and gid is a genome identifier,
 or gdb=file and gid=filename
 or gdb=file and gid=filename:genomeName

If using -orgfile -- the input file should contain one organism
specifier per line, or comment lines that begin with "#".

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
my $orgfile;
my $sixframe;

die $usage
  unless GetOptions('orgs=s{1,}' => \@orgspec,
                    'orgfile=s' => \$orgfile,
                    'out=s' => \$prefix,
                    'cache=s' => \$cachedir,
                    'febadata=s' => \$febadata,
                    'privdir=s' => \$privdir,
                    'usearch=s' => \$usearch,
                    'sixframe' => \$sixframe)
  && @ARGV == 0;
die "Must specify -out and -orgs or -orgsfile\n"
  unless $prefix && (@orgspec > 0 || defined $orgfile);
die "Cannot specify -orgs and -orgsfile\n" if @orgspec > 0 && defined $orgfile;
foreach my $dir ($cachedir, $febadata, $privdir) {
  die "Not a directory: $dir\n"
    unless -d $dir;
}
if (defined $orgfile) {
  open(my $fhOrgs, "<", $orgfile) || die "Cannot read $orgfile\n";
  while (my $line = <$fhOrgs>) {
    chomp $line;
    push @orgspec, $line
      unless $line eq "" || $line =~ m/^#/;
  }
  close($fhOrgs) || die "Error reading $orgfile\n";
  die "orgfile $orgfile contains no org specifiers\n"
    if @orgspec == 0;
}

SetFitnessBrowserPath($febadata);
SetPrivDir($privdir);

# Fetch the assemblies
FetchAssembly::setFailMode("warning");
my %gid = (); # gid => assembly
my @assemblies = ();
my @fail = ();
foreach my $orgspec (@orgspec) {
  my @parts = split /:/, $orgspec;
  die "organism specifier $orgspec should have two parts"
    unless @parts == 2 || (@parts ==3 && $parts[0] eq "file");
  my ($gdb, $gid) = @parts;
  die "Whitespace not allowed in organism specifier: $orgspec\n"
    if $gdb =~ m/\s/ || $gid =~ m/\s/;
  die "Duplicate genome identifier $gid -- each genome can be included only once\n"
    if exists $gid{$gid};
  my $assembly;
  if ($gdb eq "file") {
    my $genomeName = $parts[2] if @parts > 2;
    die "No such file: $gid\n" unless -e $gid;
    my %faa = ReadFastaDesc($gid);
    die "Invalid fasta input $gid: " . values(%faa)
      if exists $faa{error};
    my %headerToSeq = ();
    while (my ($id, $desc) = each %{ $faa{desc} }) {
      my $key = "$id $desc";
      die "Duplicate sequence for $key in $gid\n" if exists $headerToSeq{$key};
      $headerToSeq{$key} = $faa{seq}{$id} || die "Empty sequence for $id in $gid\n";
    }
    my @ids = sort keys %{ $faa{seq} };
    my $n = scalar(@ids);
    $genomeName = "$ids[0]:$ids[-1] ($n proteins)"
      unless defined $genomeName;
    $assembly = AASeqToAssembly(\%headerToSeq, $cachedir);
    $assembly->{genomeName} = $genomeName;
  } else {
    $assembly = CacheAssembly($gdb, $gid, $cachedir);
    unless($assembly) {
      push @fail, "$gdb:$gid";
      print STDERR "Warning: cannot fetch $orgspec\n";
    }
  }
  print STDERR "Loaded $orgspec\n";
  push @assemblies, $assembly;
  $gid{$gid} = $assembly;
}
die "Some assemblies failed to load: " . join(" ", @fail) . "\n" if @fail > 0;

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
close($fhOrg) || die "Error writing to $prefix.org\n";
print STDERR "Wrote $prefix.org\n";

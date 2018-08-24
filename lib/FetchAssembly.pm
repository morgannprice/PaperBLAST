# Fetch information about assemblies from NCBI

package FetchAssembly;
require Exporter;
use strict;
use LWP::Simple; # for get()
use XML::LibXML;
use Time::HiRes qw{gettimeofday};

our (@ISA,@EXPORT);
@ISA = qw(Exporter);
@EXPORT = qw(FetchAssemblyInfo FetchAssemblyFaa FetchAssemblyFna);

my $maxfetch = 20;
sub GetMaxFetch() {
  return $maxfetch;
}

# One argument: the query
# Returns a list of results (up to maxfetch) as a hash; each "assembly" includes
# uid (the assembly uid), id (an identifier like GCF_000195755.1), ftp (the URL for the ftp site),
# and organism
sub FetchAssemblyInfo($) {
  my ($query) = @_;
  # First run the query using esearch
  my $URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=assembly&term=" . $query;
  my $string = get($URL);
  my $fxml = XML::LibXML->load_xml(string => $string, recover => 1);
  my @idnodes = $fxml->findnodes("//IdList/Id");
  my @ids = map { $_->textContent } @idnodes;
  @ids = grep { defined && $_ ne "" } @ids;
  return () if @ids == 0;

  # Now fetch metadata for these ids. Limit to 20.
  @ids = $ids[0..($maxfetch-1)] if scalar(@ids) > $maxfetch;
  my $idspec= join(",",@ids);
  $URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=assembly&id=$idspec";
  $string = get($URL);
  my $sxml = XML::LibXML->load_xml(string => $string, recover => 1);
  my @obj = $sxml->findnodes("//DocumentSummary");
  die "Failed to fetch summaries for assembly ids using NCBI's eutils: @ids\n"
    unless scalar(@ids) == scalar(@obj);
  my @out = ();
  foreach my $obj (@obj) {
    my $out = {};
    foreach my $attr ($obj->attributes) {
      $out->{uid} = $attr->getValue if $attr->nodeName eq "uid";
    }
    die "No uid attribute in results from NCBI's eutils for @ids\n"
      if !exists $out->{uid};
    my @part = $obj->findnodes("AssemblyAccession");
    $out->{id} = $part[0]->textContent if @part;
    @part = $obj->findnodes("SpeciesName");
    $out->{org} = $part[0]->textContent if @part;
    @part = $obj->findnodes("Biosource/InfraspeciesList/Infraspecie/Sub_value");
    foreach my $part (@part) {
      my $extra = $part->textContent;
      # Avoid adding the suffix if it is already there
      # which happens for organisms with names like Genus sp. strain
      $out->{org} .= " " . $extra
        if $extra ne "" && $out->{org} !~ m/ $extra/;
    }
    @part = $obj->findnodes("FtpPath_RefSeq");
    $out->{ftp} = $part[0]->textContent if @part;
    if (!defined $out->{ftp} || $out->{ftp} eq "") {
      @part = $obj->findnodes("FtpPath_GenBank");
      $out->{ftp} = $part[0]->textContent if @part;
    }
    push @out, $out;
  }
  return @out;
}

# Given the assembly hash, fetch protein fasta, uncompress it, and write it to the given file
# Writes to a temporary file to avoid overwriting by another process
# Returns success or failure
sub FetchAssemblyFaa($$) {
  my ($assembly, $outfile) = @_;
  return FetchAssemblyFileGz($assembly, "protein.faa", $outfile);
}

# Given the assembly hash, fetch protein fasta, uncompress it, and write it to the given file
# Writes to a temporary file to avoid overwriting by another process
# Returns success or failure
sub FetchAssemblyFna($$) {
  my ($assembly, $outfile) = @_;
  return FetchAssemblyFileGz($assembly, "genomic.fna", $outfile);
}

sub FetchAssemblyFileGz($$$) {
  my ($assembly, $suffix, $outfile) = @_;
  my $procId = $$;
  my $timestamp = int (gettimeofday() * 1000);
  my $tmpfile = $outfile . $procId . $timestamp;
  my $tmpfilegz = $outfile . $procId . $timestamp . ".gz";
  my $URL = $assembly->{ftp};
  return 0 unless $URL =~ m!/([^/]+)$!;
  my $prefix = $1;
  $URL = join("", $URL, "/", $prefix, "_", $suffix, ".gz");
   system("wget", "-nv", "-O", $tmpfilegz, $URL) == 0
    || return 0;
  system("gunzip", $tmpfilegz) == 0 || return 0;
  return rename($tmpfile, $outfile) || return 0;
}

1;


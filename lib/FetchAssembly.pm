# Fetch information about assemblies from NCBI, JGI, UniProt, MicrobesOnline, Fitness

package FetchAssembly;
require Exporter;
use strict;
use LWP::Simple; # for get()
use LWP::UserAgent;
use XML::LibXML;
use Time::HiRes qw{gettimeofday};
use pbutils;
use URI::Escape;
use HTTP::Cookies;
use JSON;

our (@ISA,@EXPORT);
@ISA = qw(Exporter);
@EXPORT = qw(FetchNCBIInfo FetchNCBIFaa FetchNCBIFna FetchNCBIFeatureFile ParseNCBIFeatureFile
             SearchJGI CreateJGICookie FetchJGI SearchUniProtProteomes UniProtProteomeInfo);

my $maxfetchNCBI = 20;
sub GetMaxFetchNCBI() {
  return $maxfetchNCBI;
}

# One argument: the query
# Returns a list of results (up to maxfetchNCBI) as a hash; each "assembly" includes
# uid (the assembly uid), id (an identifier like GCF_000195755.1), ftp (the URL for the ftp site),
# and org
sub FetchNCBIInfo($) {
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
  @ids = $ids[0..($maxfetchNCBI-1)] if scalar(@ids) > $maxfetchNCBI;
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
sub FetchNCBIFaa($$) {
  my ($assembly, $outfile) = @_;
  return FetchNCBIFileGz($assembly, "protein.faa", $outfile);
}

# Given the assembly hash, fetch protein fasta, uncompress it, and write it to the given file
# Writes to a temporary file to avoid overwriting by another process
# Returns success or failure
sub FetchNCBIFna($$) {
  my ($assembly, $outfile) = @_;
  return FetchNCBIFileGz($assembly, "genomic.fna", $outfile);
}

# Given the assembly hash, fetch the feature table, uncompress it, and write it to the given file
# Writes to a temporary file to avoid overwriting by another process
# Returns success or failure
sub FetchNCBIFeatureFile($$) {
  my ($assembly, $outfile) = @_;
  return FetchNCBIFileGz($assembly, "feature_table.txt", $outfile);
}

sub FetchNCBIFileGz($$$) {
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

# From a file to a reference to a list of hashes
sub ParseNCBIFeatureFile($) {
  my ($file) = @_;
  my @rows = pbutils::ReadTable($file, ["# feature", "class", "start", "end", "strand", "product_accession", "non-redundant_refseq", "name", "symbol", "locus_tag", "attributes"]);
  return \@rows;
}

# Parse comma-delimited lines.
# Fields that contain commas or " characters must quoted with "
#	In that case " represented as ""
# Not set up to handle linebreaks within values
# Returns an error (or "" on success) and a reference to a list of fields
sub parse_comma_delimited_line($) {
  my ($line) = @_;
  $line =~ s/[\r\n]+$//;
  my @out = ();
  my $saw_last = 0;
  while (length($line)) {
    if ($line =~ s/^([^,"]*),//) {
      # stripped off an unquoted field
      push @out, $1;
    } elsif ($line =~ s/^"(([^"]|"")+)",//) {
      # stripped off a quoted field that may contain ""
      my $quoted = $1;
      $quoted =~ s/""/"/g;
      push @out, $quoted;
    } elsif ($line =~ s/^([^,"]*)$//) {
      # this was the last field, unquoted
      push @out, $1;
      $saw_last = 1;
    } elsif ($line =~ s/^"(([^"]|"")+)"$//) {
      # this was the last field, quoted
      my $quoted = $1;
      $quoted =~ s/""/"/g;
      push @out, $quoted;
      $saw_last = 1;
    } else {
      return "Cannot parse remaining part $line", \@out;
    }
  }
  push @out, "" unless $saw_last;
  return "", \@out;
}

# Returns a string with an error message or the empty string on success,
# and optionally a list of hashes of the query results
# The key fields in the result set are the "Project Name" (the genone name)
# and portalId
sub SearchJGI($) {
  my ($query) = @_;
  return "", [] unless defined $query && $query ne "";
  $query = uri_escape($query);
  my $hostspec = "https://genome.jgi.doe.gov";
  my $URL = join("&", "$hostspec/portal/ext-api/search-service/export?core=genome",
                 "query=$query",
                 "searchIn=JGI%20Projects",
                 "searchType=Keyword",
                 "showAll=false",
                 "externallySequenced=true",
                 "sortBy=displayNameStr",
                 "showRestricted=false",
                 "showOnlyPublished=false",
                 "showSuperseded=true",
                 "sortOrder=asc",
                 "rawQuery=true",
                 "showFungalOnly=false",
                 "programName=all",
                 "programYear=all",
                 "superkingdom=--any--",
                 "status=--any--",
                 "scientificProgram=--any--",
                 "productCategory=--any--",
                 "selectedRecords=",
                 "start=0",
                 "rows=100");
  my $ua = LWP::UserAgent->new( ssl_opts => { SSL_verify_mode => 'SSL_VERIFY_NONE' } );
  my $response = $ua->get($URL);
  return "Error contacting $hostspec" unless $response->is_success;
  my $lines = $response->decoded_content( charset => "none" );
  my @header = ();
  my @out = ();
  foreach my $line (split /\n/, $lines) {
    my ($error, $pieces) = parse_comma_delimited_line($line);
    return "Results list is not in CSV format: $error; see $URL"
      if $error ne "";
    next if @$pieces == 0;
    if (@header == 0) {
      @header = @$pieces;
    } else {
      return "Wrong #columns in CSV results", \@out
        unless scalar(@header) == scalar(@$pieces);
      my $row = {};
      foreach my $i (0..(scalar(@header)-1)) {
        $row->{ $header[$i] } = defined $pieces->[$i] ? $pieces->[$i] : "";
      }
      # Compute the portal.id field
      if (exists $row->{"Portal ID"} && $row->{"Portal ID"} =~ m/"/) {
        $row->{"Portal ID"} =~ m/,"([0-9A-Za-z_]+)"[)]$/
          || return "Cannot find the portal id in " . $row->{"Portal ID"};
        $row->{portalId} = $1;
        $row->{genomeName} = $row->{"Project Name"};
        $row->{URL} = "http://genome.jgi.doe.gov/" . $row->{portalId};
        $row->{gdb} = "IMG";
        $row->{gid} = $row->{portalId};
        push @out, $row;
      }
    }
  }
  return "", \@out;
}

# Returns an error or the empty string
sub CreateJGICookie($$$) {
  my ($username, $passwd, $cookiefile) = @_;
  my $tmpfile = "$cookiefile.$$.tmp";
  unlink($tmpfile); # in case it exists
  my $URL = "https://signon-old.jgi.doe.gov/signon/create";
  my @cmd = ("curl", "--silent", $URL,
              "--data-urlencode", "login=$username",
              "--data-urlencode", "password=$passwd",
              "-c", $tmpfile);
  open(my $fh, "-|", @cmd)
    || return "Cannot run curl to $URL";
  while(<$fh>) {
    ;
  }
  close($fh) || die "Error running curl for $URL";
  die "Cookie file not created" unless -e $tmpfile;
  rename($tmpfile, $cookiefile) || die "renaming cookie file failed";
  return "";
}

# Given a portal id and a genome name, fetch the IMG content into a new directory
# Returns an error message or ""
# Creates an additional file "fields" that contains
# filename prefix (IMG genome id), genome name, and portal prefix
sub FetchJGI($$$$) {
  my ($portalid, $genomename, $cookiefile, $dir) = @_;
  return "Directory $dir already exists" if -e $dir;
  my $down_url = "https://genome.jgi.doe.gov/portal/ext-api/downloads/get-directory?organism=$portalid";
  my $in; # for reading in from various pipes

  open($in, "-|", "curl", "--silent", $down_url, "-b", $cookiefile)
    || die "Cannot run curl on $down_url";
  my @lines = <$in>;
  close($in) || die "Error running curl on $down_url";
  my $sxml = XML::LibXML->load_xml(string => join("",@lines), recover => 1);
  return "Cannot parse file listing for $down_url"
    unless $sxml;
  my @obj = $sxml->findnodes(q{//folder[@name='IMG Data']/file});
  return "No IMG data" unless @obj;
  # The url attribute is what I want, and it should have a filename that indicates it is the tarball
  my $tar_url = undef;
  foreach my $obj (@obj) {
    foreach my $attr ($obj->attributes) {
      if ($attr->nodeName eq "url" && $attr->nodeValue =~ m/[.]tar[.]gz$/) {
        $tar_url = $attr->nodeValue;
        last;
      }
    }
    last if defined $tar_url;
  }
  return "No IMG URL for a tar.gz file" unless defined $tar_url;
  $tar_url = "https://genome.jgi.doe.gov" . $tar_url;

  my $tmp_tar = "$dir.$$.tar";
  open($in, "-|", "curl", "--silent", $tar_url, "-b", $cookiefile)
    || die "Cannot run curl on $tar_url";
  open(my $tarfh, ">", "$tmp_tar.gz") || die "Cannot write to $tmp_tar.gz";
  while(my $line = <$in>) {
    print $tarfh $line;
  }
  close($in) || return "Error reading from $tar_url";
  close($tarfh) || die "Error writing to $tmp_tar.gz";
  my $tmpdir = "$dir.$$";
  die "temporary directory $tmpdir already exists" if -e $tmpdir;

  unless(ExplodeJGI($tmp_tar, $tmpdir, $genomename, $portalid)) {
    unlink("$tmp_tar.gz");
    unlink($tmp_tar);
    system("rm", "-Rf", $tmpdir);
    return "Invalid tar.gz file from $tar_url";
  }
  unlink($tmp_tar);
  rename($tmpdir, $dir) || die "Error renaming to $tmpdir to $dir";
  return "";
}

sub ExplodeJGI($$) {
  my ($tar, $dir, $name, $portalid) = @_;
  my $in;
  return 0
    unless system("gunzip", "$tar.gz") == 0
      && open($in, "-|", "tar", "tf", $tar);
  my @files = <$in>;
  return 0 unless close($in);
  return 0 unless @files; # empty tarball

  # Check that all files are in the same subdirectory
  foreach my $file (@files) {
    chomp $file;
    $file =~ m!^[^/]+/[^/]+$!
      || return 0;
  }
  my $prefix = $files[0];
  $prefix =~ s!/.*!!;
  mkdir($dir) || die "Cannot create directory $dir";
  # explode and check that prefix.genes.faa and prefix.fna both exist
  return 0
    unless system("tar", "-xf", $tar, "--strip-components=1", "-C", $dir) == 0
      && -e "$dir/$prefix.fna"
      && -e "$dir/$prefix.genes.faa";
  my $out;
  open($out, ">", "$dir/fields") || die "Cannot write to $dir/fields";
  print $out join("\n", $prefix, $name, $portalid)."\n";
  close($out) || die "Error writing to $dir/fields";
  return 1;
}

sub SearchUniProtProteomes {
  my ($query, $limit) = @_;
  return () unless $query;
  $limit = 100 unless defined $limit;
  # Search proteome name, not keywords
  my $URL = "https://www.ebi.ac.uk/proteins/api/proteomes?name=${query}&format=json&size=${limit}";
  my $string = get($URL);
  die "No response from $URL" unless $string;
  my $json = from_json($string);
  die "Not a JSON response from $URL" unless $json;
  die "Not a list from $URL" unless ref($json) eq "ARRAY";
  my @hits = @$json;
  foreach my $hit (@hits) {
    die "No upid in element from $URL" unless $hit->{upid};
    $hit->{gdb} = "UniProt";
    $hit->{gid} = $hit->{upid};
    $hit->{genomeName} = $hit->{name};
    $hit->{URL} = "https://www.uniprot.org/proteomes/" . $hit->{upid};
  }
  return @hits;
}

sub UniProtProteomeInfo {
  my ($gid) = @_;
  die unless $gid;
  my $URL = "https://www.ebi.ac.uk/proteins/api/proteomes?upid=${gid}&format=json";
  my $string = get($URL);
  die "No response from $URL" unless $string;
  my $json = from_json($string);
  die "Not a JSON response from $URL" unless $json;
  die "Not a list from $URL" unless ref($json) eq "ARRAY";
  my ($assembly) = @$json;
  die "No upid in element from $URL" unless $assembly->{upid};
  $assembly->{gdb} = "UniProt";
  $assembly->{gid} = $gid;
  $assembly->{genomeName} = $assembly->{name};
  $assembly->{URL} = "https://www.uniprot.org/proteomes/" . $assembly->{upid};
  return $assembly;
}

1;


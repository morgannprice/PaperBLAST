# Fetch information about assemblies from NCBI

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

our (@ISA,@EXPORT);
@ISA = qw(Exporter);
@EXPORT = qw(FetchAssemblyInfo FetchAssemblyFaa FetchAssemblyFna FetchAssemblyFeatureFile ParseAssemblyFeatureFile
             SearchJGI CreateJGICookie FetchJGI);

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

# Given the assembly hash, fetch the feature table, uncompress it, and write it to the given file
# Writes to a temporary file to avoid overwriting by another process
# Returns success or failure
sub FetchAssemblyFeatureFile($$) {
  my ($assembly, $outfile) = @_;
  return FetchAssemblyFileGz($assembly, "feature_table.txt", $outfile);
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

# From a file to a reference to a list of hashes
sub ParseAssemblyFeatureFile($) {
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
  return [] unless defined $query && $query ne "";
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
        $row->{portalId} =~ m/,"([0-9A-Z_]+)[)]$/
          || return "Cannot find the portal id in " . $row->{"Portal ID"};
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
  print join("\t", "XXX", $username, $passwd)."\n";
  my @cmd = ("curl", "--silent", $URL,
              "--data-urlencode", "login=$username",
              "--data-urlencode", "password=$passwd",
              "-c", $tmpfile);
  print STDERR join("\t", "Command", @cmd)."\n";
  open(my $fh, "-|", @cmd)
    || die "Cannot run curl to $URL";
  my @lines = <$fh>;
  close($fh) || die "Error running curl for $URL";
  die "Cookie file not created" unless -e $tmpfile;
  rename($tmpfile, $cookiefile) || die "renaming cookie file failed";
  return "";
}

# Given a portal id and a genome name, fetch the IMG content into a new directory
# Returns an error message, or ""
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
  # The url attribute is what I want
  my $tar_url = undef;
  my $obj = $obj[0];
  foreach my $attr ($obj->attributes) {
    if ($attr->nodeName eq "url") {
      $tar_url = $attr->nodeValue;
      last;
    }
  }
  return "No IMG URL" unless defined $tar_url;
  my $tar_url = "https://genome.jgi.doe.gov" . $tar_url;
  my $tmp_tar = "$dir.$$.tar.gz";
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

1;


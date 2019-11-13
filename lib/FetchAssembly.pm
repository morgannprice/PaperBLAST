# Fetch information about assemblies from NCBI, JGI, UniProt, MicrobesOnline, Fitness

package FetchAssembly;
require Exporter;
use strict;
use LWP::Simple qw{get}; # for get()
use LWP::UserAgent;
use XML::LibXML;
use Time::HiRes qw{gettimeofday};
use pbutils;
use URI::Escape;
use HTTP::Cookies;
use JSON;
use DBI;
use Digest::MD5;
use CGI qw{end_html a p};

our (@ISA,@EXPORT);
@ISA = qw(Exporter);
@EXPORT = qw(GetMatchingAssemblies CacheAssembly AASeqToAssembly
             warning fail
             SetPrivDir GetPrivDir
             GetFitnessBrowserPath SetFitnessBrowserPath FitnessBrowserDbh
             FetchNCBIInfo FetchNCBIFaa FetchNCBIFna FetchNCBIFeatureFile ParseNCBIFeatureFile
             SearchJGI CreateJGICookie FetchJGI SearchUniProtProteomes UniProtProteomeInfo
             GetMaxNAssemblies
);

# page must be started already; reports any number of errors or warning
sub fail {
  warning(@_);
  print CGI::end_html;
  exit(0);
}

sub warning {
  print p({ -style => "color: red;" }, @_), "\n";
}

my $maxAssemblyList = 100;
sub GetMaxNAssemblies() {
  return $maxAssemblyList;
}

sub FetchWithRetry($$) {
  my ($URL, $nRetry) = @_;
  die if $nRetry < 0;
  for (my $i = 0; $i <= $nRetry; $i++) {
    my $out = get($URL);
    return $out if $out;
    sleep(1);
  }
  return undef;
}

# Note -- arguably should be using an API key for NCBI
# One argument: the query
# Returns a list of results (up to maxAssemblyList) as a hash; each "assembly" includes
# uid (the assembly uid), id (an identifier like GCF_000195755.1), ftp (the URL for the ftp site),
# and org
sub FetchNCBIInfo($) {
  my ($query) = @_;
  # First run the query using esearch
  my $URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=assembly&retmax=${maxAssemblyList}&term=" . $query;
  my $nRetry = 2;
  my $string = FetchWithRetry($URL, $nRetry);
  die "Failed to contact NCBI eutils at $URL\n" unless $string;
  my $fxml = XML::LibXML->load_xml(string => $string, recover => 1);
  my @idnodes = $fxml->findnodes("//IdList/Id");
  my @ids = map { $_->textContent } @idnodes;
  @ids = grep { defined && $_ ne "" } @ids;
  return () if @ids == 0;

  # Now fetch metadata for these ids.
  @ids = $ids[0..($maxAssemblyList-1)] if scalar(@ids) > $maxAssemblyList;
  my $idspec= join(",",@ids);
  $URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=assembly&id=$idspec";
  $string = FetchWithRetry($URL, $nRetry);
  die "Failed to contact NCBI eutils at $URL\n" unless $string;
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
  # Do not know how to search for genomes in IMG only; so,
  # search for 3* too many.
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
                 "rows=" . (3 * $maxAssemblyList));
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
      next unless $row->{"IMG Portal"};
      # Compute the portal.id field
      if (exists $row->{"Portal ID"} && $row->{"Portal ID"} =~ m/"/) {
        $row->{"Portal ID"} =~ m/,"([0-9A-Za-z_]+)"[)]$/
          || return "Cannot find the portal id in " . $row->{"Portal ID"};
        $row->{portalId} = $1;
        $row->{genomeName} = $row->{"Project Name"};
        $row->{URL} = "http://genome.jgi.doe.gov/portal/" . $row->{portalId};
        $row->{gdb} = "IMG";
        $row->{gid} = $row->{portalId};
        push @out, $row;
      }
      last if @out >= $maxAssemblyList;
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
  return "Sorry, there is no IMG data for this assembly." unless @obj;
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
  return "Sorry, there is no IMG tarball for this assembly" unless defined $tar_url;
  $tar_url = "https://genome.jgi.doe.gov" . $tar_url;

  my $tmp_tar = "$dir.$$.tar";
  open($in, "-|", "curl", "--silent", $tar_url, "-b", $cookiefile)
    || die "Cannot run curl on $tar_url";
  open(my $tarfh, ">", "$tmp_tar.gz") || die "Cannot write to $tmp_tar.gz";
  while(my $line = <$in>) {
    print $tarfh $line;
  }
  close($in) || return "Error reading from IMG tarball $tar_url";
  close($tarfh) || die "Error writing to $tmp_tar.gz";
  my $tmpdir = "$dir.$$";
  die "temporary directory $tmpdir already exists" if -e $tmpdir;

  unless(ExplodeJGI($tmp_tar, $tmpdir, $genomename, $portalid)) {
    unlink("$tmp_tar.gz");
    unlink($tmp_tar);
    system("rm", "-Rf", $tmpdir);
    return "Invalid IMG tar.gz file from $tar_url";
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

  # Check that all files are in a subdirectory
  foreach my $file (@files) {
    chomp $file;
    $file =~ m!^[^/]+/[^/]+$!
      || return 0;
  }
  my $prefix = $files[0];
  $prefix =~ s!/.*!!;
  mkdir($dir) || die "Cannot create directory $dir";

  # explode and check that prefix.genes.faa and prefix.fna both exist
  # Note that strip-components will ignore the 1st subdirectory indicator,
  # and above we checked that there are no other subdirectories,
  # and modern tar will not write to / or .. anyway.
  # -C means chdir to $dir before writing.
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
  my ($query) = @_;
  return () unless $query;
  # Search proteome name, not keywords
  my $URL = "https://www.ebi.ac.uk/proteins/api/proteomes?name=${query}&format=json&size=${maxAssemblyList}";
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

# Fetch matching genomes, and return list of hashes
# Each hash includes gid, genomeName, and URL
sub GetMatchingAssemblies($$) {
  my ($gdb,$gquery) = @_;
  return unless $gdb && $gquery;
  if ($gdb eq "NCBI") {
    my @hits = FetchNCBIInfo($gquery);
    foreach my $hit (@hits) {
      $hit->{gdb} = $gdb;
      $hit->{gid} = $hit->{id};
      $hit->{genomeName} = $hit->{org};
      $hit->{URL} = "https://www.ncbi.nlm.nih.gov/assembly/" . $hit->{id};
    }
    return @hits;
  } elsif ($gdb eq "MicrobesOnline") {
    my $mo_dbh = DBI->connect('DBI:mysql:genomics:pub.microbesonline.org', "guest", "guest")
      || fail("Cannot connect to MicrobesOnline:", $DBI::errstr);
    my $hits = $mo_dbh->selectall_arrayref("SELECT taxonomyId, shortName FROM Taxonomy
                                              WHERE shortName LIKE ? OR shortName LIKE ? ORDER BY shortName LIMIT $maxAssemblyList",
                                           { Slice => {} }, $gquery."%", "% ${gquery}");
    foreach my $hit (@$hits) {
      $hit->{gdb} = $gdb;
      $hit->{gid} = $hit->{taxonomyId};
      $hit->{genomeName} = $hit->{shortName};
      $hit->{URL} = "http://www.microbesonline.org/cgi-bin/genomeInfo.cgi?tId=511145";
    }
    return @$hits;
  } elsif ($gdb eq "FitnessBrowser") {
    my $fbdbh = &FitnessBrowserDbh();
    my $orginfo = $fbdbh->selectall_hashref("SELECT * FROM Organism", "orgId");
    my @hits = ();
    foreach my $hash (values %$orginfo) {
      $hash->{genomeName} = join(" ", $hash->{genus}, $hash->{species}, $hash->{strain});
      $hash->{gdb} = $gdb;
      $hash->{gid} = $hash->{orgId};
      $hash->{URL} = "http://fit.genomics.lbl.gov/cgi-bin/org.cgi?orgId=" . $hash->{orgId};
      my $lcg = lc($hash->{genomeName});
      push @hits, $hash
        if substr($lcg, 0, length($gquery)) eq lc($gquery)
          || index($lcg, " " . lc($gquery)) >= 0
          || index($lcg, "-" . lc($gquery)) >= 0;
      last if @hits >= $maxAssemblyList;
    }
    @hits = sort { $a->{genomeName} cmp $b->{genomeName} } @hits;
    return @hits;
  } elsif ($gdb eq "UniProt") {
    return SearchUniProtProteomes($gquery);
  } elsif ($gdb eq "IMG") {
    my ($imgerr, $hits) = SearchJGI($gquery);
    fail("Error searching JGI/IMG: $imgerr") if $imgerr;
    return @$hits;
  }
  # else
  fail("Database $gdb is not supported");
}

my $saved_privdir = "";
sub SetPrivDir($) {
  my ($dir) = @_;
  fail("Not a directory: $dir") unless -d $dir;
  $saved_privdir = $dir;
}

sub GetPrivDir() {
  return $saved_privdir;
}

sub CacheAssembly($$$) {
  my ($gdb, $gid, $dir) = @_;
  return undef unless $gdb && $gid;
  die "Not a directory: $dir\n" unless -d $dir;
  if ($gdb eq "local") {
    die "Invalid gid $gid\n" unless $gid =~ m/^[0-9A-Fa-f]+$/;
    my $faaFile = "$dir/$gid/faa";
    die "No such file: $faaFile\n" unless -e $faaFile;
    # Read the first line to get a more useful genome name
    open (my $fh, "<", $faaFile) || die "Error reading $faaFile";
    my $header = <$fh> || die "$faaFile is empty";
    chomp $header;
    $header =~ s/^>//;
    my $id = $header; $id =~ s/ .*//;
    close($fh) || die "Error reading $faaFile";
    my $assembly = { 'gdb' => $gdb, 'gid' => $gid,
                     'faafile' => $faaFile, 'URL' => "",
                     'genomeName' => "Proteome with $id..." };
    my $fnaFile = "$dir/$gid/fna";
    $assembly->{fnafile} = $fnaFile if -e $fnaFile;
    return $assembly;
  } elsif ($gdb eq "NCBI") {
    my @hits = FetchNCBIInfo($gid);
    fail("Do not recognize NCBI assembly $gid")
      unless @hits;
    my $assembly = $hits[0];
    # note redundancy with code above
    $assembly->{gid} = $assembly->{id};
    $assembly->{genomeName} = $assembly->{org};
    $assembly->{URL} = "https://www.ncbi.nlm.nih.gov/assembly/" . $assembly->{id};

    my $faafile = "$dir/refseq_" . $assembly->{gid} . ".faa";
    my $fnafile = "$dir/refseq_" . $assembly->{gid} . ".fna";
    my $featurefile = "$dir/refseq_" . $assembly->{id} . ".features.tab";
    unless (-e $faafile && -e $fnafile && -e $featurefile) {
      print "<P>Fetching assembly $assembly->{gid}\n";
    }
    unless (-e $faafile) {
      fail("Sorry, failed to fetch the protein fasta file for this assembly ($!). This assembly might not have any predicted proteins.")
        unless FetchNCBIFaa($assembly, $faafile);
    }
    unless (-e $fnafile) {
      fail("Sorry, failed to fetch the nucleotide assembly for this assembly: $!")
        unless FetchNCBIFna($assembly, $fnafile);
    }
    unless (-e $featurefile) {
      fail("Sorry, failed to fetch the feature file for this assembly: $!")
        unless &FetchNCBIFeatureFile($assembly, $featurefile);
    }
    my $features = ParseNCBIFeatureFile($featurefile);
    $assembly->{features} = $features;
    $assembly->{prot} = {};
    $assembly->{oldid} = {};
    foreach my $row (@$features) {
      next unless $row->{class} eq "with_protein";
      $assembly->{prot}{$row->{product_accession}} = $row;
      $assembly->{prot}{$row->{"non-redundant_refseq"}} = $row;
    }
    foreach my $row (@$features) {
      if ($row->{class} eq "protein_coding"
          && $row->{locus_tag}
          && $row->{attributes} =~ m/old_locus_tag=([A-Za-z0-9_]+)/) {
        $assembly->{oldid}{$row->{locus_tag}} = $1;
      }
    }
    $assembly->{faafile} = $faafile;
    $assembly->{fnafile} = $fnafile;
    $assembly->{gdb} = $gdb;
    return $assembly;
  } elsif ($gdb eq "MicrobesOnline") {
    my $mo_dbh = DBI->connect('DBI:mysql:genomics:pub.microbesonline.org', "guest", "guest")
      || fail("Cannot connect to MicrobesOnline:", $DBI::errstr);
    my $taxId = $gid;
    my ($genomeName) = $mo_dbh->selectrow_array(qq{ SELECT shortName FROM Taxonomy WHERE taxonomyId = ? },
                                                 {}, $taxId);
    fail("Unknown taxonomy $taxId") unless defined $genomeName;
    my $assembly = { gdb => $gdb,
                     gid => $gid,
                     genomeName => $genomeName,
                     URL => "http://www.microbesonline.org/cgi-bin/genomeInfo.cgi?tId=$taxId",
                     faafile => "$dir/mogenome_${taxId}.faa",
                     fnafile => "$dir/mogenome_${taxId}.fna" };
    unless (-e $assembly->{faafile}) {
      print p("Loading the proteome of $genomeName from", a({-href => "http://www.microbesonline.org/" }, "MicrobesOnline")), "\n";
      my $desc = $mo_dbh->selectall_hashref(qq{ SELECT locusId, description
                                                    FROM Locus JOIN Scaffold USING (scaffoldId)
                                                    JOIN Description USING (locusId,version)
                                                    WHERE taxonomyId = ? AND isActive=1 AND priority=1 AND Locus.type=1 },
                                            "locusId", {}, $taxId);
      my $sysNames = $mo_dbh->selectall_hashref(qq{ SELECT locusId, name
                                                    FROM Locus JOIN Scaffold USING (scaffoldId)
                                                    JOIN Synonym USING (locusId, version)
                                                    WHERE taxonomyId = ? AND isActive=1 AND priority=1
                                                          AND Locus.type=1 AND Synonym.type = 1 },
                                                "locusId", {}, $taxId);
      my $genes = $mo_dbh->selectall_arrayref(qq{ SELECT locusId, sequence
                                                  FROM Locus JOIN Scaffold USING (scaffoldId)
                                                  JOIN AASeq USING (locusId, version)
                                                  WHERE taxonomyId = ? AND isActive=1 AND priority=1 AND type=1; },
                                              { Slice => {} }, $taxId);
      my $tmpfile = $assembly->{faafile} . ".$$.tmp";
      open(my $fh, ">", $tmpfile) || fail("Cannot write to $tmpfile");
      foreach my $gene (@$genes) {
        my $locusId = $gene->{locusId};
        my $sysName = $sysNames->{$locusId}{name} || "";
        my $desc = $desc->{$locusId}{description} || "";
        print $fh ">$locusId $sysName $desc\n$gene->{sequence}\n";
      }
      close($fh) || fail("Error writing to $tmpfile");
      rename($tmpfile, $assembly->{faafile})
        || fail("Rename $tmpfile to $assembly->{faafile} failed");
    }
    unless (-e $assembly->{fnafile}) {
      print p("Loading the genome of $genomeName from", a({-href => "http://www.microbesonline.org/" }, "MicrobesOnline")), "\n";
      my $sc = $mo_dbh->selectall_arrayref(qq{ SELECT scaffoldId, sequence
                                               FROM Scaffold JOIN ScaffoldSeq USING (scaffoldId)
                                               WHERE taxonomyId = ? AND isActive = 1 },
                                           {}, $taxId);
      fail("Cannot fetch genome sequence for $taxId from MicrobesOnline")
        unless @$sc > 0;
      my $tmpfile = $assembly->{fnafile} . ".$$.tmp";
      open(my $fh, ">", $tmpfile) || fail("Cannot write to $tmpfile");
      foreach my $row (@$sc) {
        my ($scaffoldId, $seq) = @$row;
        print $fh ">${scaffoldId}\n$seq\n";
      }
      close($fh) || fail("Error writing to $tmpfile");
      rename($tmpfile, $assembly->{fnafile})
        || fail("Rename $tmpfile to $assembly->{fnafile} failed");
    }
    return $assembly;
  } elsif ($gdb eq "FitnessBrowser") {
    my $dbh = &FitnessBrowserDbh();
    my $assembly = $dbh->selectrow_hashref("SELECT * FROM Organism WHERE orgId = ?",
                                         {}, $gid);
    fail("Genome $gid in the Fitness Browser is not known")
      unless defined $assembly->{orgId};
    $assembly->{gdb} = $gdb;
    $assembly->{gid} = $gid;
    $assembly->{faafile} = "$dir/fbrowse_${gid}.faa";
    $assembly->{fnafile} = "$dir/fbrowse_${gid}.fna";
    $assembly->{genomeName} = join(" ", $assembly->{genus}, $assembly->{species}, $assembly->{strain});
    $assembly->{URL} = "http://fit.genomics.lbl.gov/cgi-bin/org.cgi?orgId=$gid";
    unless (-e $assembly->{faafile}) {
      my $tmpfile = $assembly->{faafile} . ".$$.tmp";
      open(my $fh, ">", $tmpfile) || fail("Cannot write to $tmpfile");
      my $genes = $dbh->selectall_hashref("SELECT * from Gene WHERE orgId = ?", "locusId", {}, $gid);
      my $aaseqsFile = GetFitnessBrowserPath() . "/aaseqs";
      open(my $aafh, "<", "$aaseqsFile") || die "Cannot read $aaseqsFile";
      my $state = {};
      while (my ($header, $sequence) = ReadFastaEntry($aafh,$state)) {
        my ($org2, $locusId) = split /:/, $header;
        die $header unless defined $locusId;
        next unless $org2 eq $gid;
        my $gene = $genes->{$locusId}
          || die "Unrecognized locusId $locusId";
        print $fh ">$locusId $gene->{sysName} $gene->{desc}\n$sequence\n";
      }
      close($aafh) || die "Error reading $aaseqsFile";
      rename($tmpfile, $assembly->{faafile}) || die "Rename $tmpfile to $assembly->{faafile} failed";
    }
    unless (-e $assembly->{fnafile}) {
      my $tmpfile = $assembly->{fnafile} . ".$$.tmp";
      open(my $fh, ">", $tmpfile) || fail("Cannot write to $tmpfile");
      my $sc = $dbh->selectall_arrayref("SELECT scaffoldId, sequence FROM ScaffoldSeq
                                           WHERE orgId = ?",
                                          {}, $gid);
      foreach my $row (@$sc) {
        my ($scaffoldId,$sequence) = @$row;
        print $fh ">${scaffoldId}\n${sequence}\n";
      }
      close($fh) || fail("Error writing to $tmpfile");
      rename($tmpfile, $assembly->{fnafile}) || die "Rename $tmpfile to $assembly->{fnafile} failed";
    }
    return $assembly;
  } elsif ($gdb eq "UniProt") {
    my $assembly = UniProtProteomeInfo($gid);
    $assembly->{faafile} = "$dir/uniprot_${gid}.faa";
    unless (-e $assembly->{faafile}) {
      # First try getting from https://www.uniprot.org/uniprot/?query=proteome:UP000002886&format=fasta
      # This is a bit slow, so I considered using links like
      # ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Archaea/UP000000242_399549.fasta.gz
      # but those are also surprisingly slow, and would need to figure out which section to look in.
      print p("Loading proteome",
              a({-href => "https://www.uniprot.org/proteomes/$gid"}, $gid),
              "from UniProt"), "\n";
      my $refURL = "https://www.uniprot.org/uniprot/?query=proteome:${gid}&format=fasta";
      my $uniparcURL = "https://www.uniprot.org/uniparc/?query=proteome:${gid}&format=fasta";
      my $faa = get($refURL);
      $faa = get($uniparcURL) unless $faa;
      fail("Proteome $gid seems to be empty, see", a({href => $uniparcURL}, "here"))
        unless $faa;
      my $tmpfile = $assembly->{faafile} . ".$$.tmp";
      open(my $fh, ">", $tmpfile) || fail("Cannot write to $tmpfile");
      print $fh $faa;
      close($fh) || die "Error writing to $tmpfile";
      rename($tmpfile, $assembly->{faafile}) || die "Rename $tmpfile to $assembly->{faafile} failed";
    }
    return $assembly;
  } elsif ($gdb eq "IMG") {
    my $privdir = GetPrivDir();
    fail("PrivDir not set") unless $privdir;
    fail("Invalid project identifier $gid")
      unless $gid =~ m/^[a-zA-Z0-9_]+$/;
    my $subdir = "$dir/$gid";
    unless (-d $subdir) {
      # Query again to find out what the genome is
      my ($err, $hits) = SearchJGI($gid);
      my @hits = grep { $_->{portalId} eq $gid} @$hits;
      fail("Identifier $gid not known in JGI portal") unless @hits > 0;
      my $genomeName = $hits[0]{genomeName};
      print p("Fetching",
              a({ -href => $hits[0]{URL} }, $genomeName),
              "from the JGI portal"), "\n";
      die "No information in PrivDir" unless -e "$privdir/.JGI.info";
      open(my $fh_info, "<", "$privdir/.JGI.info")
        || die "Error creating JGI cookie";
      my $uname = <$fh_info>; chomp $uname;
      my $pwd = <$fh_info>; chomp $pwd;
      close($fh_info) || die "Error creating JGI cookie";
      my $cfile = "$privdir/JGI.$$.cookie";
      $err = CreateJGICookie($uname, $pwd, $cfile);
      fail("Error creating JGI cookie: $err") if $err;
      $err = FetchJGI($gid, $hits[0]{genomeName}, $cfile, $subdir);
      fail($err) if $err;
      unlink($cfile);
    }
    my ($imgId, $genomeName, $fnafile, $faafile);
    open (my $fhfields, "<", "$subdir/fields")
      || fail("Cannot read $subdir/fields");
    my $portalId2;
    ($imgId, $genomeName, $portalId2) = <$fhfields>;
    close($fhfields) || fail("Error reading $subdir/fields");
    chomp $imgId;
    chomp $genomeName;
    chomp $portalId2;
    fail("Invalid $subdir/fields file") unless defined $portalId2 && $portalId2 eq $gid;
    $faafile = "$subdir/$imgId.genes.faa";
    $fnafile = "$subdir/$imgId.fna";
    fail("No such file: $faafile") unless -e $faafile;
    fail("No such file: $fnafile") unless -e $fnafile;
    my $assembly = { gdb => $gdb, gid => $gid, portalId => $gid, genomeName => $genomeName,
                     imgId => $imgId,
                     URL => "http://img.jgi.doe.gov/genome.php?id=" . $imgId,
                     URL2 => "http://genome.jgi.doe.gov/" . $gid,
                     fnafile => $fnafile, faafile => $faafile };
    return $assembly;
  }
  fail("Database $gdb is not supported");
}

# Given a hash of header to protein sequence, compute the CRC, and save it as an assembly
# Note that "local" also supports fna files, but that is not supported
sub AASeqToAssembly($$) {
  my ($aaseq, $dir) = @_;
  my $n = scalar(keys %$aaseq);
  die "Empty input to AASeqToAssembly is not allowed\n" unless $n > 0;
  my $md5 = Digest::MD5->new;
  foreach my $key (sort keys %$aaseq) {
    $md5->add($key, "\n", $aaseq->{$key}, "\n");
  }
  my $hex = $md5->hexdigest;
  -d "$dir/$hex" || mkdir("$dir/$hex") || die "Cannot make directory $dir/$hex";
  my $faaFile = "$dir/$hex/faa";
  unless (-e $faaFile) {
    my $tmpFile = $faaFile . ".$$.tmp";
    open(my $fhA, ">", $tmpFile) || die "Cannot write to $tmpFile";
    foreach my $key (sort keys %$aaseq) {
      print $fhA ">" . $key . "\n" . $aaseq->{$key} . "\n";
    }
    close($fhA) || die "Error writing to $tmpFile";
    rename($tmpFile, $faaFile) || die "Cannot rename $tmpFile to $faaFile";
  }
  return CacheAssembly('local', $hex, $dir);
}

my $fbdb_path = "";
sub SetFitnessBrowserPath($) {
  my ($path) = @_;
  fail("No such directory: $path") unless -d $path;
  $fbdb_path = $path;
}

sub GetFitnessBrowserPath() {
  return $fbdb_path;
}

my $fbdbh; # set just once
sub FitnessBrowserDbh() {
  return $fbdbh if $fbdbh;
  my $path = GetFitnessBrowserPath();
  fail("Fitness Browser path is not set") unless $path;
  my $dbfile = "$path/feba.db";
  fail("Fitness Browser database $dbfile is missing") unless -e $dbfile;
  $fbdbh = DBI->connect("dbi:SQLite:dbname=$dbfile","","",{ RaiseError => 1 })
    || fail("Cannot connect to Fitness Browser database: " . $DBI::errstr);
  return $fbdbh;
}


1;


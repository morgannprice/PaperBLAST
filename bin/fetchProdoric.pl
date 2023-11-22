#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use LWP::UserAgent;
use LWP::Simple qw(get);
use JSON;
use FindBin qw{$RealBin};
use lib "$RealBin/../lib";
use pbutils qw{ReadTable};

my $usage = <<END
loadProdoric.pl -out directory > curated_parsed

Fetches all protein information, references, and weight matrices from
prodoric using the API (https://www.prodoric.de/api). Writes a
tab-delimited curated-parsed file to stdout. Also writes the weight
matrices in the output directory as tab-delimited files with names
like MX000160.0.dat as well as a table noseq.tsv describing matrices
with no uniprotIds.

Optional arguments:
-uniprot uniprotIds.tsv -- uniprotId assignments for matrices that are
   missing them. Must be tab-delimited with fields organism, strain,
   name (these should match prodoric) and uniprotIds, subunits.
   For multi-subunit TFs, uniprotIds and subunits are comma-delimited and
   subunits names each component. Otherwise, subunits should be empty.
-matrix 1 2 -- specific a list of matrices to fetch (useful for testing).
END
;

my @matrixSpec;
my $outDir;
my $uniprotTsv;
die $usage
  unless GetOptions('out=s' => \$outDir,
                    'uniprot=s' => \$uniprotTsv,
                    'matrix=s{1,}' => \@matrixSpec)
  && @ARGV == 0
  && defined $outDir;
die "Not a directory: $outDir\n" unless -d $outDir;

my %uniprot; # organism::strain::name => row
if (defined $uniprotTsv) {
  my @uniprots = ReadTable($uniprotTsv, ["organism","strain","name","uniprotIds","subunits"]);
  foreach my $row (@uniprots) {
    my $key = join("::", $row->{organism}, $row->{strain}, $row->{name});
    die "Duplicate row for $key in $uniprotTsv\n" if exists $uniprot{$key};
    $uniprot{$key} = $row;
  }
}

my $ua = LWP::UserAgent->new(ssl_opts => { verify_hostname => 0 });

if (@matrixSpec == 0) {
  my $url = "https://www.prodoric.de/api/matrix";
  my $response = $ua->get($url);
  die "Cannot access $url: " . $response->status_line . "\n"
    unless $response->is_success;
  my $list = decode_json($response->decoded_content);
  foreach my $hash (@$list) {
    die "No mx entry" unless exists $hash->{mx};
    push @matrixSpec, $hash->{mx};
  }
}

print STDERR "Fetching " . scalar(@matrixSpec) . " matrices from PRODORIC\nWriting to $outDir\n";

open(my $fhNoseq, ">", "$outDir/noseq.tsv")
  || die "Cannot write to $outDir/noseq.tsv";
# Fill out the 1st three fields using data from PRODORIC
print $fhNoseq join("\t", qw{organism strain name uniprotIds subunits})."\n";

my $nNoUniprot = 0;
my $nDeprecated = 0;
foreach my $matrix (@matrixSpec) {
  die "Invalid matrix: $matrix\n" unless $matrix =~ m/^\d+$/;
  my $url = "https://www.prodoric.de/api/matrix/MX".$matrix;
  my $response = $ua->get($url);
  die "Cannot access $url: " . $response->status_line . "\n"
    unless $response->is_success;
  my $hash = decode_json($response->decoded_content);
  if ($hash->{deprecated} eq "true") {
    $nDeprecated++;
    next;
  }
  foreach my $key (qw{mx type name organism pwm xrefs}) {
    die "No value for $key in $url\n" unless exists $hash->{$key};
  }
  my @uniprotIds;
  my @subunits;
  foreach my $xref (@{ $hash->{xrefs} }) {
    $uniprotIds[0] = $xref->{acc} if $xref->{name} eq "uniprot";
  }
  if (@uniprotIds == 0) {
    # Try to find it using the %uniprot table
    my @fields = ($hash->{organism}{name}, $hash->{organism}{strain}, $hash->{name});
    my $key = join("::", @fields);
    if (exists $uniprot{$key}) {
      my $row = $uniprot{$key};
      @uniprotIds = split /,/, $row->{uniprotIds};
      foreach my $uniprotId (@uniprotIds) {
        die "Invalid uniprotId $uniprotId for $key from $uniprotTsv\n"
          unless $uniprotId =~ m/[A-Z][A-Z0-9_]+$/;
      }
      if (@uniprotIds > 1) {
        @subunits = split /,/, $row->{subunits};
        die "Wrong number of subunits for $key in $uniprotTsv\n"
          unless scalar(@uniprotIds) == scalar(@subunits);
      }
    } else { # no match for $key
      $nNoUniprot++;
      print $fhNoseq join("\t", @fields, "", "")."\n";
      next;
    }
  }
  my @pmIds;
  foreach my $ref (@{ $hash->{references} }) {
    if (exists $ref->{xrefs}) {
      foreach my $xref (@{ $ref->{xrefs} }) {
        push @pmIds, $xref->{acc} if $xref->{name} eq "pubmed";
      }
    }
  }
  die unless @uniprotIds > 0;
  foreach my $i (0..(scalar(@uniprotIds)-1)) {
    my $uniprotId = $uniprotIds[$i];
    die "Invalid uniprotId $uniprotId for $matrix"
      unless $uniprotId =~ m/^[A-Z][A-Z0-9_]+$/;
    my $fasta = get("https://rest.uniprot.org/uniprotkb/$uniprotId.fasta");
    die "Invalid uniprotId $uniprotId\n" unless defined $fasta && $fasta =~ m/^>/;
    my @fastaLines = split /\n/, $fasta;
    shift @fastaLines;
    my $seq = join("", @fastaLines);
    my $desc;
    if (lc($hash->{type}) eq "transcription factor") {
      $desc = "Transcription factor " . $hash->{name};
    } elsif (lc($hash->{type}) eq "response regulator") {
      $desc = "Response regulator " . $hash->{name};
    } elsif (lc($hash->{type}) eq "sigma factor") {
      $desc = "Sigma factor " . $hash->{name};
    } else {
      print STDERR "Unknown type $hash->{type} for matrix $matrix\n";
      $desc = $hash->{name};
    }
    $desc .= ", $subunits[$i] subunit" if $subunits[$i];
    if ($hash->{regulation} eq "+") {
      $desc .= " (activator)";
    } elsif ($hash->{regulation} eq "-") {
      $desc .= " (repressor)";
    } elsif ($hash->{regulation} eq "+-") {
      $desc .= " (activator/repressor)";
    }
    my $id = sprintf("MX%06d.$i", $matrix, $i);
    print join("\t", "prodoric", $id, $uniprotIds[$i],
               $hash->{name},
               $desc,
               join(" ", $hash->{organism}{name}, $hash->{organism}{strain}),
               $seq,
               "", # comment (should I summarize the known sites?)
               join(",", @pmIds))."\n";
    # And save the weight matrix, in transfac format
    if (exists $hash->{pwm}) {
      my $nA = $hash->{pwm}{A};
      my $nC = $hash->{pwm}{C};
      my $nG = $hash->{pwm}{G};
      my $nT = $hash->{pwm}{T};
      my $n = scalar(@$nA);
      die "Empty pwm for $id" unless $n > 0;
      die "Inconsistent lengths for pwm"
        unless scalar(@$nC) == $n
          && scalar(@$nG) == $n
          && scalar(@$nT) == $n;
      my $datFile = "$outDir/$id.dat";
      open(my $fhWM, ">", $datFile) || die "Cannot write to $datFile";
      print $fhWM join("\t", qw{A C G T})."\n";
      for (my $i = 0; $i < $n; $i++) {
        print $fhWM join("\t", $nA->[$i], $nC->[$i], $nG->[$i], $nT->[$i])."\n";
      }
      close($fhWM) || die "Error writing to $datFile\n";
    }
  } # end loop over subunits
} # end loop over matrices
close($fhNoseq) || die "Error writing to $outDir/noseq.tsv\n";

print STDERR "Skipped $nNoUniprot entries with no uniprot and $nDeprecated deprecated entries\n";

#!/usr/bin/perl -w
# Compare a genome to an HMM query using HMMer

# Required CGI parameters:
# gdb -- one of the genome sources supported by FetchAssembly
# gid -- which organism or assembly in that database
# hmmId -- what hmm to search for
#
# Optional CGI parameters:
# byeval -- ignore trusted cutoff and use (per-domain) E < 0.01 instead

use strict;
use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use Time::HiRes qw{gettimeofday};
use IO::Handle; # for autoflush
use URI::Escape; # for uri_escape()
use lib "../lib";
use pbutils;
use FetchAssembly;
use pbweb;

my $nCPU = 4;
my $hmmsearch = "../bin/hmmsearch";
die "No such executable: $hmmsearch" unless -x $hmmsearch;
my $cgi = CGI->new;
my $gdb = $cgi->param('gdb') || die "must specify gdb";
my $gid = $cgi->param('gid');
die "must specify gid" unless defined $gid && $gid ne "";

my $hmmId = $cgi->param('hmmId');
die "Must specify hmmId" unless $hmmId;
$hmmId =~ s/ +$// if defined $hmmId;
$hmmId =~ s/^ +// if defined $hmmId;
my $hmmFile = HmmToFile($hmmId);
if (!defined $hmmFile) {
  my $id2 = HmmNameToId($hmmId);
  if (defined $id2) {
    $hmmId = $id2;
    $hmmFile = HmmToFile($hmmId);
  }
}

# To allow access to fitness browser genomes
my $fbdata = "../fbrowse_data"; # path relative to the cgi directory
SetFitnessBrowserPath($fbdata);
SetPrivDir("../private");

my $tmpDir = "/tmp";
my $procId = $$;
my $timestamp = int (gettimeofday() * 1000);
my $baseFile = $tmpDir . "/genomeHmm." . $procId . $timestamp;

my $assembly = CacheAssembly($gdb, $gid, "../tmp/downloaded")
   || fail("Cannot fetch assembly $gid from database $gdb");

&start_page('title' => "Genome vs. $hmmId");
autoflush STDOUT 1; # show preliminary results

if (!defined $hmmFile) {
  print p("Sorry, HMM", HTML::Entities::encode($hmmId), " is not known");
  exit(0);
}
my $hmmName = `egrep '^NAME' $hmmFile`;
$hmmName =~ s/[\r\n].*//;
$hmmName =~ s/^NAME +//;
my $showId = $hmmId;

print h3("Search for", $hmmId,
         "(".HTML::Entities::encode($hmmName).")",
         "in",
         HTML::Entities::encode($assembly->{genomeName})),
  GetMotd();

my $faaFile = $assembly->{faafile} || fail("Sorry, this assembly has no protein predictions");
open (my $fhFaa, "<", $faaFile) || die "Cannot read $faaFile\n";
my %seqs = ();
my %descs = ();
my $state = {};
while (my ($header, $seq) = ReadFastaEntry($fhFaa, $state)) {
  my $desc = $1 if $header =~ m/ (.*)$/;
  $header =~ s/ .*//;
  $seqs{$header} = $seq;
  $descs{$header} = $desc if defined $desc && $desc ne "";
}
close ($fhFaa) || die "Error reading $faaFile\n";

my $byE = $cgi->param('byeval') ? 1 : 0; # ignore trusted cutoff if set
my @tc;
open(my $fh, "<", $hmmFile) || die "Cannot read $hmmFile";
while(my $line = <$fh>) {
  my @parts = split /\s+/, $line;
  if (@parts == 3 && $parts[0] eq "TC") {
    shift @parts;
    @tc = @parts;
    last;
  }
}
close($fh) || die "Error reading $hmmFile";
$byE = 1 if @tc == 0; # Always use by E if no trusted cutoff

my $resultsFile = "$baseFile.hits";
my @cmd;
if ($byE) {
  print p("Searching for hits with E &le; 0.01");
  @cmd = ($hmmsearch, "-E", 0.01, "-o", "/dev/null", "--domtblout", $resultsFile, $hmmFile, $faaFile);
} else {
  print p("Searching for hits above the trusted cutoff");
  @cmd = ($hmmsearch, "--cut_tc", "-o", "/dev/null", "--domtblout", $resultsFile, $hmmFile, $faaFile);
}
print "\n";
system(@cmd) == 0 || die "hmmsearch failed: $!";

# Parse the results. For each unique sequence id, a hash
# of hmmAcc, hmmlen, evalue, score, and hits, which is a list of hashes,
# each containing domscore, hmmfrom, hmmto, envfrom, envto [envelope coordinates]
my %hits = ();

open(my $fhIn, "<", $resultsFile) || die "Cannot read $resultsFile";
while(my $line = <$fhIn>) {
  next if $line =~ m/^#/;
  my ($protId, undef, $seqlen, undef, $hmmAcc, $hmmlen,
      $evalue, $score, undef,
      undef, undef, undef, undef, $domscore, undef, $hmmfrom, $hmmto, undef, undef,
      $envfrom, $envto) = split / +/, $line;
  die "Cannot parse line $line"
    unless defined $envto
      && $seqlen =~ m/^\d+$/ && $hmmlen =~ m/^\d+$/
      && $hmmfrom =~ m/^\d+$/ && $hmmto =~ m/^\d+$/
      && $envfrom =~ m/^\d+$/ && $envto =~ m/^\d+$/;
  $hits{$protId} = { evalue => $evalue, score => $score, hmmAcc => $hmmAcc, hmmlen => $hmmlen, hits => [] }
    unless exists $hits{$protId};
  push @{ $hits{$protId}{hits} },
    { domscore => $domscore, hmmfrom => $hmmfrom, hmmto => $hmmto, envfrom => $envfrom, envto => $envto };
}
close($fhIn) || die "Error reading $resultsFile";
unlink($resultsFile);

my $URL = pbweb::HMMToURL($hmmId);
my $hmmLink = a({-href => $URL}, $hmmId);
if (scalar(keys %hits) > 0) {
  my $evalueSection = $byE ? "with E &le; 0.01" :
    "above the trusted cutoff of " . escapeHTML(join(",",@tc)) . " bits";
  print p($hmmLink, "hits", scalar(%hits), "sequences in the genome", $evalueSection); 
} else {
  print p("Sorry, no hits for", $hmmLink);
}

my @ids = sort { $hits{$b}{score} <=> $hits{$a}{score} } (keys %hits);
foreach my $id (@ids) {
  my $obj = $hits{$id};
  my $desc = "";
  $desc = b($descs{$id}) if exists $descs{$id};
  my @alnLines = ();
  my $seqlen = length($seqs{$id});
  foreach my $hit (@{ $hits{$id}{hits} }) {
    push @alnLines, "$hit->{envfrom}:$hit->{envto}/$seqlen matches $hit->{hmmfrom}:$hit->{hmmto}/$obj->{hmmlen} of $hmmId";
  }
  my $seqLink = "litSearch.cgi?query=".uri_escape(">id\n". $seqs{$id});
  print p(a({-href => $seqLink}, $id), $desc,
          br(),
          "$obj->{score} bits, E = $obj->{evalue}",
          br(),
          join(br(), @alnLines));
}

print p("Try another search:",
        start_form(-name => 'input', -method => 'GET', -action => 'genomeHmm.cgi'),
        qq{<input type="hidden" name="gdb" value="$gdb">},
        qq{<input type="hidden" name="gid" value="$gid">},
        "HMM",
        textfield(-name => "hmmId", value => '', -size => 12),
        submit('Search'),
        end_form,
        br(),
        "or",
        a({-href => "genomeSearch.cgi?gdb=${gdb}&gid=${gid}"}, "Curated BLAST"),
        "or",
        a({-href => "genomeBlast.cgi?gdb=${gdb}&gid=${gid}"}, "BLAST this genome"));

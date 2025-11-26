#!/usr/bin/perl -w
# Compare a genome to a protein query using usearch

# Required CGI parameters:
# gdb -- one of the genome sources supported by FetchAssembly
# gid -- which organism or assembly in that database
#
# Optional CGI parameters:
# query -- a protein query, which can be specified a variety of ways
#    (as in the main PaperBLAST page)

use strict;
use CGI qw(:standard Vars start_ul end_ul);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use Time::HiRes qw{gettimeofday};
use URI::Escape; # for uri_escape()
use HTML::Entities;
use IO::Handle; # for autoflush
use lib "../lib";
use pbutils qw{ReadFastaEntry};
use FetchAssembly; # for CacheAssembly()
use pbweb;

my $maxEvalue = 0.001;

my $usearch = "../bin/usearch";
die "No such executable: $usearch" unless -x $usearch;

my $nCPU = 4;
my $base = "../data";
my $sqldb = "$base/litsearch.db";
my $dbh = DBI->connect("dbi:SQLite:dbname=$sqldb","","",{ RaiseError => 1 }) || die $DBI::errstr;

my $cgi = CGI->new;
my $gdb = $cgi->param('gdb');
die "Invalid genome database $gdb\n" unless defined $gdb && $gdb =~ m/^[a-zA-Z0-9]+$/;
my $gid = $cgi->param('gid');
die "Invalid genome identifier $gid\n" unless defined $gid && $gid =~ m/^[0-9A-Za-z_:.-]+$/;
my $query = $cgi->param('query') || "";

&start_page('title' => "Genome BLAST",
           'bannerURL' => "genomeBlast.cgi");
autoflush STDOUT 1; # show preliminary results

# To allow access to fitness browser genomes
my $fbdata = "../fbrowse_data"; # path relative to the cgi directory
SetFitnessBrowserPath($fbdata);
SetPrivDir("../private");

my $tmpDir = "../tmp/tmp";
my $procId = $$;
my $timestamp = int (gettimeofday() * 1000);
my $basefile = $tmpDir . "/" . $procId . $timestamp;

my $assembly = CacheAssembly($gdb, $gid, "../tmp/downloaded")
   || fail("Cannot fetch assembly $gid from database $gdb");

my ($def, $seq) = parseSequenceQuery(-query => $query,
                                     -dbh => $dbh,
                                     -blastdb => "$base/uniq.faa",
                                     -fbdata => $fbdata);
my $hasDef = defined $def && $def ne "";
if ($seq) {
  $def = sequenceToHeader($seq) if ! $hasDef;
  $query = ">$def\n$seq\n";
}

if (!$ seq) {
  print join("\n",
    h3("Compare a protein to the predicted proteome of",
       HTML::Entities::encode($assembly->{genomeName})),
    GetMotd(),
    start_form( -name => 'input', -method => 'GET', -action => 'genomeBlast.cgi'),
    qq{<input type="hidden" name="gdb" value="$gdb">},
    qq{<input type="hidden" name="gid" value="$gid">},
    p(b("Enter a protein sequence in FASTA or Uniprot format,",
        br(),
        "or an identifier from UniProt, RefSeq, or MicrobesOnline:"),
      br(),
      textarea( -name  => 'query', -value => '', -cols  => 70, -rows  => 10 ),
      p(submit('Blast'), reset())),
    end_form), "\n";
  finish_page();
} else {
  print h3("Comparing", HTML::Entities::encode($def),
           "to", br(),
           HTML::Entities::encode($assembly->{genomeName})), "\n";
  my $faaFile = $assembly->{faafile} || fail("Sorry, this assembly has no protein predictions");
  open (my $fhFaa, "<", $faaFile) || die "Cannot read $faaFile\n";
  my %seqs = ();
  my $state = {};
  while (my ($header, $s) = ReadFastaEntry($fhFaa, $state)) {
    $header =~ s/ .*//;
    $seqs{$header} = $s;
  }
  close ($fhFaa) || die "Error reading $faaFile\n";
  my $queryFile = "$basefile.faa";
  open(my $fhQuery, ">", $queryFile);
  print $fhQuery ">$def\n" . $seq . "\n";
  close($fhQuery) || die "Error writing to $basefile.faa\n";
  my $queryLen = length($seq);

  print join("\n",
             "<SCRIPT>",
             qq(ShowQuery = function() {
                    document.getElementById("showSequenceLink").style.display = "none";
                    document.getElementById("querySequence").style.display = "block";
                    return false;
                 }),
             "</SCRIPT>"), "\n";
  my @pieces = $seq =~ /.{1,60}/g;
  print p({ -id => 'showSequenceLink', -style => "font-size:90%;" },
          a({ -href => "#", -onclick => "return ShowQuery()" },
            "Show query sequence"));
  print p({ -id => 'querySequence',
            -style => "font-family: monospace; display:none; font-size:90%; padding-left: 2em;"},
          join(br(), ">" . HTML::Entities::encode($def), @pieces));
  print "\n";

  print p("Running ublast with E < $maxEvalue"), "\n";
  my $hitsFile = "$basefile.hits";
  my $cmd = "$usearch -ublast $queryFile -db $faaFile -evalue $maxEvalue --threads $nCPU --quiet -blast6out $hitsFile";
  system($cmd) == 0 || die "ublast failed: $cmd\n$!";
  open(my $fhHits, "<", $hitsFile) || die "Cannot read $hitsFile\n";
  my @hits = ();
  while (my $line = <$fhHits>) {
    chomp $line;
    my @F = split /\t/, $line;
    push @hits, \@F;
  }
  close($fhHits) || die "Error reading $hitsFile\n";
  unlink($queryFile);
  unlink($hitsFile);

  print p("Found", scalar(@hits), @hits == 1 ? "hit" : "hits"), "\n";
  if (@hits > 0) {
    my @header = ("Protein",
                  a({-title => "Percent identity", -style => "color:black;"}, '%Id'),
                  a({-title => "Percent coverage of query", -style => "color:black;"}, '%Cov'),
                  "Description",
                  a({-title => "Expected number of hits this strong or stronger"}, "E"));
    my @th = map th($_), @header;
    print "<TABLE cellpadding=2 cellspacing=0 border=1 >", Tr({align => "left"}, @th), "\n";
    foreach my $hit (@hits) {
      my (undef, $subject, $identity, $alen, $mm, $gap, $qbeg, $qend, $sbeg, $send, $evalue, $bits) = @$hit;
      my $id = $subject;
      $id =~ s/ .*//;
      my $desc = $subject;
      $desc =~ s/^\S+ +//;
      my $cov = int(0.5 + 100 * ($qend-$qbeg+1)/$queryLen);
      my $idShow = $id;
      my $hitSeq = $seqs{$id} || ""; # should always exist
      my $hitLen = length($hitSeq);
      $idShow = a({-title => "PaperBLAST",
                   -href => "litSearch.cgi?query=".uri_escape(">$id\n" . $seqs{$id})},
                  $id);
      my @td = (td({align => "left"}, $idShow),
                td({align => "right"},
                   a({-title => "$mm mismatches and $gap gaps out of $alen"}, $identity)),
                td({align => "right"},
                   a({-title => "$qbeg:$qend/$queryLen of query aligns to $sbeg:$send/$hitLen of $id"}, $cov)),
                td({align => "left"}, HTML::Entities::encode($desc)),
                td({align => "right"},
                   a({-title => "$bits bits"}, $evalue)));
      print Tr(@td), "\n";
    }
    print "</TABLE>\n";
  }
  print p("Or try",
          a({-href => "litSearch.cgi?query=" . uri_escape($query)}, "PaperBLAST"),
          "on the query protein or",
          a({-href => "genomeSearch.cgi?gdb=$gdb&gid=$gid"}, "Curated BLAST for Genomes"),
         "on the genome."),
        "\n";
  finish_page();
}

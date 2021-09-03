#!/usr/bin/perl -w
# psort.cgi -- given a single sequence, run psortb and show the results
#
# Required CGI parameters:
# seq -- the sequence
# name -- a name or description of the query (such as the definition line from a fasta file)
# type -- negative (for Gramm-negative), positive, or archaea

use strict;
use CGI qw{:standard Vars};
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use HTML::Entities;
use Time::HiRes qw{gettimeofday};
use FindBin qw{$RealBin};
use IO::Handle; # for autoflush
use lib "../lib";
use pbweb qw{finish_page};

my $psortb = "$RealBin/../bin/psortb";
die "No such executable: $psortb" unless -x $psortb;

my $cgi=CGI->new;
my $seq = $cgi->param('seq') || "";
$seq =~ s/[ *\r\n]//g;
die "No sequence specified" unless $seq =~ m/^[A-Za-z]+$/;
die "Maximum sequence length is 100K" if length($seq) > 100*1000;

my $name = $cgi->param('name');
die "No name specified" unless defined $name && $name ne "";
$name =~ s/[\r\n]+$//;
die "Invalid name" if $name =~ m/[\r\n]/;

my $type = $cgi->param('type');
die "No type specified" unless $type;
die "Unknown type" unless $type eq "negative" || $type eq "positive" || $type eq "archaea";

my $typeShow = $type;
$typeShow = "Gram-" . $type unless $type eq "archaea";
my $title = "PSORTb v3.0 ($typeShow)";
print
  header(-charset => 'utf-8'),
  start_html(-title => $title),
  h2($title),
  p("Running",
          a({-href => "https://pubmed.ncbi.nlm.nih.gov/20472543/" },
                  "PSORTb v3.0"),
          "on",
          HTML::Entities::encode($name),
          "(" . length($seq) . " amino acids)"),
  "\n";
autoflush STDOUT 1; # show preliminary results

my $procId = $$;
my $timestamp = int (gettimeofday() * 1000);
my $filename = $procId . $timestamp;
my $seqFile = "/tmp/psortb.$filename.fasta";
open(my $fh, ">", $seqFile) || die "Cannot write to $seqFile";
print $fh (">", HTML::Entities::encode($name), "\n", $seq, "\n");
close($fh) || die "Error writing to $seqFile";

print "<PRE>\n";
system($psortb, "--" . $type, $seqFile);
print "</PRE>\n";

unlink($seqFile);
my $newline = "%0A";

my @pieces = $seq =~ /.{1,60}/g;

print
  h3("Other sequence analysis tools"),
  p(a({-href => "http://papers.genomics.lbl.gov/cgi-bin/litSearch.cgi?query=>${name}$newline$seq"},
      "PaperBLAST"),
    "(search for papers about homologs of this protein)"),
  p(a({-href => "http://www.ncbi.nlm.nih.gov/Structure/cdd/wrpsb.cgi?seqinput=>${name}$newline$seq"},
      "Search CDD"),
    "(the Conserved Domains Database, which includes COG and superfam)"),

  # See documentation of HMMer web server (version 1.0) at
  # https://hmmer-web-docs.readthedocs.io/en/latest/searches.html
  # Set E= and domE= to cause weak hits to appear (they are still labeled insignificant)
  p(start_form(-name => "PfamForm", -id => "PfamForm",
               -method => "POST", -action => "https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan"),
    hidden('seq', ">$name\n$seq"),
    hidden('hmmdb', 'pfam'),
    hidden('E', '1'),
    hidden('domE', '1'),
    submit(-style => "display: none;", -id => 'PfamButton'),
    end_form,
    a({-href => "javascript:document.getElementById('PfamForm').submit()",
       -title => "Pfam search at EBI (must be run in this window)"},
       "Search PFam"),
    "(including for weak hits, up to E = 1)"),

  p(a({ -href => "http://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/FindSequence.pl?pasted=$seq",
        -title => "Find similar proteins with known structures (PDBsum)"},
      "Search structures")),

  p("Predict transmembrane helices:",
    a({-href => "http://www.cbs.dtu.dk/cgi-bin/webface2.fcgi?"
             . join("&",
                    "configfile=/usr/opt/www/pub/CBS/services/TMHMM-2.0/TMHMM2.cf",
                    "outform=-noshort",
                    "SEQ=>${name}$newline$seq")},
        "TMHMM")),

  p("Check the current SEED with",
    # %0A encodes "\n" so that it looks like fasta input.
    a({-href => "http://pubseed.theseed.org/FIG/seedviewer.cgi?page=FigFamViewer&fasta_seq=>${name}%0A$seq"},
      "FIGfam search")),
  p("Find homologs in the",
    a({-href => "https://iseq.lbl.gov/genomes/seqsearch?sequence=>${name}%0A$seq"},
      "ENIGMA genome browser")),

  h3("Input"),
  join("\n",
       "<PRE>",
       ">" . HTML::Entities::encode($name),
       @pieces,
       "</PRE>"),
  "\n";
finish_page();

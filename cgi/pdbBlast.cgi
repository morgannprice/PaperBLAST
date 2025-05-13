#!/usr/bin/perl -w

# Optional CGI arguments:
# query -- the query sequence in FASTA or UniProt or plain format, or, an identifier (as in litSearch.cgi)


use strict;
use strict;
use CGI qw(:standard Vars start_ul);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use Time::HiRes qw{gettimeofday};
use DBI;
use LWP::Simple qw{get};
use HTML::Entities;
use IO::Handle; # for autoflush
use lib "../lib";
use pbweb;
use URI::Escape;

my $cgi=CGI->new;
my $query = $cgi->param('query') || "";

my $base = "../data";
my $blastdb = "$base/uniq.faa";
my $sqldb = "$base/litsearch.db";
my $dbh = DBI->connect("dbi:SQLite:dbname=$sqldb","","",{ RaiseError => 1 }) || die $DBI::errstr;
die "No such file: $blastdb.pin\n" unless -e "$blastdb.pin";
my $fbdata = "../fbrowse_data"; # path relative to the cgi directory
my $tmpDir = "../tmp";

my $pdbDb = "$base/pdbClust.faa";
my $blastall = "../bin/blast/blastall";
die "No such executable: $blastall\n" unless -x $blastall;

autoflush STDOUT 1; # show preliminary results
my ($def, $seq) = parseSequenceQuery(-query => $query,
                                     -dbh => $dbh,
                                     -blastdb => $blastdb,
                                     -fbdata => $fbdata);
my $hasDef = defined $def && $def ne "";
if ($seq) {
  $def = sequenceToHeader($seq) if ! $hasDef;
  $query = ">$def\n$seq\n";
}

start_page('title' => 'BLAST against PDB');
print GetMotd(), "\n";
if (!defined $seq) {
  print
    start_form( -name => 'input', -method => 'GET', -action => 'pdbBlast.cgi'),
    p(br(),
      b("Enter a protein sequence in FASTA or Uniprot format,<BR>or an identifier from UniProt, RefSeq, or MicrobesOnline: "),
      br(),
      textarea( -name  => 'query', -value => '', -cols  => 70, -rows  => 10 )),
    p(submit('Search'), reset()),
    end_form;
} else {
  die "No sequence to search" unless $seq;
  my $initial = substr($seq, 0, 10);
  my $seqlen = length($seq);
  $initial .= "..." if $seqlen > 10;
  $initial = "($seqlen a.a., $initial)";
  $initial = "" if ! $hasDef;
  print
    h3("Comparing PDB to", HTML::Entities::encode($def), $initial),
    "\n";

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

  my $procId = $$;
  my $timestamp = int (gettimeofday() * 1000);
  my $filename = $procId . $timestamp;
  my $seqFile = "$tmpDir/$filename.fasta";
  my $hitsFile = "$tmpDir/$filename.pdbhits";
  open(my $fh, ">", $seqFile) || die "Cannot write to $seqFile";
  print $fh ">query\n", $seq, "\n";
  close($fh) || die "Error writing to $seqFile";
  my $cmd = qq{$blastall -p blastp -i $seqFile -d $pdbDb -e 1e-3 -F "m S" -m 8 -o $hitsFile > /dev/null 2>&1};
  system($cmd) ==0 || die "blast failed: $cmd\n$!";
  open($fh, "<", $hitsFile) || die "Cannot read $hitsFile";
  my @hits = ();
  while  (my $line = <$fh>) {
    my @F = split /\t/, $line;
    push @hits, \@F;
  }
  close($fh) || die "Error reading $hitsFile";
  unlink($seqFile);
  unlink($hitsFile);
  my $dbDesc = "non-redundant PDB (which was clustered at 80% identity)";
  if (@hits == 0) {
    print p("Sorry, no significant homologs were found against  using BLAST against $dbDesc, with E &le; 1e-3");
  } else {
    print p(scalar(@hits), "hits were found against $dbDesc");
    foreach my $hit (@hits) {
      my (undef, $subject, $identity, $alen, $mm, $gap, $qbeg, $qend, $sbeg, $send, $evalue, $bits) = @$hit;
      $subject =~ m/^([0-9a-zA-Z]+)_([0-9a-zA-Z]+)$/ || die "invalid subject $subject";
      my ($id, $chain) = ($1,$2);
      my ($desc) = $dbh->selectrow_array("SELECT desc FROM PdbClustInfo WHERE id = ? AND chain = ?",
                                         {}, $id, $chain);
      $desc = "" if !defined $desc;
      my $cov = ($qend-$qbeg+1)/length($seq);
      print p(a({-href => "https://rcsb.org/structure/$id",
                 -style => "text-decoration: none;",
                 -title => "see structure at RCSB"}, $id),
              a({-href => "https://rcsb.org/sequence/$id#$chain",
                 -style => "text-decoration: none;",
                 -title => "see sequence at RCSB"}, $chain),
              b(HTML::Entities::encode($desc)),
              br(),
              small(sprintf(qq{<B>%0.f%%</B> identity, <A title="$qbeg:$qend of query aligns to $sbeg:$send of $id:$chain">}
                            . qq{<B>%0.f%%</B> coverage of query</A>, %.1f bits, E = %g},
                            $identity, $cov*100, $bits, $evalue)));
      print "\n";
    }
    print p(qq{To find additional similar structures, try the "Find Similar Assemblies" links at RCSB.});
  }
  print p("Or try", a({-href => "pdbBlast.cgi"}, "another search")), "\n";
}
finish_page();


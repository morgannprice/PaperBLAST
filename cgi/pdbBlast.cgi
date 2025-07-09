#!/usr/bin/perl -w

# Optional CGI arguments:
# query -- the query sequence in FASTA or UniProt or plain format, or, an identifier (as in litSearch.cgi)
# subject -- should be something like "7quc_A". If specified, show the alignment of the query to this sequence.

use strict;
use strict;
use CGI qw(:standard Vars start_ul);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use Time::HiRes qw{gettimeofday};
use DBI;
use LWP::Simple qw{get};
use HTML::Entities;
use IO::Handle; # for autoflush
use URI::Escape;
use Bio::SearchIO; # parsing alignments from bl2seq
use lib "../lib";
use pbweb;

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
my $fastacmd = "../bin/blast/fastacmd";
my $bl2seq = "../bin/blast/bl2seq";
foreach my $x ($blastall, $fastacmd, $bl2seq) {
  die "No such executable: $x\n" unless -x $x;
}

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

my $subjectSpec = $cgi->param('subject') || "";
die "Invalid subject" unless $subjectSpec =~ m/^[0-9a-zA-Z_]*$/;

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
} else { # has query
  die "No sequence to search" unless $seq;
  my $initial = substr($seq, 0, 10);
  my $seqlen = length($seq);
  $initial .= "..." if $seqlen > 10;
  $initial = "($seqlen a.a., $initial)";
  $initial = "" if ! $hasDef;

  my $subjectShow = $subjectSpec;
  if ($subjectSpec) {
    $subjectShow =~ s/_//g;
    print h3("Align", HTML::Entities::encode($def), $initial, "to $subjectShow"),
      "\n";
  } else { # has query and subject
    print
      h3("Comparing PDB to", HTML::Entities::encode($def), $initial),
      "\n";
  }
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
  open(my $fh, ">", $seqFile) || die "Cannot write to $seqFile";
  print $fh ">query\n", $seq, "\n";
  close($fh) || die "Error writing to $seqFile";

  my $hitsFile = "$tmpDir/$filename.pdbhits";
  if ($subjectSpec) {
    my $subjectFile = "$tmpDir/$filename.subject.fasta";
    my $fastacmd = qq{$fastacmd -s $subjectSpec -d $pdbDb -o $subjectFile};
    system($fastacmd) == 0 || die "fastacmd failed: $fastacmd\n$!";

    # Get information about the subject
    open(my $fh, "<", $subjectFile) || die "Cannot read $subjectFile";
    my $headerS = <$fh>;
    chomp $headerS;
    my $sDesc = $headerS; $sDesc =~ s/^\S+ +//;
    my ($id, $chain) = split /_/, $subjectSpec;
    my $qSeq = "";
    while (my $line = <$fh>) {
      chomp $line;
      $qSeq .= $line;
    }
    close($fh) || die "Error reading $subjectFile";

    print p("Subject: ",
            a({ -href => "https://rcsb.org/sequence/$id#$chain", -style => "text-decoration: none;" },
              HTML::Entities::encode($sDesc))), "\n";
    my $cmd = qq{$bl2seq -p blastp -i $seqFile -j $subjectFile -e 1e-5 -F "m S" -o $hitsFile > /dev/null 2>&1};
    system($cmd) == 0 || die "bl2seq failed: $cmd\n$!";
    my $searchio = Bio::SearchIO->new(-format => 'blast', -file => $hitsFile)
      || die "Failed to read $hitsFile";
    my $result = $searchio->next_result;
    my $hit = $result->next_hit if $result;
    if ($hit) {
      my $hsp = $hit->hsp;
      my ($queryBeg, $hitBeg) = $hsp->start('list');
      my ($queryEnd, $hitEnd) = $hsp->end('list');
      my $aln = $hsp->get_aln();
      my $alnQ = $aln->get_seq_by_pos(1)->seq;
      my $alnS = $aln->get_seq_by_pos(2)->seq;
      my $seqQ = $alnQ; $seqQ =~ s/-//g;
      die "BLAST parsing error, query sequences do not line up:\n$seqQ\n"
        . substr($seq,$queryBeg-1,$queryEnd-$queryBeg+1)
        unless $seqQ eq substr($seq, $queryBeg-1, $queryEnd-$queryBeg+1);
    print p(sprintf("%.1f%% identity over %d:%d / %d of query and %d:%d / %d of subject (%.1f bits)",
                    $hsp->percent_identity,
                    $queryBeg, $queryEnd, length($seq),
                    $hitBeg, $hitEnd, length($qSeq), $hsp->bits)), "\n";

      my $qAt = $queryBeg; # 1-based
      my $sAt = $hitBeg; # 1-based
      print qq{<PRE style="line-height: 1.4;">}, "\n";
      my $width = 60;
      for (my $i = 0; $i < length($alnQ); $i += $width) {
        my $qLine = "Query  " . sprintf("%4d", $qAt) . " ";
        my $sLine = "Subject" . sprintf("%4d", $sAt) . " ";
        for (my $j = 0; $j < $width && $i+$j < length($alnQ); $j++) {
          my $qChar = substr($alnQ, $i+$j, 1);
          my $sChar = substr($alnS, $i+$j, 1);
          $qLine .= coloredAAChar($qChar, $qAt, "query");
          $sLine .= coloredAAChar($sChar, $sAt, "subject");
          $qAt++ unless $qChar eq "-";
          $sAt++ unless $sChar eq "-";
        }
        print $qLine, sprintf("%4d", $qAt-1), "\n",
              $sLine, sprintf("%4d", $sAt-1), "\n",
              "\n";
      }
      print "</PRE>", "\n";
      print p(qq{The residue numbers for the subject reflect the actual position in the sequence from RCSB/PDB.
                       The author-provided "auth" coordinates, which are often used in papers and by structure viewers,
                       may not match. The bar for the UNIPROT line on the
                       <A HREF="https://rcsb.org/sequence/$id#$chain">sequence view</A> at RCSB will
                       usually show how the two relate.}), "\n";
    } else { # no hit
      print p("Sorry, BLAST did not find a significant alignment (at E &le; 1e-5)") . "\n";
    }
    if (param('debug')) {
      print "<PRE>\n"; system("cat $hitsFile"); print "</PRE>\n";
    }
    unlink($hitsFile);
    unlink($subjectFile);
    # end compare to one subject
  } else { # database search
    my $cmd = qq{$blastall -p blastp -i $seqFile -d $pdbDb -e 1e-3 -F "m S" -m 8 -o $hitsFile > /dev/null 2>&1};
    system($cmd) == 0 || die "blast failed: $cmd\n$!";
    open($fh, "<", $hitsFile) || die "Cannot read $hitsFile";
    my @hits = ();
    while (my $line = <$fh>) {
      my @F = split /\t/, $line;
      push @hits, \@F;
    }
    close($fh) || die "Error reading $hitsFile";

    my $dbDesc = "non-redundant PDB (which was clustered at 80% identity)";
    if (@hits == 0) {
      print p("Sorry, no significant homologs were found using BLAST against $dbDesc with E &le; 1e-3");
    } else { # has hits
      print p(scalar(@hits), "hits were found against $dbDesc");
      foreach my $hit (@hits) {
        my (undef, $subject, $identity, $alen, $mm, $gap, $qbeg, $qend, $sbeg, $send, $evalue, $bits) = @$hit;
        $subject =~ m/^([0-9a-zA-Z]+)_([0-9a-zA-Z]+)$/ || die "invalid subject $subject";
        my ($id, $chain) = ($1,$2);
        my ($desc) = $dbh->selectrow_array("SELECT desc FROM PdbClustInfo WHERE id = ? AND chain = ?",
                                           {
                                           }, $id, $chain);
        $desc = "" if !defined $desc;
        my $cov = ($qend-$qbeg+1)/length($seq);
        my $alnURL = "pdbBlast.cgi?query=" . uri_escape($query) . "&subject=$subject";
        my $alnInfo = sprintf(qq{<B>%0.f%%</B> identity, <B>%0.f%%</B> coverage of query,}
                              . qq{ %.1f bits, E = %g},
                              $identity, $cov*100, $bits, $evalue);
        print p(a({-href => "https://rcsb.org/structure/$id",
                   -style => "text-decoration: none;",
                   -title => "see structure at RCSB"}, $id),
                a({-href => "https://rcsb.org/sequence/$id#$chain",
                   -style => "text-decoration: none;",
                   -title => "see sequence at RCSB"}, $chain),
                b(HTML::Entities::encode($desc)),
                br(),
                small(a({-href => $alnURL, -style =>"text-decoration: none;",
                         -title => "$qbeg:$qend of query aligns to $sbeg:$send of $id:$chain"},
                        $alnInfo)));
        print "\n";
      }
      print p(qq{To find additional similar structures, try the "Find Similar Assemblies" links at RCSB.});
    } # end else has hits
  } # end database search
  unlink($hitsFile);
  unlink($seqFile);
  print p("Or try", a({-href => "pdbBlast.cgi"}, "another search")), "\n";
}
finish_page();


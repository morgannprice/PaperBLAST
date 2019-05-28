#!/usr/bin/perl -w
#######################################################
## showAlign.cgi
##
## Copyright (c) 2017 University of California
##
## Authors: Morgan Price
#######################################################
#
# Required CGI parameters:
# seq1 -- 1st sequence
# Either seq2 -- 2nd sequence or acc2 -- accession for 2nd sequence, to look it up in the database
#
# Optional CGI garameters:
# def1 -- name of 1st sequence
# def2 -- name of 2nd sequence
# debug

use strict;
use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use Time::HiRes qw{gettimeofday};
use DBI;
use lib "../lib";
use pbutils; # for ReadTable()
use pbweb; # for FetchFasta()

my $cgi=CGI->new;
my $def1 = $cgi->param('def1') || "query";
my $def2 = $cgi->param('def2') || "subject";
my $debug = $cgi->param('debug');

my $base = "../data";
my $sqldb = "$base/litsearch.db";
my $dbh = DBI->connect("dbi:SQLite:dbname=$sqldb","","",{ RaiseError => 1 }) || die $DBI::errstr;

my $title = "Align $def1 vs. $def2";
start_page('title' => $title);

my $tmpDir = "../tmp";
my $bl2seq = "../bin/blast/bl2seq";
die "No such executable: $bl2seq" unless -x $bl2seq;

my $procId = $$;
my $timestamp = int (gettimeofday() * 1000);
my $prefix = "$tmpDir/$procId$timestamp";

my $seq1 = $cgi->param('seq1') || die "No seq1 argument";
my $seq2 = $cgi->param('seq2');
my $acc2;
if (!defined $seq2) {
    $acc2 = $cgi->param('acc2') || die "Either seq2 or acc2 must be specified";
    # Find it in the database using fastacmd
    my $db = "$base/uniq.faa";
    die "No such file: $db" unless -e $db;
    $seq2 = FetchFasta($dbh, $db, $acc2)
       || die "Cannot fetch sequence for $acc2";
}
$seq1 =~ m/^[A-Z*]+$/ || die "Invalid seq1";
$seq1 =~ s/[*]//g;
$seq2 =~ m/^[A-Z*]+$/ || die "Invalid seq2";
$seq2 =~ s/[*]//g;
die "Empty seq1" if $seq1 eq "";
die "Empty seq2" if $seq2 eq "";

my $len1 = length($seq1);
my $len2 = length($seq2);
my $cdd_base = "http://www.ncbi.nlm.nih.gov/Structure/cdd/wrpsb.cgi?seqinput=>";
my $newline = "%0A";
print
    p("$def1: $len1 amino acids", br(),
      "$def2: $len2 amino acids"),
    p("Or search the Conserved Domains Database with",
      a({ -href => "$cdd_base$def1$newline$seq1" }, $def1),
      "or",
      a({ -href => "$cdd_base$def2$newline$seq2" }, $def2));

open(FAA1, ">", "$prefix.faa1");
print FAA1 ">$def1\n$seq1\n";
close(FAA1) || die "Error writing to $prefix.faa1";

open(FAA2, ">", "$prefix.faa2");
print FAA2 ">$def2\n$seq2\n";
close(FAA2) || die "Error writing to $prefix.faa2";

print p("Running $bl2seq with prefix $prefix") if $debug;
system($bl2seq,
       "-p", "blastp",
       "-e", 0.01,
       "-F", "m S",
       "-i", "$prefix.faa1",
       "-j", "$prefix.faa2",
       "-o", "$prefix.bl2seq") == 0 || die "bl2seq failed";
my @lines = ();
open(IN, "<", "$prefix.bl2seq") || die "Cannot read $prefix.bl2seq";
@lines = <IN>;
@lines = map { chomp; $_; } @lines;
close(IN) || die "Error reading $prefix.bl2seq";

unlink("$prefix.faa1");
unlink("$prefix.faa2");

# remove lines up to the first alignment ("Score")
do {
    shift @lines
} until $lines[0] =~ m/Score/i || $lines[0] =~ m/No hits/i || @lines == 1;
my @out = ();
foreach my $line (@lines) {
    if ($line =~ m/^Lambda\s+K\s+H$/) {
        last;
    } else {
        push @out, $line;
    }
}

print
    h3("Protein alignments from bl2seq"),
    pre(join("\n", @out));

$| = 1; # flush STDOUT

print p("acc2 = $acc2; out[0] = $out[0]")."\n" if $debug && defined $acc2;
if (defined $acc2 && $out[0] !~ m/No hits/i) {
  print p("Looking for uniprot id for $acc2")."\n" if $debug;
  # Make a feature report
  # First, find the UniProt identifier, if it is known
  my $uniprotId = undef;
  if ($acc2 =~ m/SwissProt::(.*)$/) {
    $uniprotId = $1;
  } else {
    # try to find acc2 via SeqToDuplicate
    my @dups = $dbh->selectrow_array("SELECT duplicate_id FROM SeqToDuplicate WHERE sequence_id = ?",
                                        {}, $acc2);
    unshift @dups, $acc2;
    foreach my $dup (@dups) {
      if ($dup =~ m/^[A-Z][0-9A-Z]+$/ && $dup !~ m/^VIMSS/) {
        $uniprotId = $dup;
      } elsif ($dup =~ m/(SwissProt|BRENDA)::(.*)$/) {
        $uniprotId = $2;
      } elsif ($dup =~ m/(metacyc|CharProtDB)::(.*)$/) {
        my ($db,$id) = ($1,$2);
        my ($id2) = $dbh->selectrow_array("SELECT id2 FROM CuratedGene WHERE db = ? AND protId = ?",
                                          {}, $db, $id);
        if ($db eq "metacyc" && $id2) {
          $uniprotId = $id2;
        } elsif ($db eq "CharProtDB" && $id2 =~ m/SP[|]/) {
          my %id2 = split /[|]/, $id2;
          $uniprotId = $id2{"SP"} if $id2{"SP"};
        }
      }
      last if $uniprotId;
    }
  }
  print p("uniprotId " . ($uniprotId || "undef")) if $debug;
  my %typeNames = ("ACT_SITE" => "Active site",
                   "BINDING" => "Binding site",
                   "CA_BIND" => "Calcium binding",
                   "ZN_FING" => "Zinc finger",
                   "METAL" => "Metal binding",
                   "DNA_BIND" => "DNA binding",
                   "NP_BIND" => "Nucleotide binding",
                   "DOMAIN" => "Domain",
                   "MOTIF" => "Motif",
                   "REGION" => "Other region",
                   "REPEAT" => "Low-complexity",
                   "PROPEP" => "Pro-peptide",
                   "SIGNAL" => "Signal sequence",
                   "TRANSIT" => "Transit peptide",
                   "TRANSMEM" => "Membrane-spanning",
                   "INTRAMEM" => "Intra-membrane",
                   "NON_STD" => "Non-standard residue",
                   "MOD_RES" => "Modification",
                   "CARBOHYD" => "Glycosylation",
                   "LIPID" => "Lipidation",
                   "DISULFID" => "Disulfide bond",
                   "CROSSLNK" => "Covalent cross-link",
                   "MUTAGEN" => "Mutagenesis",
                   "HELIX" => "Alpha helix",
                   "STRAND" => "Beta strand",
                   "COILED" => "Coiled coil",
                   "TURN" => "Hydrogen-bonded turn",
                   "SITE" => "Other residue");

  if ($uniprotId) {
    die "Invalid uniprotId $uniprotId" unless $uniprotId =~ m/^[A-Z0-9_]+$/;
    if (system("../bin/alnreport.pl $uniprotId $prefix.bl2seq > $prefix.alnreport") != 0) {
      print p(small("No feature table because alnreport.pl failed for $uniprotId -- is bioperl installed?"));
      unlink("$prefix.alnreport");
    } else {
      my @alnreport = ReadTable("$prefix.alnreport");
      unlink("$prefix.alnreport");
      if (@alnreport > 0) {
        my @trows = Tr({ -valign => "bottom", -align => "center" },
                       th( ['Type', 'Comment', '# Aligned', '# Conserved', 'Sbjct Pos.', 'Sbjct Seq.', 'Query Pos.', 'Query Seq.'] ));
        foreach my $row (@alnreport) {
          my ($qpos, $spos);
          if ($row->{sbegin} eq $row->{send}) {
            $spos = $row->{sbegin};
            $qpos = $row->{qbegin};
          } else {
            my $joinby = $row->{type} eq "DISULFID" ? "," : "..";
            $spos = $row->{sbegin} . $joinby . $row->{send};
            $qpos = $row->{qbegin} eq "" ? "" : $row->{qbegin} . $joinby . $row->{qend};
          }
          my $slen = $row->{send} - $row->{sbegin} + 1;
          $slen = 2 if $slen > 1 && $row->{type} eq "DISULFID";
          my $short = length($row->{subjectSeq}) <= 10;
          $row->{nAligned} = 0 if $row->{nAligned} eq "";
          $row->{nMatch} = 0 if $row->{nMatch} eq "";
          my $comment = $row->{comment};
          $comment =~ s![(]by similarity[)]!<small>(by similarity)</small>!g;
          $comment =~ s! ?/FTId=.*$!!;
          $comment =~ s/[.] *$// if $comment =~ m/^[^.]+[.] *$/;
          my $match = $row->{nAligned} > 0 ? "$row->{nMatch} / $slen" : "--";
          push @trows, Tr({ -valign => 'top', -align => "left" },
                          td($typeNames{$row->{type}} || $row->{type}),
                          td($comment),
                          td("$row->{nAligned} / $slen"),
                          td($match),
                          td($spos),
                          td($short ? $row->{subjectSeq} : "..."),
                          td($qpos),
                          td($short ? $row->{querySeq} : "..."));
        }
        print h3("Conservation of Functional Sites and Regions"),
          table({cellspacing => 0, cellpadding => 3, border => 1}, @trows),
            p(small("Conservation of each site or region is based on the highest-scoring alignment only.")),
              p(a({ -href => "http://uniprot.org/uniprot/$uniprotId" }, "Also see UniProt entry for $def2" ));
      } else {
        print p(small("No functional sites for the subject",
                    a({ -href => "http://www.uniprot.org/uniprot/$uniprotId" }, $uniprotId)));
      }
    }
  } else {
    print p(small("No functional sites because the subject $def2 is not linked to a UniProt identifier in PaperBLAST's database"));
  }
}
unlink("$prefix.bl2seq");


#!/usr/bin/perl -w
#######################################################
## curatedSearch.cgi -- curated proteins that are
##      relevant to a text query
##
## Copyright (c) 2018 University of California
##
## Authors: Morgan Price
#######################################################
#
# Optional CGI parameters:
# query -- what term to search for
# word -- if non-empty, report whole word matches only
# faa -- if set, returns the fasta sequences instead of showing a table of descriptions

use strict;
use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use Time::HiRes qw{gettimeofday};
use DBI;
use lib "../lib";
use pbutils; # for ReadFastaEntry(), Curated search functions
use pbweb qw{start_page finish_page loggerjs AddCuratedInfo};
use URI::Escape;
use HTML::Entities;

my $base = "../data";
my $sqldb = "$base/litsearch.db";
my $dbh = DBI->connect("dbi:SQLite:dbname=$sqldb","","",{ RaiseError => 1 }) || die $DBI::errstr;
my $blastdir = "../bin/blast";
die "No such directory: $blastdir" unless -d $blastdir;
my $blastdb = "$base/uniq.faa";
die "No such file: $blastdb" unless -e $blastdb;

my $cgi=CGI->new;
my $query = $cgi->param('query');
my $word = $cgi->param('word');
my $faa_mode = $cgi->param('faa');
my $table_mode = ! $faa_mode;

&start_page('title' => "Search for Curated Proteins",
            'banner' => "Curated BLAST for Genomes",
            'bannerURL' => "genomeSearch.cgi")
  if $table_mode;

if ($query) {
  my $maxhits = 1000;
  my $chits = CuratedMatch($dbh, $query, $maxhits+1);
  my $quotedquery = HTML::Entities::encode($query);
  if ($table_mode && @$chits > $maxhits) {
    print p(qq{Sorry, too many curated entries match the query '$quotedquery'. Please try},
            a({ -href => "curatedSearch.cgi" }, "another query").".");
    finish_page();
  }
  if ($table_mode && @$chits == 0) {
    print p(qq{None of the curated entries in PaperBLAST's database match '$quotedquery'. Please try},
            a({ -href => "curatedSearch.cgi" }, "another query") . ".");
    finish_page();
  }
  if ($word) {
    $chits = CuratedWordMatch($chits, $query);
    if ($table_mode && @$chits == 0) {
      print p(qq{None of the curated entries in PaperBLAST's database match '$query' as complete words. Please try},
              a({ -href => "curatedSearch.cgi" }, "another query") . ".");
      finish_page();
    }
  }
  my $wordstatement = $word ? " as complete word(s)" : "";
  print p("Found", scalar(@$chits),
          qq{curated entries in PaperBLAST's database that match '$query'${wordstatement}.},
          "Try",
          a({-href => "curatedSearch.cgi"}, "another search"))
    if $table_mode;

  my %idToChit = (); # sequence identifier to curated gene hit(s)
  foreach my $hit (@$chits) {
    my $seqid = $hit->{db} . "::" . $hit->{protId};
    my $uniqid = IdToUniqId($dbh, $seqid);
    push @{ $idToChit{$uniqid} }, $hit;
  }
  my @uniqids = sort keys %idToChit;

  if ($faa_mode) {
    my $tmpdir = $ENV{TMP} || "/tmp";
    die "Not a directory: $tmpdir" unless -d $tmpdir;
    my $tmpfile = "$tmpdir/curated.$$.faa";
    FetchSeqs($blastdir, $blastdb, \@uniqids, $tmpfile);
    print header('text/x-fasta'); # content type
    # Should improve this to show deflines!

    open(my $fh, "<", $tmpfile) || die "Cannot read $tmpfile";
    my $state = {};
    while (my ($header, $sequence) = ReadFastaEntry($fh, $state)) {
      $header =~ m/^lcl[|](\S+) /
        || die "Cannot parse fasta header $header";
      my $uniqid = $1;
      die "Unexpected fasta identifier $uniqid from $header"
        unless exists $idToChit{$uniqid};
      my @ids = ();
      my @descs = ();
      foreach my $chit (@{ $idToChit{$uniqid}}) {
        push @ids, $chit->{db} . "::" . $chit->{protId};
        push @descs, $chit->{desc};
      }
      print ">" . join(" ", @ids) . " " . join(";; ", @descs) . "\n" . $sequence . "\n";
    }
    close($fh) || die "Error reading $tmpfile";
    unlink($tmpfile);
    exit(0);
  } else {
    # table_mode
    print p("These curated entries have", scalar(@uniqids), "distinct",
            a({-href => join("&", "curatedSearch.cgi?faa=1",
                             "query=" . uri_escape($query),
                             "word=" . ($word ? 1 : 0)) },
              "sequences") . ".");
    my %uniqDesc = map { $_ => $idToChit{$_}[0]{desc} } @uniqids;
    my @sorted = sort { $uniqDesc{$a} cmp $uniqDesc{$b} } @uniqids;
    foreach my $uniqid (@sorted) {
      my @descs = ();
      foreach my $chit (@{ $idToChit{$uniqid} }) {
        AddCuratedInfo($chit);
        my @showOrgWords = split / /, $chit->{organism};
        @showOrgWords = @showOrgWords[0..1] if @showOrgWords > 2;
        push @descs, a({-href => $chit->{URL},
                        -title => "from " . $chit->{db},
                        -onmousedown => loggerjs("curated", $chit->{subjectId}) },
                       $chit->{showName}) . ": " . $chit->{desc}
                         . " " . small("from " . i(join(" ", @showOrgWords)));
      }
      print p(join("<BR>", @descs));
    }
    finish_page(); # finish showing the hits
  }
} else {
  # show form
  print
    start_form(-method => 'get', -action => 'curatedSearch.cgi'),
    p(small("Example:", a({ -href => "curatedSearch.cgi?query=perchlorate"}, "perchlorate"))),
    p("Search for:", textfield(-name => 'query', -value => '', -size => 50, -maxlength => 200)),
    p({-style => "margin-left: 3em;" },
      checkbox(-name => "word", -label => "Match whole words only?")),
    p(submit(-name => 'Search')),
    end_form;
}
finish_page();

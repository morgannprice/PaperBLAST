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
sub TransporterMatch($$$);

my $base = "../data";
my $sqldb = "$base/litsearch.db";
my $dbh = DBI->connect("dbi:SQLite:dbname=$sqldb","","",{ RaiseError => 1 }) || die $DBI::errstr;
my $blastdir = "../bin/blast";
die "No such directory: $blastdir" unless -d $blastdir;
my $blastdb = "$base/uniq.faa";
die "No such file: $blastdb" unless -e $blastdb;

my $staticDir = "../static";

my $cgi=CGI->new;
my $query = $cgi->param('query');
my $queryShow = $query;
my $wordMode = $cgi->param('word') || 0;
my $transporterMode;
my $faa_mode = $cgi->param('faa');
my $table_mode = ! $faa_mode;

&start_page('title' => "Search for Curated Proteins",
            'banner' => "Curated BLAST for Genomes",
            'bannerURL' => "genomeSearch.cgi")
  if $table_mode;

if ($query) {
  my $maxhits = 1000;
  my $chits;
  if ($query =~ m/^transporter:(.+)$/) {
    $queryShow = $1;
    $wordMode = 0;
    $transporterMode = 1;
    $chits = TransporterMatch($dbh, $staticDir, $queryShow);
    $queryShow =~ s!:! / !g;
  } else {
    $chits = CuratedMatch($dbh, $query, $maxhits+1);
  }
  my $quotedquery = HTML::Entities::encode($queryShow);
  my $transporterStatement = $transporterMode ? " as transporters" : "";
  if ($table_mode && @$chits > $maxhits) {
    print p(qq{Sorry, too many curated entries match the query '$quotedquery'${transporterStatement}. Please try},
            a({ -href => "curatedSearch.cgi" }, "another query").".");
    finish_page();
  }
  if ($table_mode && @$chits == 0) {
    print p(qq{None of the curated entries in PaperBLAST's database match '$quotedquery'${transporterStatement}. Please try},
            a({ -href => "curatedSearch.cgi" }, "another query") . ".");
    finish_page();
  }
  if ($wordMode) {
    $chits = CuratedWordMatch($chits, $query);
    if ($table_mode && @$chits == 0) {
      print p(qq{None of the curated entries in PaperBLAST's database match '$quotedquery' as complete words. Please try},
              a({ -href => "curatedSearch.cgi" }, "another query") . ".");
      finish_page();
    }
  }
  my $wordStatement = "";
  $wordStatement = " as complete word(s)" if $wordMode;
  $wordStatement = " as transporters" if $transporterMode;
  print p("Found", scalar(@$chits),
          qq{curated entries in PaperBLAST's database that match '$quotedquery'${wordStatement}.},
          "Or try",
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
                             "word=" . ($wordMode ? 1 : 0)) },
              "sequences") . ".",
              a({-href => "curatedClusters.cgi?query=$quotedquery&word=$wordMode"},
                "Cluster these sequences"),
              "(recent entries may not be included in clustering results).");

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

sub TransporterMatch($$$) {
  my ($dbh, $staticDir, $compoundSpec) = @_;
  my @compoundList = split /:/, $compoundSpec;
  my %compoundList = map { $_ => 1, lc($_) => 1 } @compoundList;

  # metacyc files linking reactions to compounds and to protein(s)
  my $reactionCompoundsFile = "$staticDir/metacyc.reaction_compounds";
  my $reactionProteinsFile = "$staticDir/metacyc.reaction_links";

  my %rxnCompounds = (); # rxnId => compounds with side, coefficient, compartment, compoundId, and compoundName
  my %rxnLocation = ();
  open(my $fhc, "<", $reactionCompoundsFile) || die "Cannot read $reactionCompoundsFile";
  while (my $line = <$fhc>) {
    chomp $line;
    my @F = split /\t/, $line;
    my $rxnId = shift @F;
    my $rxnLoc = shift @F;
    $rxnLocation{$rxnId} = $rxnLoc;
    my @cmp = ();
    foreach my $f (@F) {
      my @pieces = split /:/, $f;
      my $side = shift @pieces;
      my $coefficient = shift @pieces;
      my $compartment = shift @pieces;
      my $compoundId = shift @pieces;
      my $compoundName = join(":", @pieces) || "";
      push @cmp, { 'side' => $side, 'coefficient' => $coefficient,
                   'compartment' => $compartment,
                   'compoundId' => $compoundId, 'compoundName' => $compoundName };
    }
    $rxnCompounds{$rxnId} = \@cmp;
  }
  close($fhc) || die "Error reading $reactionCompoundsFile";

  my %rxnProt = (); # rxnId => list of subunits of the form metacyc::id
  my %rxnName = ();
  open (my $fhp, "<", $reactionProteinsFile) || die "Cannot read $reactionProteinsFile";
  while (my $line = <$fhp>) {
    chomp $line;
    my @F = split /\t/, $line;
    my $rxnId = shift @F;
    my $rxnName = shift @F;
    push @{ $rxnProt{$rxnId} }, @F;
    $rxnName{$rxnId} = $rxnName;
  }
  close($fhp) || die "Error reading $reactionProteinsFile";

  my @out = (); # list of rows from CuratedGene
  my %metacycSaved = (); # protId => 1 for entries already in @out

  # Find the MetaCyc reactions that match the compound, either as compoundId or compoundName,
  # and determine if they are transport reactions, that is, the compound appears
  # with two different locations.
  # (For now ignore PTS)
  foreach my $rxnId (keys %rxnCompounds) {
    my $cmps = $rxnCompounds{$rxnId};
    my @cmpMatch = grep { exists $compoundList{$_->{compoundId}}
                            || exists $compoundList{$_->{compoundName}}
                              || exists $compoundList{lc($_->{compoundName})} } @$cmps;
    my %compartments = map { $_->{compartment} => 1 } @cmpMatch;
    my @compartments = grep { $_ ne "" } (keys %compartments);
    my %compartmentsAll = map { $_->{compartment} => 1 } @$cmps;
    my @compartmentsAll = grep { $_ ne "" } (keys %compartmentsAll);
    # The same compound on both sides, or, this compound on one side and something else on the other
    my $use = (@cmpMatch > 1 && @compartments > 1)
      || (@cmpMatch > 0 && @compartments > 0 && @compartmentsAll > 1);
    if ($use) {
      # Compound of interest in more than one compartment
      # So find all the proteins for this reaction
      my $protids = $rxnProt{$rxnId};
      foreach my $protid (@$protids) {
        my $id = $protid; $id =~ s/^metacyc:://;
        next if exists $metacycSaved{$id};
        $metacycSaved{$id} = 1;
        my $chit = $dbh->selectrow_hashref(qq{ SELECT * FROM CuratedGene WHERE db="metacyc" AND protId=? },
                                           {}, $id);
        if (defined $chit) {
          push @out, $chit;
        } else {
          print p("Skipped unknown metacyc id $id");
        }
      }
    }
  }

  my $tcdbs = $dbh->selectall_arrayref("SELECT * FROM CuratedGene WHERE db = 'TCDB';",
                                       { Slice => {} });
  foreach my $tcdb (@$tcdbs) {
    my @comments = split /_:::_/, $tcdb->{comment};
    my $keep = 0;
    foreach my $comment (@comments) {
      if ($comment =~ m/^SUBSTRATES: (.*)$/) {
        my @substrates = split /, /, $1;
        foreach my $substrate (@substrates) {
          $keep = 1 if exists $compoundList{lc($substrate)};
        }
      }
    }
    push @out, $tcdb if $keep;
  }

  # Search the description for words transport, transporter, exporter, porter, permease, import
  my $chits = $dbh->selectall_arrayref(
    qq{SELECT * FROM CuratedGene WHERE db IN ("SwissProt","CharProtDB","BRENDA","reanno","ecocyc")
       AND (desc LIKE "%transport%"
            OR desc LIKE "%porter%"
            OR desc LIKE "%import%"
            OR desc LIKE "%permease%"
            OR desc like "%PTS system%")},
    { Slice => {} });
  foreach my $chit (@$chits) {
    my $desc = $chit->{desc};
    my $keep = 0;
    foreach my $cmp (@compoundList) {
      my $pattern = quotemeta($cmp);
      # doing a word-based check this is tricky because substrates may be separated by "/", without spaces
      # or may appear as "glucose-binding"
      # or because another term like "6-phosphate" could be present as the next word.
      # That last issue is not handled.
      my $b = "[,/ ]"; # pattern for word boundary
      if ($desc =~ m!^$pattern$b!i # at beginning
          || $desc =~ m!$b$pattern$b!i # in middle
          || $desc =~ m!$b$pattern$!i # at end
          || $desc =~ m!^$pattern-(binding|specific)!i # at beginning
          || $desc =~ m!$b$pattern-(binding|specific)!i) { # in middle
        $keep = 1;
      }
    }
    push @out, $chit if $keep;
  }

  return \@out;
}

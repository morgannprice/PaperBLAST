#!/usr/bin/perl -w
# treeSites.cgi -- show key sites along a tree
use strict;
use CGI qw(:standard Vars start_ul);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use URI::Escape;
use HTML::Entities;
use Digest::MD5 qw{md5_hex};
use List::Util qw{sum min max};
use IO::Handle; # for autoflush
use lib "../lib";
use pbutils qw{ReadFastaEntry ParseClustal ParseStockholm seqPosToAlnPos};
use pbweb;
use MOTree;
use DBI;

# In rendering mode, these options are required:
# alnFile or alnId -- alignment in fasta format, or the file id (usually, md5 hash)
#   In header lines, anything after the initial space is assumed to be a description.
#   Either "." or "-" are gap characters.
# treeFile or treeId -- tree in newick format, or the file id (usually, md5 hash)
#   Labels on internal nodes are ignored.
#   Multi-furcations are allowed. It is treated as rooted.
# Optional arguments:
# anchor -- sequence to choose positions from
# pos -- comma-delimited list of positions (in anchor if set, or else, in alignment)
#	ranges like 90:96 may also be used
# tsvFile or tsvId -- descriptions for the ids (as uploaded file or file id)
# (This tool automatically saves uploads using alnId, treeId, or tsvId, under their md5 hash)
# zoom -- an internal node (as numbered by MOTree) to zoom into
#
# In pattern search mode, required options are:
# alnId
# pattern -- the pattern to search for
# All the other arguments from rendering mode are preserved if present
#	(and if anchor is set, hits to that sequence are shown first)
#
# Tree building mode -- an alignment is provided, but not a tree, and buildTree is set
# Optional arguments:
# trimGaps -- set if positions that are >=50% gaps should be trimmed off
# trimLower -- set if positions that are >=50% lower case should be trimmed off
#
# Tree options mode -- an alignment is provided, but not a tree, and buildTree is not set
#
# Alignment building mode:
# seqsFile or seqsId -- sequences (not aligned) in fasta format
# addSeq -- sequence identifier(s) to add
# buildAln -- if set, does the alignment
#
# If neither an alignment nor sequences are provided, it shows an input form.
#
# Find homologs mode:
# query -- the query (fasta format, uniprot format, raw sequence, or an identifier)

sub handleTsvLines; # handle tab-delimited description lines

# id and type to an open filehandle, usually from ../tmp/aln/id.type,
# but it also looks in ../static/
# Also it verifies that the id is harmless (alphanumeric, _, or - characters only)
sub openFile($$);
sub findFile($$); # similar but returns the filename

# Given a list of lines and a type, saves it to a file in $tmpDir, if necessary,
# and returns the hash id
sub savedHash($$);

# maximum size of posted data, in bytes
my $maxMB = 25;
$CGI::POST_MAX = $maxMB*1024*1024;
my $maxN = 2000; # maximum number of sequences in alignment or leaves in tree
my $maxNAlign = 200; # maximum number of sequences to align
my $maxLenAlign = 10000;

my $tmpDir = "../tmp/aln";
my $tmpPre = "$tmpDir/treeSites.$$";
my $base = "../data";
my $blastdb = "$base/uniq.faa";
my $hassitesdb = "$base/hassites.faa";
my $nCPU = 4;
my $blastall = "../bin/blast/blastall";
my $fastacmd = "../bin/blast/fastacmd";
my $sqldb = "$base/litsearch.db";
my $dbh = DBI->connect("dbi:SQLite:dbname=$sqldb","","",{ RaiseError => 1 }) || die $DBI::errstr;

# A symbolic link to the Fitness Browser data directory is used (if it exists)
# to allow quick access to proteins from the fitness browser.
# That directory must include feba.db (sqlite3 database) and aaseqs (in fasta format)
my $fbdata = "../fbrowse_data"; # path relative to the cgi directory

# The "taylor" color scheme (no longer used) is based on
# "Residual colours: a proposal for aminochromography" (Taylor 1997)
# and https://github.com/omarwagih/ggseqlogo/blob/master/R/col_schemes.r
# I added a dark-ish grey for gaps.
my %taylor = split /\s+/,
  qq{D #FB9A99 E #E31A1C N #B2DF8A Q #95CF73 K #78C05C R #58B044 H #33A02C F #FDBF6F W #FFA043 Y #FF7F00 P #FFFF99 M #DBAB5E C #B15928 G #A6CEE3 A #8AB8D7 V #6DA2CC L #4D8DC0 I #1F78B4 S #CAB2D6 T #6A3D9A - #555555};

# Used color brewer to get 6 sets of color pairs and interpolate between them within groups of related a.a.
# The 5 colors within green (GAVLI) and blue (NQKRH) were too difficult to tell apart, so changed
# the lightest end of these ranges to be more pale.
# It's still a bit difficult to distinguish arg/his or leu/ile.
my %colors = split /\s+/,
  qq{D #FB9A99 E #E31A1C N #CCE6FF Q #A7C9EC K #81ADD9 R #5892C7 H #1F78B4 F #FDBF6F W #FFA043 Y #FF7F00 P #FFFF99 M #DBAB5E C #B15928 G #E6FFE6 A #BDE7B7 V #93D089 L #68B85C I #33A02C S #CAB2D6 T #6A3D9A - #555555};

print
  header(-charset => 'utf-8'),
  start_html(-title => "Sites on a Tree",
             -script => [{ -type => "text/javascript", -src => "../static/treeSites.js"}]),
  h2("Sites on a Tree"),
  "\n";
autoflush STDOUT 1; # show preliminary results

my $query = param('query');

if (defined $query && $query ne "") {
  # Homologs mode, with 1-sequence input
  my ($id, $seq) = parseSequenceQuery(-query => $query,
                                      -dbh => $dbh,
                                      -blastdb => $blastdb,
                                      -fbdata => $fbdata);
  fail("No sequence") unless defined $seq;

  # split off the description
  my $desc = "";
  if ($id =~ m/^(\S+) (.*)/) {
    $id = $1;
    $desc = $2;
  }
  fail("Sorry, ,(): are not allowed in identifiers")
    if $id =~ m/[():;]/;

  $id = sequenceToHeader($seq)
    if defined $seq && $id eq "";

  my $maxHits = 100;
  print p("Searching for up to $maxHits curated homologs for",
          encode_entities($id),
          encode_entities($desc),
          "(". length($seq) . " a.a.)"), "\n";

  my $tmpFaa = "$tmpPre.faa";
  open(my $fhFaa, ">", $tmpFaa) || die "Cannot write to $tmpPre.faa";
  print $fhFaa ">", $id, "\n", $seq, "\n";
  close($fhFaa) || die "Error writing to $tmpFaa";

  die "No such executable: $blastall" unless -x $blastall;
  my @hits = ();

  my $tmpHits = "$tmpPre.hits";
  foreach my $db ($blastdb, $hassitesdb) {
    die "No such file: $db" unless -e $db;
    my @cmd = ($blastall, "-p", "blastp", "-d", $db, "-i", $tmpFaa, "-o", $tmpHits,
               "-e", 0.001, "-m", 8, "-a", $nCPU, "-F", "m S");
    system(@cmd) == 0 || die "Error running blastall: $!";
    open(my $fhHits, "<", $tmpHits) || die "Cannot read $tmpHits";
    while (<$fhHits>) {
      chomp;
      my @F = split /\t/, $_;
      push @hits, \@F;
    }
    close($fhHits) || die "Error reading $tmpHits";
  }
  unlink($tmpHits);
  unlink($tmpFaa);

  # Sort hits by bit score
  @hits = sort { $b->[11] <=> $a->[11] } @hits;

  # Keep non-duplicate hits to curated
  my @keep = (); # list of rows including subject, identity, qbeg, qend, sbeg, send, bits, and evalue
  my %seen = ();
  foreach my $hit (@hits) {
    my ($queryId2, $subject, $identity, $alen, $mm, $gaps, $qbeg, $qend, $sbeg, $send, $evalue, $bits) = @$hit;
    die unless defined $bits;
    next if $identity < 30 || ($qend - $qbeg + 1) < length($seq) * 0.7;

    if ($subject =~ m/SwissProt:(.*)/) {
      # Figure out if this subject has a corresponding uniq id and use that instead
      my $uniprotId = $1;
      my $id = "SwissProt::$uniprotId";
      my ($len2) = $dbh->selectrow_array("SELECT protein_length FROM CuratedGene WHERE db='SwissProt' AND protId=?",
                                         {}, $uniprotId);
      my $uniqId;
      if (defined $len2) {
        ($uniqId) = $dbh->selectrow_array("SELECT sequence_id FROM SeqToDuplicate WHERE duplicate_id = ?",
                                             {}, "SwissProt::$uniprotId");
        $uniqId = "SwissProt::$uniprotId" if !defined $uniqId; # maps to itself if no record of duplicate
        $subject = $uniqId;
      }
    }

    next if exists $seen{$subject};

    $seen{$subject} = 1;

    # Figure out if this subject is curated, if it's not a hassites or curated id
    unless ($subject =~ m/^[a-zA-Z]+:/) {
      my $dupIds = $dbh->selectcol_arrayref("SELECT duplicate_id FROM SeqToDuplicate WHERE sequence_id = ?",
                                            {}, $subject);
      my $keep = 0;
      foreach my $dupId (@$dupIds) {
        $keep = 1 if $dupId =~ m/:/;
      }
      next unless $keep;
    }

    push @keep, { 'subject' => $subject, 'identity' => $identity,
                  'qbeg' => $qbeg, 'qend' => $qend,
                  'sbeg' => $sbeg, 'send' => $send,
                  'evalue' => $evalue, 'bits' => $bits };
  }

  fail("Sorry, no hits to curated proteins at above 30% identity and 70% coverage")
    if scalar(@keep) == 0;

  print p("Found hits above 30% identity and 70% coverage to", scalar(@keep), " curated proteins");

  my $minIdentity = min(map $_->{identity}, @keep);
  if ($minIdentity == 100) {
    fail("All hits are identical to the query");
  } elsif ($minIdentity >= 95) {
    print p("All hits are nearly identical to the query");
  }

  my $nOld = scalar(@keep);
  @keep = grep $_->{identity} < 100, @keep;
  print p("Removed hits that are identical to the query, leaving", scalar(@keep))
    if scalar(@keep) < $nOld;

  # For now, just show a table of results, and add a link to align them
  foreach my $row (@keep) {
    print p($row->{identity}."% identity", $row->{subject}, );
  }

  # Fetch the sequences
  foreach my $row (@keep) {
    my $subject = $row->{subject};
    my $uniqId = $subject;
    my @fasta = ();
    my $subjectDesc;
    if ($uniqId =~ m/^([a-zA-Z]+):([a-zA-Z0-9].*)$/) {
      my ($db, $id) = ($1,$2);
      my $chain = "";
      my @pieces = split /:/, $id;
      if (@pieces == 2) {
        ($id, $chain) = @pieces;
      } elsif (@pieces == 1) {
        ;
      } else {
        die "Cannot handle identifier $subject";
      }

      # From hassites, not from the main database
      my $tmpFile = "$tmpPre.fastacmd";
      die "No such command: $fastacmd" unless -x $fastacmd;
      system($fastacmd, "-s", $subject, "-o", $tmpFile, "-d", $hassitesdb) == 0
        || die "fastacmd failed to find $subject in $hassitesdb";
      open(my $fh, "<", $tmpFile) || die "Cannot read $tmpFile";
      @fasta = <$fh>;
      close($fh) || die "Error reading $tmpFile";
      unlink($tmpFile);
      die "fastacmd failed for $subject in $hassitesdb" unless @fasta > 1;
      my $info = $dbh->selectrow_hashref("SELECT * from HasSites WHERE db = ? AND id = ? AND chain = ?",
                                         {}, $db, $id, $chain);
      $subjectDesc = $info->{desc};
    } else {
      my $fasta = DBToFasta($dbh, $blastdb, $uniqId);
      die "No sequence for $uniqId in $blastdb" unless defined $fasta;
      @fasta = split /\n/, $fasta;
      # compute description
      my $dupIds = $dbh->selectcol_arrayref("SELECT duplicate_id FROM SeqToDuplicate WHERE sequence_id = ?",
                                            {}, $subject);
      my @ids = ( $subject );
      push @ids, @$dupIds;
      my @subjectDescs = ();
      my $org;
      foreach my $id (@ids) {
        if ($id =~ m/^([a-zA-Z]+)::(.*)$/) {
          my ($db, $protId) = ($1,$2);
          my $info = $dbh->selectrow_hashref("SELECT * FROM CuratedGene WHERE db = ? AND protId = ?",
                                             {}, $db, $protId);
          die "Unknown curated item $db::$protId for $subject"
            unless defined $info;
          push @subjectDescs, $info->{desc};
          $org = $info->{organism} if $info->{organism} ne "";
        }
      }
      warning("No descriptions for $subject") if @subjectDescs == 0;
      $subjectDesc = join("; ", @subjectDescs);
      $subjectDesc .= " ($org)" if $org;
    }
    my $subject2 = $subject; $subject2 =~ s/:/_/g;
    # turn : in identifiers into _, for compatibility with newick format
    $fasta[0] = ">" . $subject2 . " " . $subjectDesc;
    $row->{fasta} = join("\n", @fasta);
  }
  my $headerLine = ">$id";
  $headerLine .= " $desc" if $desc ne "";
  my $fasta = join("\n", $headerLine, $seq, map $_->{fasta}, @keep)."\n";
  $fasta =~ s/\n+/\n/g; # remove blank lines (@fasta is not always chomped)
  my @lines = split /\n/, $fasta;
  my $seqsId = savedHash(\@lines, "seqs");

  print p(a({-href => "treeSites.cgi?seqsId=$seqsId" },
            "Build an alignment for", encode_entities($id), "and", scalar(@keep), "homologs")),
        p(start_form(-name => 'addSeqs', -method => 'POST', -action => 'treeSites.cgi'),
          hidden(-name => 'seqsId', -default => $seqsId, -override => 1),
          "Or add sequences from UniProt, PDB, RefSeq, or MicrobesOnline (separate identifiers with commas or spaces):",
          br(),
          textfield(-name => "addSeq", -default => "", -size => 50, -maxLength => 1000),
          br(),
          submit(-name => "Add"),
          end_form),
        p("Or",
          a({-href => findFile($seqsId, "seqs")}, "download"),
          "the sequences"),
        p("Or", a({-href => "treeSites.cgi"}, "start over")),
          end_html;
  exit(0);
} # end homologs mode

#else

my $seqsSet = param('seqsFile') || param('seqsId');
my $alnSet = param('alnFile') || param('alnId');
my $treeSet = param('treeFile') || param('treeId');

if ($seqsSet) {
  my $seqsId;
  my $fhSeqs;
  if (param('seqsId')) {
    $seqsId = param('seqsId');
    $fhSeqs = openFile($seqsId, 'seqs');
  } elsif (param('seqsFile')) {
    $fhSeqs = param('seqsFile')->handle;
  } else {
    die "Unreachable";
  }

  my %seqs = ();
  my %seqsDesc = ();
  my $state = {};

  while (my ($id, $seq) = ReadFastaEntry($fhSeqs, $state, 1)) {
    if ($id =~ m/^(\S+)\s(.*)/) {
      $id = $1;
      $seqsDesc{$id} = $2;
    }
    fail("Duplicate sequence for " . encode_entities($id))
      if exists $seqs{$id};
    fail(". or - is not allowed in unaligned sequences")
      if $seq =~ m/[.-]/;
    $seq =~ s/[*]//;
    fail("Illegal characters in sequence for " . encode_entities($id))
      unless $seq =~ m/^[a-zA-Z]+$/;
    fail("Sequence identifier " . encode_entities($id)
         . " contains an invalid character, one of  :(),;")
      if $id =~ m/[:(),;]/;
    fail("Sequences for " . encode_entities($id) . " is too long")
      if length($seq) > $maxLenAlign;
    $seqs{$id} = $seq;
  }
  fail(encode_entities($state->{error})) if defined $state->{error};
  fail("Too many sequences to align, the limit is $maxNAlign")
    if scalar(keys %seqs) > $maxNAlign;
  fail("Must have at least 2 sequences")
    if scalar(keys %seqs) < 2;

  # Handle addseq if set. If any do get set, set seqsId back to empty
  my $addSeq = param('addSeq');
  if (defined $addSeq) {
    $addSeq =~ s/^\s*//;
    $addSeq =~ s/\s*$//;
  }
  if (defined $addSeq && $addSeq ne "") {
    my %addSeq = ();
    my %addDesc = ();
    foreach my $id (split /[ \t,;]/, $addSeq) {
      next if $id eq "" || exists $addSeq{$id};
      my $id2 = $id; $id2 =~ s/:/_/g; # convert some ids to safe form
      if (exists $seqs{$id2}) {
        warning("Skipping", encode_entities($id), "which already has a sequence");
        next;
      }
      my ($def, $seq) = parseSequenceQuery(-query => $id,
                                           -dbh => $dbh,
                                           -blastdb => $blastdb,
                                           -fbdata => $fbdata);
      if (!defined $seq) {
        warning("Could not find sequence for", encode_entities($id));
      } else {
        print p(encode_entities($id), "is", encode_entities($def));
        $addSeq{$id2} = $seq;
        if ($def =~ m/ /) {
          $def =~ s/^\S+ +//;
          $addDesc{$id2} = $def;
        }
      }
    }
    # known sequences => 1
    my %bySeq = map { $seqs{$_} => $_ } keys %seqs;
    if (keys %addSeq > 0) { # adding sequences
      $seqsId = undef; # so it is rebuilt below
      my $formatdb = "../bin/blast/formatdb";
      die "No such executable: $formatdb" unless -x $formatdb;
      my $tmpDb = "$tmpPre.db";
      # Warn if the new sequence is not homologous to any of the previously existing sequences
      my $fh;
      open($fh, ">", $tmpDb) || die "Cannot write to $tmpDb";
      foreach my $id (sort keys %seqs) {
        print $fh ">$id\n" . $seqs{$id} . "\n";
      }
      close($fh) || die "Error writing to $tmpDb";
      my $formatCmd = "$formatdb -p T -i $tmpDb >& /dev/null";
      system($formatCmd) == 0 || die "Command failed: $formatCmd";

      while (my ($id,$seq) = each %addSeq) {
        my $seq = $addSeq{$id};
        warning("Warning:", encode_entities($id), "has the same sequence as", encode_entities($bySeq{$seq}))
          if exists $bySeq{$seq};
        my $seqFile = "$tmpPre.seq";
        open($fh, ">", $seqFile) || die "Cannot write to $seqFile";
        print $fh ">seq\n" . $seq . "\n";
        close($fh) || die "Error writing to $seqFile";
        my $hitsFile = "$tmpPre.hits";
        my $blastCmd = "$blastall -p blastp -e 0.01 -d $tmpDb -i $seqFile -o $hitsFile -m 8 >& /dev/null";
        system($blastCmd) ==0 || die "blast failed: $blastCmd";
        open($fh, "<", $hitsFile) || die "Cannot read $hitsFile from $blastCmd";
        my @hitLines = <$fh>;
        close($fh) || die "Error reading $hitsFile";
        unlink($seqFile);
        unlink($hitsFile);
        warning("Warning:", encode_entities($id), "is not similar to any of the initial sequences",
                "(no BLASTp hit with E < 0.01)")
          if @hitLines == 0;
        $seqs{$id} = $seq;
        $bySeq{$seq} = $id;
        $seqsDesc{$id} = $addDesc{$id} if exists $addDesc{$id};
      }
      unlink($tmpDb);
      foreach my $suffix (qw{phr pin psq}) {
        unlink("$tmpDb.$suffix");
      }
    } # end if adding sequences
  }

  if (!defined $seqsId) {
    my @lines = ();
    foreach my $key (sort keys %seqs) {
      my $header = ">$key";
      $header .= " $seqsDesc{$key}" if exists $seqsDesc{$key};
      push @lines, $header;
      push @lines, $seqs{$key};
    }
    $seqsId = savedHash(\@lines, "seqs");
  }

  if (param('buildAln')) {
    # Alignment building mode
    die "buildAln without seqsId" unless $seqsId;
    autoflush STDOUT 1; # show preliminary results
    my $muscle = "../bin/muscle3";
    die "No such executable: $muscle" unless -x $muscle;
    print p("Running",
            a({ -href => "https://doi.org/10.1093/nar/gkh340" }, "MUSCLE 3"),
            "on", scalar(keys %seqs), "sequences"), "\n";
    # also considered -diags but didn't necessarily speed things up much
    my $muscleOptions = "-maxiters 2 -maxmb 1000";
    print p("For speed, MUSCLE is running with $muscleOptions");
    my $tmpIn = "$tmpPre.in";
    open (my $fh, ">", $tmpIn) || die "Cannot write to $tmpIn";
    foreach my $id (sort keys %seqs) {
      print $fh ">", $id, "\n", $seqs{$id}, "\n";
    }
    close($fh) || die "Error writing to $tmpIn";
    my $tmpAln = "$tmpPre.aln";
    my $cmd = "$muscle -in $tmpIn $muscleOptions > $tmpAln";
    system($cmd) == 0 || die "$cmd\nfailed: $!";
    print p("MUSCLE succeeded"),"\n";

    open(my $fhAln, "<", $tmpAln) || die "Cannot read $tmpAln";
    my %aln = ();
    $state = {};
    while (my ($id, $seq) = ReadFastaEntry($fhAln, $state, 1)) {
      die "Unexpected sequence $id" unless exists $seqs{$id};
      $aln{$id} = $seq;
    }
    fail(encode_entities($state->{error})) if defined $state->{error};
    close($fhAln) || die "Error reading $tmpAln";
    unlink($tmpIn);
    unlink($tmpAln);
    my @lines = ();
    foreach my $id (sort keys %aln) {
      die "No aligned sequences for $id" unless exists $aln{$id};
      my $header = ">".$id;
      $header .= " " . $seqsDesc{$id} if exists $seqsDesc{$id};
      push @lines, $header, $aln{$id};
    }
    my $alnId = savedHash(\@lines, "aln");
    print p("Next",
            a({-href => "treeSites.cgi?alnId=$alnId"}, "build a tree").".");
    print p("Or download",
            a({-href => findFile($seqsId, "seqs") }, "sequences"),
            "or",
            a({-href => findFile($alnId, "aln") }, "alignment"));
  } else {
    # Show option to build alignment
    print
      p("Have", a({-href => findFile($seqsId,"seqs")}, scalar(keys %seqs), "sequences")),
      start_form(-name => 'input', -method => 'POST', -action => 'treeSites.cgi'),
      hidden(-name => 'seqsId', -default => $seqsId, -override => 1),
      p(submit(-name => "buildAln", -value => "Align with MUSCLE")),
      end_form;
  }
  print p("Or", a{-href => "treeSites.cgi"}, "start over");
  print end_html;
  exit(0);
}

# else
if (!$alnSet) {
  # Show the initial form to upload an alignment or sequences
  print
    p("View a phylogenetic tree along with selected sites from a protein alignment.",
      a({-href => "treeSites.cgi?alnId=DUF1080&treeId=DUF1080&tsvId=DUF1080&anchor=BT2157&pos=134,164,166",
         -title => "putative active site of the 3-ketoglycoside hydrolase family (formerly DUF1080)" },
        "See example.")),
    p("The first step is to search for characterized homologs of your sequence, or to upload your sequences."),
    start_form(-name => 'query', -method => 'POST', -action => 'treeSites.cgi'),
    p(b("Enter a protein sequence in FASTA or Uniprot format,",
        br(),
        "or an identifier from UniProt, RefSeq, or MicrobesOnline: ")),
    p({ -style => "margin-left: 2em;" },
      textarea( -name  => 'query', -value => '', -cols  => 70, -rows  => 10 )),
    p({ -style => "margin-left: 2em;" }, submit('Search'), reset()),
    end_form,
    start_form(-name => 'aln', -method => 'POST', -action => 'treeSites.cgi'),
    p(b("Or upload an alignment in fasta, clustal, or stockholm format."),
      "Limited to $maxN sequences or $maxMB megabytes."),
    p({-style => "margin-left: 2em;" }, filefield(-name => 'alnFile', -size => 50),
      br(),
      submit('Upload')),
    end_form,
    p(b("Or upload up to unaligned protein sequences in fasta format."),
      "Limited to $maxNAlign sequences."),
    start_form(-name => 'seqs', -method => 'POST', -action => 'treeSites.cgi'),
    p({ -style=> "margin-left: 2em;" }, filefield(-name => 'seqsFile', -size => 50),
      br(),
      submit('Upload')),
    end_form,
    p("Or try",
      a({-href => "sites.cgi" }, "SitesBLAST").":",
      "find homologs with known functional residues and see if they are conserved");
  print end_html;
  exit(0);
}

# else load the alignment
my @alnLines = ();
my $alnId;
if (param('alnFile')) {
  my $fhAln = param('alnFile')->handle;
  fail("alnFile not a file") unless $fhAln;
  @alnLines = <$fhAln>;
} elsif (param('alnId')) {
  $alnId = param('alnId');
  my $fhAln = openFile($alnId, "aln");
  @alnLines = <$fhAln>;
  close($fhAln) || die "Error reading $alnId.aln";
} else {
  fail("No alignment specified");
}

my %alnSeq; # with - as gaps (converted from "." if necessary, but potentially with lower-case)
my %alnDesc;
my %idInfo = (); # id to color or URL to value (from the input tsv)

my $hash;
if ($hash = ParseClustal(@alnLines)) {
  %alnSeq = %$hash;
} elsif ($hash = ParseStockholm(@alnLines)) {
  %alnSeq = %$hash;
} else {
  my $alnString = join("\n", @alnLines);
  open (my $fh, "<", \$alnString);
  my $state = {};
  while (my ($id, $seq) = ReadFastaEntry($fh, $state, 1)) {
    if ($id =~ m/^(\S+)\s(.*)/) {
      $id = $1;
      $alnDesc{$id} = $2;
    }
    fail("Duplicate sequence for " . encode_entities($id))
      if exists $alnSeq{$id};
    $alnSeq{$id} = $seq;
  }
  fail(encode_entities($state->{error})) if defined $state->{error};
}
fail("No sequences in the alignment")
  if (scalar(keys %alnSeq) == 0);

fail("Too many sequences in the alignment") if scalar(keys %alnSeq) > $maxN;

# Convert any . characters to -
while (my ($id, $seq) = each %alnSeq) {
  $seq =~ m/^[A-Za-z.-]+$/ || fail("Invalid sequence for " . encode_entities($id));
  $seq =~ s/[.]/-/g;
  $alnSeq{$id} = $seq;
}

my $alnLen;
while (my ($id, $seq) = each %alnSeq) {
  fail("Sequence identifier $id in the alignment contains an invalid character, one of  :(),;")
    if $id =~ m/[:(),;]/;
  $alnLen = length($seq) if !defined $alnLen;
  fail("Inconsistent sequence length for " . encode_entities($id))
    unless length($seq) == $alnLen;
}

# Save the alignment, if necessary
if (!defined $alnId) {
  $alnId = savedHash(\@alnLines, "aln");
}

my ($moTree, $treeId);

if (! $treeSet && param('buildTree')) {
  # Tree building mode
  my $trimGaps = param('trimGaps') ? 1 : 0;
  my $trimLower = param('trimLower') ? 1 : 0;
  my $ft = "../bin/FastTree";
  die "No such executable: $ft" unless -x $ft;

  # Trim the alignment
  my @keep = (); # positions to keep
  my $nSeq = scalar(keys %alnSeq);
  for (my $i = 0; $i < $alnLen; $i++) {
    my $nGaps = 0;
    my $nLower = 0;
    foreach my $seq (values %alnSeq) {
      my $char = substr($seq, $i, 1);
      if ($char eq '-') {
        $nGaps++;
      } elsif ($char eq lc($char)) {
        $nLower++;
      }
    }
    my $nUpper = $nSeq - $nGaps - $nLower;
    push @keep, $i
      unless ($trimGaps && $nGaps >= $nSeq/2)
        || ($trimLower && $nLower >= $nUpper);
  }

  print p("Removed positions that are at least 50% gaps.") if $trimGaps;
  print p("Removed positions that have as many lower-case as upper-case values.") if $trimGaps;
  print p("Trimmed to",scalar(@keep),"positions");

  if (scalar(@keep) < 10) {
    fail("Sorry: less than 10 alignment positions remained after trimming");
  }

  my $tmpTrim = "$tmpPre.trim";
  open (my $fhTrim, ">", $tmpTrim) || die "Cannot write to $tmpTrim";
  foreach my $id (sort keys %alnSeq) {
    my $seq = $alnSeq{$id};
    my $trimmed = join("", map substr($seq, $_, 1), @keep);
    print $fhTrim ">", $id, "\n", $trimmed, "\n";
  }
  close($fhTrim) || die "Error writing to $tmpTrim";

  my $tmpTreeFile = "$tmpPre.tree";
  autoflush STDOUT 1; # show preliminary results
  print p("Running",
          a({ -href => "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2835736/" },
            "FastTree 2"),
          "on an alignment with",
          scalar(@keep), "positions and $nSeq sequences.",
          "This could take a few minutes."), "\n";
  system("$ft -quiet $tmpTrim > $tmpTreeFile") == 0
    || die "FastTree failed on $tmpTrim : $!";

  $moTree = MOTree::new('file' => $tmpTreeFile)
    || die "Error parsing $tmpTreeFile";
  $moTree->rerootMidpoint();
  unlink($tmpTrim);
  unlink($tmpTreeFile);
  my $newick = $moTree->toNewick();
  $treeId = savedHash([$newick], "tree");
  print
    p("FastTree succeeded."),
    p("Rerooted the tree to minimize its depth (midpoint rooting)."),
    p(a{ -href => "treeSites.cgi?alnId=$alnId&treeId=$treeId"},
      "View the tree");
  print end_html;
  exit(0);
} # end tree building mode

die unless $alnId;

# Try to parse the table to get descriptions
my $tsvId;
if (param('tsvFile')) {
  my $fh = param('tsvFile')->handle;
  fail("tsvFile is not a file") unless $fh;
  my @lines = <$fh>;
  my $n = handleTsvLines(@lines);
  if ($n == 0) {
    warning("No descriptions found for matching ids in the uploaded table");
  } else {
    print p("Found $n descriptions in the uploaded table"),"\n";
    $tsvId = savedHash(\@lines, "tsv");
  }
} elsif (param('tsvId')) {
  $tsvId = param('tsvId');
  my $fh = openFile($tsvId, "tsv");
  my @lines = <$fh>;
  close($fh) || die "Error reading $tsvId.tsv";
  my $n = handleTsvLines(@lines);
  warning("No descriptions found for matching ids in the table") if $n == 0;
}

# Save alignment position mapping (id => 1-based position => 1-based aligned position)
my %seqPosToAlnPos = map { $_ => seqPosToAlnPos($alnSeq{$_}) } keys %alnSeq;
my %alnPosToSeqPos = (); # id => 1-based aligned non-gap position => 1-based sequence position
while (my ($id, $hash) = each %seqPosToAlnPos) {
  my $seq = $alnSeq{$id}; $seq =~ s/-//g;
  while (my ($seqPos, $alnPos) = each %$hash) {
    $alnPosToSeqPos{$id}{$alnPos} = $seqPos
      unless substr($seq, $seqPos-1, 1) eq "-";
  }
}

if (param('pattern')) {
  # pattern search mode
  my $pattern = uc(param('pattern'));
  fail("Invalid pattern: only amino acid characters or X or . are allowed")
    unless $pattern =~ m/^[.A-Z]+$/;
  print p("Searching for sequence matches to $pattern");
  $pattern =~ s/X/./g;
  my $anchorId = param('anchor');
  $anchorId = "" if !defined $anchorId;

  my %hits = (); # id to list of starting indexes
  my $nFound = 0;
  while (my ($id, $aln) = each %alnSeq) {
    my $seq = $aln; $seq =~ s/-//g;
    while ($seq =~ m/$pattern/gi) {
      # pos() returns the 0-based position past end of match; convert to 1-based beginning of match
      my $pos = pos($seq) - length($pattern) + 1;
      push @{ $hits{$id} }, $pos;
      $nFound++;
    }
  }
  my $nLimit = 100;
  if ($nFound == 0) {
    print p("Sorry, no matches found");
  } else {
    print p("Found $nFound matches in",
            scalar(keys %hits), "of", scalar(keys %alnSeq), "sequences.");
    if ($nFound > $nLimit) {
      print p("Showing the first $nLimit matches.");
    }

    my @ids = sort keys %hits;
    if ($anchorId ne "" && exists $hits{$anchorId}) {
      # push anchor sequence to front
      @ids = grep $_ ne $anchorId, @ids;
      unshift @ids, $anchorId;
    }

    # Show a table of matches.
    my @rows = ();
    push @rows, Tr(th("Protein"), th("Position"), th("Matching sequence"));

    my $nShow = 0;
    foreach my $id (@ids) {
      my $alnSeq = $alnSeq{$id};
      my $seqPosToAlnPos = $seqPosToAlnPos{$id};
      my $seq = $alnSeq; $seq =~ s/-//g;
      foreach my $pos (@{ $hits{$id} }) {
        if (++$nShow <= $nLimit) {
          my $idShow = $id;
          $idShow = a({ -title => $alnDesc{$id} }, $id)
            if exists $alnDesc{$id};
          my @matchShow = ();
          foreach my $offset (0..(length($pattern)-1)) {
            my $i = $pos + $offset;
            my $alignI = $seqPosToAlnPos->{$i-1} + 1;
            my $char = substr($seq, $i-1, 1);
            my $color = $colors{uc($char)} || "white";
            push @matchShow, a({-style => "background-color:$color; font-family:monospace; padding: 0.2em;", -title => "position $i ($alignI in alignment) is $char"}, $char);
          }
          push @rows, Tr(td($idShow,
                            td("$pos"."..".($pos+length($pattern)-1)),
                            td(join("",@matchShow))));
        }
      }
    }
    print table(@rows);
  } # end else has matches

  my $baseURL = "treeSites.cgi?alnId=$alnId";
  foreach my $attr (qw{treeId tsvId anchor pos zoom"}) {
    my $value = param($attr);
    $baseURL .= "&" . $attr . "=" . uri_escape($value)
      if defined $value && $value ne "";
  }
  print p(a({-href => $baseURL}, param('treeId') ? "Back to tree" : "Back to alignment"));
  print p("Or", a{-href => "treeSites.cgi"}, "start over");
  print end_html;
  exit(0);
} # end pattern search mode

if (! $treeSet) {
  print p("The alignment has", scalar(keys %alnSeq), "sequences of length $alnLen");

  # Form to compute a tree
  print
    start_form(-name => 'buildTree', -method => 'POST', -action => 'treeSites.cgi'),
    hidden(-name => 'alnId', -default => $alnId, -override => 1),
    p("Compute a tree:"),
    p(checkbox(-name => 'trimGaps', -checked => 1, -value => 1, -label => ''),
      "Trim columns that are &ge;50% gaps"),
    p(checkbox(-name => 'trimLower', -checked => 1, -value => 1, -label => ''),
      "Trim columns with more lowercase than uppercase"),
    p(submit(-name => "buildTree", -value => "Run FastTree")),
    end_form;

  print
    start_form(-name => 'inputTree', -method => 'POST', -action => 'treeSites.cgi'),
    hidden(-name => 'alnId', -default => $alnId, -override => 1),
    p("Or upload a rooted tree in newick format:", filefield(-name => 'treeFile', -size => 50)),
    p(submit(-name => "input", -value => "Upload")),
    end_form;
  print end_html;
  exit(0);
}

# else rendering mode

if (param('treeFile')) {
  my $fh = param('treeFile')->handle;
  fail("treeFile not a file") unless $fh;
  eval { $moTree = MOTree::new('fh' => $fh) };
} elsif (param('treeId')) {
  $treeId = param('treeId');
  my $fh = openFile($treeId, "tree");
  eval { $moTree = MOTree::new('fh' => $fh) };
  close($fh) || die "Error reading $treeId.tree";
} else {
  fail("No tree specified") unless defined $treeId;
}
fail("Could not parse tree") unless $moTree;

# Fail if there are any leaves in the tree without sequences
my $nodes = $moTree->depthfirst(); # depth-first ordering of all nodes
my $root = $moTree->get_root_node;
my @leaves = grep $moTree->is_Leaf($_), @$nodes;
my %idToLeaf = ();
foreach my $leaf (@leaves) {
  my $id = $moTree->id($leaf);
  $idToLeaf{$id} = $leaf;
  fail("Leaf named " . encode_entities($id) . " is not in the alignment")
    unless exists $alnSeq{$id};
}

# Issue a warning error for any sequences not in the tree, if this
# is the first time they were used together
if (!defined $alnId || !defined $treeId) {
  foreach my $id (keys %alnSeq) {
    warning("Sequence " . encode_entities($id) . " is not in the tree")
      unless exists $idToLeaf{$id};
  }
}

# Save the tree, if necessary
if (!defined $treeId) {
  my $newick = $moTree->toNewick();
  $treeId = savedHash([$newick], "tree");
}

# Finished loading input

my %branchLen = (); # node to branch length
my $missingLen = 0; # set if any leaf or internal node had no branch length

# Convert any negative branch lengths to 0
# Convert any missing branch lengths to 1
foreach my $node (@$nodes) {
  next if $node == $root;
  my $len = $moTree->branch_length($node);
  if ($len eq "") {
    print p($node);
    $missingLen = 1;
    $len = 1;
  }
  $len = 0 if $len < 0;
  $branchLen{$node} = $len;
}
warning("Missing branch lengths were set to 1") if $missingLen;

my $anchorId = param('anchor');
$anchorId = "" if !defined $anchorId;
fail("Unknown anchor id " . encode_entities($anchorId))
  if $anchorId ne "" && !exists $alnSeq{$anchorId};
my ($anchorAln, $anchorSeq, $anchorLen);
if ($anchorId ne "") {
  $anchorAln = $alnSeq{$anchorId};
  $anchorSeq = $anchorAln; $anchorSeq =~ s/[-]//g;
  $anchorLen = length($anchorSeq);
}

my @anchorPos;               # 1-based, and in the anchor if it is set
if (defined param('pos') && param('pos') ne "") {
  my $pos = param('pos');
  $pos =~ s/\s//g;
  my @posSpec = split /,/, $pos;
  foreach my $spec (@posSpec) {
    if ($spec =~ m/^(\d+):(\d+)$/ && $1 <= $2) {
      push @anchorPos, ($1..$2);
    } else {
      fail("Invalid position " . encode_entities($spec))
        unless $spec =~ m/^\d+$/;
      push @anchorPos, $spec;
    }
  }
  foreach my $i (@anchorPos) {
    fail("Invalid position $i") unless $i >= 1 && $i <= $alnLen;
    fail("position $i is past end of anchor " . encode_entities($anchorId))
      if $anchorId ne "" && $i > $anchorLen;
  }
}

my @alnPos = ();                # 0-based, and in the alignment
if ($anchorId eq "") {
  @alnPos = map { $_ - 1 } @anchorPos;
} else {
  my %anchorToAln = ();    # 1-based in anchor to 0-based in alignment
  my $at = 0;
  for (my $i = 0; $i < $alnLen; $i++) {
    my $c = substr($anchorAln, $i, 1);
    if ($c ne "-") {
      $at++;
      $anchorToAln{$at} = $i;
    }
  }
  @alnPos = map $anchorToAln{$_}, @anchorPos;
}

my $nodeZoom = param('zoom');
my @showLeaves = ();
if (defined $nodeZoom && $nodeZoom =~ m/^\d+$/) {
  die "Invalid node id $nodeZoom"
    if $nodeZoom == $root
      || !defined $moTree->ancestor($nodeZoom)
      || $moTree->is_Leaf($nodeZoom);
  @showLeaves = @{ $moTree->all_leaves_below($nodeZoom) };
  # update $nodes to only include this node and its descendents
  my @nodes = @{ $moTree->all_descendents($nodeZoom) };
  unshift @nodes, $nodeZoom;
  $nodes = \@nodes;
} else {
  $nodeZoom = undef;
  @showLeaves = @leaves;
}
my $rootUse = defined $nodeZoom ? $nodeZoom : $root;

my $baseURL = join("",
                   "treeSites.cgi?anchor=", uri_escape($anchorId),
                   "&pos=", join(",",@anchorPos),
                   "&treeId=", $treeId,
                   "&alnId=", $alnId,
                   "&tsvId=", $tsvId || "");

my @drawing;
if (defined $nodeZoom) {
  push @drawing, "Zoomed into a clade of " . scalar(@showLeaves) . " proteins, or see",
    a({ -href => $baseURL }, "all", scalar(@leaves), "proteins").".";
} else {
  push @drawing, "Drawing a tree for " . scalar(@leaves) . " proteins.";
}
if (@alnPos > 0 && $anchorId ne "") {
  push @drawing, "Position numbering is from " . encode_entities($anchorId) . ".";
}
print p(@drawing);
my @downloads = ();
push @downloads, a({ -href => findFile($treeId, "tree") }, "tree");
push @downloads, a({ -href => findFile($alnId, "aln") }, "alignment");
push @downloads, a({ -href => findFile($tsvId, "tsv") }, "table of descriptions")
  if $tsvId;

my $selfURL = $baseURL;
$selfURL .= "&zoom=$nodeZoom" if defined $nodeZoom;

print p("Download", join(" or ", @downloads),
        "or see", a({ -href =>  $selfURL }, "permanent link"),
        "to this page, or",
        a({ -href => "treeSites.cgi" }, "start over").".");

print p(start_form(-method => 'GET', -action => 'treeSites.cgi'),
        hidden( -name => 'alnId', -default => $alnId, -override => 1),
        hidden( -name => 'treeId', -default => $treeId, -override => 1),
        hidden( -name => 'tsvId', -default => $tsvId, -override => 1),
        hidden( -name => 'zoom', -default => defined $nodeZoom ? $nodeZoom : "", -override => 1),
        "Select positions",
        textfield(-name => "pos", -default => join(",",@anchorPos), -size => 30, -maxlength => 200),
        "in",
        textfield(-name => "anchor", -default => $anchorId, -size => 20, -maxlength => 200),
        submit(-value => "Go"),
        end_form);

print p(start_form(-method => 'POST', -action => 'treeSites.cgi'),
        hidden( -name => 'alnId', -default => $alnId, -override => 1),
        hidden( -name => 'treeId', -default => $treeId, -override => 1),
        hidden( -name => 'anchor', -default => $anchorId, -override => 1),
        hidden( -name => 'pos', -default => join(",",@anchorPos), -override => 1),
        hidden( -name => 'zoom', -default => defined $nodeZoom ? $nodeZoom : "", -override => 1),
        "Upload descriptions:", filefield(-name => 'tsvFile', -size => 50),
        submit(-value => "Go"),
        br(),
        small("The table should be tab-delimited with the sequence identifier in the 1st column",
              "and the description in the 2nd column. Optionally, add fields named",
              qq{<A HREF="https://www.december.com/html/spec/colorsvg.html">color</A>},
              "and URL."),
        end_form);

my $renderLarge = defined $nodeZoom || scalar(@leaves) <= 20;
my $renderSmall = !defined $nodeZoom && scalar(@leaves) > 100;
my @acts;
push @acts, "Hover or click on a leaf for information about that sequence." unless $renderLarge;
push @acts, "Click on an internal node to zoom in to that group." if $renderSmall;

print p(start_form( -onsubmit => "return leafSearch();" ),
        "Search for proteins to highlight:",
        textfield(-name => 'query', -id => 'query', -size => 20),
        button(-name => 'Search', -onClick => "leafSearch()"),
        button(-name => 'Clear', -onClick => "leafClear()"),
        br(),
        div({-style => "font-size: 80%; height: 1.5em;", -id => "searchStatement"}, ""),
        end_form),
  p(start_form,
    hidden( -name => 'alnId', -default => $alnId, -override => 1),
    hidden( -name => 'treeId', -default => $treeId, -override => 1),
    hidden( -name => 'tsvId', -default => $tsvId, -override => 1),
    hidden( -name => 'anchor', -default => $anchorId, -override => 1),
    hidden( -name => 'pos', -default => join(",",@anchorPos), -override => 1),
    hidden( -name => 'zoom', -default => defined $nodeZoom ? $nodeZoom : "", -override => 1),
    "Search for a sequence pattern:",
    textfield(-name => 'pattern', -size => 20),
    button(-name => 'Search'),
    end_form),

  p(@acts);

# Build an svg
# Layout:
# y axis (0 at top):
# 1 blank row at top, of height $padTop
# 1 row per leaf, of height $rowHeight
# space for the scale bar at bottom, of height padBottom
# x axis: $padLeft to $padLeft + $treeWidth has the tree
# spacer of width $pdMiddle
# then 1 column for each position
# and padRight
my $padTop = 70;
my $rowHeight = $renderSmall ? 3 : ($renderLarge ? 20 : 8);
my $minShowHeight = 20; # minimum height of a character to draw
my $padBottom = 70;
my $padLeft = 10;
my $treeWidth = 250;
my $padMiddle = 50;
my $padRight = $renderLarge ? 600 : 40; # space for labels
my $posWidth = 30;
my $svgHeight = $padTop + scalar(@showLeaves) * $rowHeight + $padBottom;
my $svgWidth = $padLeft + $treeWidth + $padMiddle + scalar(@alnPos) * $posWidth + $padRight;

my %rawY; # Unscaled y for each node (0 to nLeaves-1)
for (my $i = 0; $i < scalar(@showLeaves); $i++) {
  $rawY{ $showLeaves[$i] } = $i;
}
foreach my $node (reverse @$nodes) {
  if (!exists $rawY{$node}) {
    my @values = map $rawY{$_}, $moTree->children($node);
    die unless @values > 0;
    foreach my $value (@values) {
      die "rawY not set yet for child of $node" if !defined $value;
    }
    $rawY{$node} = sum(@values) / scalar(@values);
  }
}
my $maxY = max(values %rawY);
$maxY = 1 if $maxY == 0;
my %nodeY;
while (my ($node, $rawY) = each %rawY) {
  $nodeY{$node} = $padTop + $rowHeight * (0.5 + (scalar(@showLeaves)-1) * $rawY / $maxY);
}

my %rawX = ($rootUse => 0); # Unscaled y, with root at 0
foreach my $node (@$nodes) {
  next if $node eq $rootUse;
  my $parentX = $rawX{ $moTree->ancestor($node) };
  die $node unless defined $parentX;
  die $node unless defined $branchLen{$node};
  $rawX{$node} = $parentX + $branchLen{$node};
}
my %nodeX;
my $maxX = max(values %rawX);
$maxX = 0.5 if $maxX == 0;
while (my ($node, $rawX) = each %rawX) {
  $nodeX{$node} = $padLeft + $treeWidth * $rawX / $maxX;
}

my %leafHas = (); # all the shown positions for that leaf
foreach my $leaf (@showLeaves) {
  my $id = $moTree->id($leaf);
  my $seq = $alnSeq{$id} || die;
  my @val = map substr($seq, $_, 1), @alnPos;
  $leafHas{$leaf} = join("", @val);
}

my %nodeTitle = ();
foreach my $node (@$nodes) {
  my $id = $moTree->id($node);
  my $title = encode_entities($id);
  if ($moTree->is_Leaf($node)) {
    $title .= ": " . encode_entities($alnDesc{$id}) if exists $alnDesc{$id};
    $title .= " (has $leafHas{$node})" if @alnPos > 0;
  }
  $nodeTitle{$node} = $title;
}

my @svg = (); # lines in the svg

# For leaves, add an invisible horizontal bar with more opportunities for popup text
# These need to be output first to ensure they are behind everything else
for (my $i = 0; $i < @showLeaves; $i++) {
  my $leaf = $showLeaves[$i];
  my $x1 = 0; # $nodeX{$leaf} - 2;
  my $x2 = $svgWidth;
  my $width = $x2 - $x1;
  my $y1 = $nodeY{$leaf} - $rowHeight/2;
  push @svg, qq{<rect x="$x1" y="$y1" width="$x2" height="$rowHeight" fill="white" stroke="none" >};
  push @svg, qq{<TITLE>$nodeTitle{$leaf}</TITLE>\n</rect>};
}

# Show selected alignment positions (if any)
my $pos0X = $padLeft + $treeWidth + $padMiddle;
for (my $i = 0; $i < @alnPos; $i++) {
  my $pos = $alnPos[$i];
  my $left = $pos0X + $posWidth * $i;
  my $x = $left + $posWidth/2;
  my $labelY = $padTop - 3;
  my $labelChar = "#";
  $labelChar = substr($anchorAln, $pos, 1) if $anchorId ne "";
  my $colLabel = $labelChar . $anchorPos[$i];
  my $pos1 = $pos+1;
  my $title = "Alignment position $pos1";
  $title = "$anchorId has $labelChar at position $anchorPos[$i] (alignment position $pos1)" if $anchorId ne "";
  my $titleTag = "<TITLE>$title</TITLE>";
  push @svg, qq{<text text-anchor="left" transform="translate($x,$labelY) rotate(-45)">$titleTag$colLabel</text>};
  if (@showLeaves >= 20) {
    # show alignment position labels at bottom as well
    my $labelY2 = $svgHeight - $padBottom + 3;
    push @svg, qq{<text transform="translate($x,$labelY2) rotate(90)">${titleTag}$colLabel</text>"};
  }

  # draw boxes for every position
  foreach my $leaf (@showLeaves) {
    my $id = $moTree->id($leaf);
    my $seq = $alnSeq{$id} || die;
    my $char = uc(substr($seq, $pos, 1));
    my $top = $nodeY{$leaf} - $rowHeight/2;
    my $heightUse = $rowHeight + 0.2; # extra height to prevent thin white lines
    my $color = exists $colors{$char} ? $colors{$char} : "grey";
    my $encodedId = encode_entities($id);
    my $boxLeft = $left + $posWidth * 0.1;
    my $boxWidth = $posWidth * 0.8;
    push @svg, qq{<rect x="$boxLeft" y="$top" width="$boxWidth" height="$heightUse" style="fill:$color; stroke-width: 0;" >};
    my @val = map substr($seq, $_, 1), @alnPos;
    # Use spaces to emphasize this position, and show what position it is in that sequence
    $val[$i] .= ($alnPosToSeqPos{$id}{$pos}+1)
      if exists $alnPosToSeqPos{$id}{$pos};
    $val[$i] = " " . $val[$i] . " ";
    my $has = join("", @val);
    push @svg, qq{<TITLE>$encodedId has $has</TITLE>};
    push @svg, qq{</rect>};
  }

  # compute conservation of the position up the tree
  my %conservedAt = (); # node => value if it is conserved within this subtree, or ""
  foreach my $node (reverse @$nodes) {
    if ($moTree->is_Leaf($node)) {
      my $id = $moTree->id($node);
      die unless exists $alnSeq{$id};
      $conservedAt{$node} = substr($alnSeq{$id}, $pos, 1);
    } else {
      my @children = $moTree->children($node);
      my $char;
      foreach my $child (@children) {
        die unless exists $conservedAt{$child};
        if ($conservedAt{$child} eq "") {
          $char = "";
          last;
        } elsif (!defined $char) {
          $char = $conservedAt{$child};
        }  elsif ($char ne $conservedAt{$child}) {
          $char = "";
          last;
        }
      }
      $conservedAt{$node} = $char;
    }
  }

  # draw the character for each conserved clade, if there is space
  foreach my $node (@$nodes) {
    my $ancestor = $moTree->ancestor($node);
    if ($conservedAt{$node} ne "" && ($node == $rootUse || $conservedAt{$ancestor} eq "")) {
      # Check if the height of this subtree is at least $minShowHeight
      my @leavesBelow;
      if ($moTree->is_Leaf($node)) {
        @leavesBelow = ($node);
      } else {
        @leavesBelow = @{ $moTree->all_leaves_below($node) };
      }
      # Hover text to report how large the clade is and give an example id,
      # if the clade has more than one member; otherwise just
      # show the id and its sequence across selected positions
      my @leavesBelowY = map $nodeY{$_}, @leavesBelow;
      my $height = $rowHeight + max(@leavesBelowY) - min(@leavesBelowY);
      next unless $height >= $minShowHeight;
      my $midY = (max(@leavesBelowY) + min(@leavesBelowY))/2;
      my $title = "";
      my $leafUse = $leavesBelow[0];
      my $id = $moTree->id($leafUse);
      my $idShow = encode_entities($id);
      my $n1 = scalar(@leavesBelow) - 1;
      my $charTitle = $conservedAt{$leafUse};
      $charTitle .= ($alnPosToSeqPos{$id}{$pos}+1)
        if $n1 == 0 && exists $alnPosToSeqPos{$id}{$pos};
      $title = ($n1 > 0 ? "$idShow and $n1 similar proteins have" : "$idShow has") . " " . $charTitle
        . " aligning to " . $anchorPos[$i];
      push @svg, qq{<text text-anchor="middle" dominant-baseline="middle" x="$x" y="$midY"><TITLE>$title</TITLE>$conservedAt{$node}</text>};
    }
  }
} # End loop over positions

# Draw the tree after drawing the positions, so that text for leaf names (if displayed)
# goes on top of the color bars
foreach my $node (@$nodes) {
  next if $node == $rootUse;
  my $parent = $moTree->ancestor($node);
  die unless defined $parent;

  # draw lines left and then up or down to ancestor
  push @svg, qq{<line x1="$nodeX{$node}" y1="$nodeY{$node}" x2="$nodeX{$parent}" y2="$nodeY{$node}" stroke="black" />};
  push @svg, qq{<line x1="$nodeX{$parent}" y1="$nodeY{$node}" x2="$nodeX{$parent}" y2="$nodeY{$parent}" stroke="black" />};

  # draw node with popup info, if any
  # If it is a leaf, also make an (invisible) label; the group is to join these together
  my $radius;
  my $dotStyle = "";
  my $color = "";

  if ($moTree->is_Leaf($node) && $moTree->id($node) eq $anchorId) {
    $radius = $renderSmall ? 3 : 4;
  } else {
    # radius 2 is a bit small to click on but any bigger looks funny if $renderSmall
    $radius = $renderSmall ? 2 : 3;
  }

  if ($moTree->is_Leaf($node)) {
    my $id = $moTree->id($node);
    if (exists $idInfo{$id}{color}
        && $idInfo{$id}{color} ne ""
        && $idInfo{$id}{color} ne "black") {
      $radius++;
      $color = $idInfo{$id}{color};
      $color =~ s/[^a-zA-Z0-9#_-]//g; # remove problematic characters
      $color = "" unless $color =~ m/^[#a-zA-Z]/;
    }
    $color = "red" if $color eq "" && $id eq $anchorId;
    $dotStyle = qq{style="fill:$color;"} if $color;
  }

  push @svg, "<g>";
  my $onclick = $moTree->is_Leaf($node) ? qq{onclick="leafClick(this)"} : "";
  my $circle = qq{<circle cx="$nodeX{$node}" cy="$nodeY{$node}" r="$radius" $dotStyle $onclick >};
  $circle .= "<TITLE>$nodeTitle{$node}</TITLE>" if $nodeTitle{$node} ne "";
  $circle .= "</circle>";
  $circle = qq{<A xlink:href="${baseURL}&zoom=${node}">$circle</A>} if
    ! $moTree->is_Leaf($node) && $node ne $rootUse;
  push @svg, $circle;
  if ($moTree->is_Leaf($node)) {
    my $xLabel = $nodeX{$node} + $radius + 2;
    my $id = $moTree->id($node);
    my $idShow = encode_entities($id);
    $idShow = a({-href => $idInfo{$id}{URL}, -target => "_blank" }, $id) if $idInfo{$id}{URL};
    my $textStyle = "display:none; font-size:80%;";
    push @svg, qq{<text dominant-baseline="middle" x="$xLabel" y="$nodeY{$node}" text-anchor="left" style="$textStyle" >$idShow<TITLE>$nodeTitle{$node}</TITLE></text>};
  }
  push @svg, "</g>";
}

# Draw labels at right if $renderLarge
if ($renderLarge) {
  my $xLabel = $svgWidth - $padRight + 8;
  push @svg, qq{<g id="gLabels">};
  foreach my $node (@showLeaves) {
    my $id = $moTree->id($node);
    my $idShow = encode_entities($id);
    my $desc = "";
    $desc = encode_entities($alnDesc{$id})
      if exists $alnDesc{$id} && $alnDesc{$id} ne "";
    $desc .= " (has $leafHas{$node})" if @alnPos > 0;
    $idShow = qq{<tspan>$idShow</tspan><tspan style="font-size:80%;"> $desc</tspan>};
    $idShow = a({-href => $idInfo{$id}{URL}, -target => "_blank" }, $idShow) if $idInfo{$id}{URL};
    push @svg, qq{<text text-anchor="left" dominant-baseline="middle" x="$xLabel" y="$nodeY{$node}" ><title>$nodeTitle{$node}</title>$idShow</text>};
  }
  push @svg, "</g>";
}

# Scale bar
if (! $missingLen) {
  my @scales = reverse qw{0.001 0.002 0.005 0.01 0.02 0.05 0.1 0.2 0.5 1};
  while ($scales[0] > 0.8 * $maxX && @scales > 1) {
    shift @scales;
  }
  my $scaleSize = $scales[0];
  my $scaleLeft = $padLeft;
  my $scaleRight = $padLeft + $treeWidth * $scaleSize/$maxX;
  my $scaleY = $padTop + scalar(@showLeaves) * $rowHeight + $padBottom * 0.8;
  push @svg, qq{<line x1="$scaleLeft" y1="$scaleY" x2="$scaleRight" y2="$scaleY" stroke="black" />};
  my $scaleMid = ($scaleLeft+$scaleRight)/2;
  my $scaleY2 = $scaleY - 4;
  push @svg, qq{<text text-anchor="middle" x="$scaleMid" y="$scaleY2">$scaleSize /site</text>};
}

print join("\n",
  "<DIV>",
  qq{<SVG width="$svgWidth" height="$svgHeight" style="position: relative; left: 1em;">},
  qq{<g transform="scale(1)">},
  @svg,
  "</g>",
  "</SVG>",
  "</DIV>");

print end_html;
exit(0);

sub handleTsvLines {
  my @lines = @_;
  my ($iColor, $iURL);
  my $n = 0;
  if (@lines > 1) {
    my $header = $lines[0];
    $header =~ s/[\r\n]+//;
    my @fields = split /\t/, $header;
    ($iColor) = grep $fields[$_] eq "color", (0..(@fields-1));
    ($iURL) = grep $fields[$_] eq "URL", (0..(@fields-1));
  }

  foreach my $line (@lines) {
    $line =~ s/[\r\n]+$//;
    my @fields = split /\t/, $line;
    my ($id, $desc) = @fields;
    if (exists $alnSeq{$id} && defined $desc && $desc =~ m/\S/) {
      $alnDesc{$id} = $desc;
      $idInfo{$id}{color} = $fields[$iColor]
        if defined $iColor && defined $fields[$iColor];
      $idInfo{$id}{URL} = $fields[$iURL]
        if defined $iURL && defined $fields[$iURL];
      $n++;
    }
  }
  return $n;
}

sub findFile($$) {
  my ($id, $type) = @_;
  die "Undefined input to openFile()" unless defined $id && defined $type;
  die "Invalid type" unless $type=~ m/^[a-zA-Z_]+$/;
  die "Invalid id of type $type" unless $id =~ m/^[a-zA-Z0-9_-]+$/;
  my $file = "../static/$id.$type";
  $file = "$tmpDir/$id.$type" unless -e $file;
  die "No such file: $file" unless -e $file;
  return $file;
}

sub openFile($$) {
  my ($id, $type) = @_;
  my $file = findFile($id, $type);
  my $fh;
  open($fh, "<", $file) || die "Cannot read $file";
  return $fh;
}

sub savedHash($$) {
  my ($lines, $type) = @_;
  die unless defined $lines && defined $type;
  my $id = md5_hex(@$lines);
  my $file = "$tmpDir/$id.$type";

  if (! -e $file) {
    open(my $fh, ">", $file) || die "Cannot write to $file";
    foreach my $line (@$lines) {
      $line =~ s/[\r\n]+$//;
      print $fh $line."\n";
    }
    close($fh) || die "Error writing to $file";
  }
  return $id;
}

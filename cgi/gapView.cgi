#!/usr/bin/perl -w
# View pathway reconstructions

# Optional parameters:
# orgs -- either an orgId as in FetchAssembly or a multi-organism directory
# set -- which group of pathways
# orgId -- which genome
# path -- which pathway
# showdef -- show the pathway definition
# step -- which step in that pathway
# locusId -- which candidate gene for that step
# gdb -- which genome database to search in
# gquery -- an organism name to search for
# findgene -- a search term for finding proteins
# file -- uploaded fasta file (used to create orgs)
#	(This usually leads to running the analysis, or showing
#	the overview page)
# pathByScore -- sort pathways by score instead of by name
#	(in all pathways for this organism mode)
#
# Modes of this viewer:
# If gdb and gquery are set, search for relevant genomes
# Otherwise, if orgs is missing, either build it from gdb and gid,
#	or, show the front page
# Otherwise, require orgs (set defaults to "aa")
# If the analysis has not been run yet, it tries to run the analysis
#	If the orgs directory does not exist, it tries to fetch the genome first
# If the analysis was started recently, it waits
# Otherwise, the analysis exists:
# No other arguments -- list the genomes and pathways
#	(if there is just one organism, shows that overview instead)
# gaps -- list the gaps for all organisms or, if orgId is set, for 1 organism
# orgId -- overview of the organism
# path -- overview of the pathway across organisms
#	(if there is just one organism, shows that organism/pathway page instead)
# path & showdef -- detailed pathway definition (mostly, show the .steps file verbatim)
# orgId & path -- the pathway in the organism, with lists of rules and top candidates for each step
# orgId & path & step -- all candidates for the step, and the detailed definition of the step
# orgId & locusId -- show information about the gene, which step it is a candidate for, etc.
# orgId & path & step & locusId -- show relevant alignments
# orgId & findgene -- find genes matching

use strict;
use CGI qw(:standard Vars start_ul);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use HTML::Entities;
use IO::Handle qw{autoflush};
use DBI;
use lib "../lib";
use Steps;
use pbutils;
use pbweb qw{start_page GetMotd};
use FetchAssembly qw{CacheAssembly AASeqToAssembly GetMatchingAssemblies GetMaxNAssemblies};
use List::Util qw{max};
use File::stat;
use DB_File;
use URI::Escape;

# This uses the DB_File in $org.faa.db if that exists,
# or else tries to use fastacmd against the blast database.
# (Added DB_File because fastacmd fails for some uploaded fasta files.)
sub GetOrgSequence($$$$);
sub FaaToDb($$); # fasta file to db (from DB_File)

sub ProcessUpload($); # file handle to hashtable with gdb and gid, or with error

# handle step scores of 0, 1, 2, or empty
sub ScoreToStyle($);
sub ScoreToLabel($);
sub ShowScoreShort($);
sub RuleToMinScore($); # row in RuleScore to 0, 1, or 2
sub HMMToURL($);
sub GeneURL($$); # orgId (that is in %orgs), locusId
sub OrgToAssembly($); # orgId to assembly object
sub OrgIdToURL($);
sub ShowCandidatesForStep($); # row in StepScore => HTML for the top two candidates
sub LegendForColorCoding();
sub ShowWarnings($$); # orgId (or undef) and pathSpec (or undef)

sub Finish(); # show "About GapMind" and exit

# a a row from StepScore to a row from KnownGap (similar or same orgId), or undef.
# Also fills out orgId that the gap is known in.
# If this is from marker-based similarity (orgId != stepScore->orgId),
# then it also adds the identity and nMarkers fields.
# It ignores cases with different orgId if score is 1.
sub StepScoreToKnownGap($);

# Arguments: list of stepIds and a hash of stepId => row from StepScore
# Returns: formatted HTML for this path
sub PathToHTML($$);

# Arguments: stepId and row from StepScore
# Returns: short formatted HTML for this step
sub StepToShortHTML($$);

sub UniqueLoci(@);
sub SplitLoci(@);

# Given curatedIds (or uniprot: or curated2:), pathwayId, and stepId
# (if needed to link to the step definition),
# build a link
sub CuratedToLink($$$);
sub CuratedToSeq($); # fetch sequence

# Given a row from the Candidate table, format HTML for the description, %identity, and bit score
# for the other hit (if any). Returns a list of 3 HTML elements.
sub CandToOtherHTML($);

# Given a comment on a step or pathway or pathway instance, format HTML links.
sub LinkifyComment($);

sub DataForStepParts($$); # pathwayId & stepId to hashed information
sub FormatStepPart($$$); # the data, the step part row (in the database), and the orgId

sub RulesToHTML($$$); # the steps object, pathwayId, and orgId or ""

# Global variables
my $dataDir = "../tmp"; # where all the assemblies and GapMind results live
my %orgs = (); # orgId => hash including gdb, gid, genomeName
my $nCPU = $ENV{CPU_USE} || 12;
my $tmp = "/tmp/gapView.$$";

# maximum size of posted data, in bytes
my $maxMB = 100;
$CGI::POST_MAX = $maxMB*1024*1024;
my $maxNSeqsK = 100;
my $maxNSeqs = $maxNSeqsK * 1000;
my $charsInId = "a-zA-Z0-9:._-"; # only these characters are allowed in protein ids
# (Some other characters, like "|", break fastacmd, and
# ";" could creates issues with links because it affects CGI parameter parsing
# and the ids are currently not quoted)

# These are global because they are used by many subroutines to build URLs and HTML
my ($set, $orgsSpec, $setDesc);
# database handles for the curated database, the steps database, and the gaps database
my ($dbhC, $dbhS, $dbhG);

# Not-too-big data from the database
my %markerSim = (); # orgId to list of rows
my %knownGaps = (); # orgId => pathwayId => stepId => row (with orgId filled in)
my %stepDesc = (); # pathwayId => stepId => desc

my $transporterStyle = " background-color: gainsboro; padding:0.05em; border-radius: 0.25em;";

{
  FetchAssembly::SetFitnessBrowserPath("../fbrowse_data");
  FetchAssembly::SetPrivDir("../private");

  $set = param("set") || "aa";
  $set =~ m/^[a-zA-Z0-9._-]+$/ || die "Invalid set $set";
  my $stepsDir = "../tmp/path.$set"; # with curated.db and steps.db
  die "Invalid set $set: no $stepsDir directory" unless -d $stepsDir;
  die "set $set needs to be rebuilt"
    unless NewerThan("$stepsDir/steps.db", "$stepsDir/curated.db");
  $dbhC = DBI->connect("dbi:SQLite:dbname=${stepsDir}/curated.db","","",{ RaiseError => 1 }) || die $DBI::errstr;
  $dbhS = DBI->connect("dbi:SQLite:dbname=${stepsDir}/steps.db","","",{ RaiseError => 1 }) || die $DBI::errstr;
  my ($stepsVersion) = $dbhS->selectrow_array("SELECT stepsVersion FROM Version");

  ($setDesc) = $dbhS->selectrow_array("SELECT desc FROM Pathway WHERE pathwayId = 'all'");
  die "No pathway named all" unless defined $setDesc;
  my $pathInfo = $dbhS->selectall_arrayref("SELECT * from Pathway WHERE pathwayId <> 'all'", { Slice => {} });
  my %pathDesc = map { $_->{pathwayId} => $_->{desc} } @$pathInfo;

  autoflush STDOUT 1; # show preliminary results
  my $banner = "GapMind for $setDesc";

  my @gdbs = ("NCBI", "IMG", "UniProt", "MicrobesOnline", "FitnessBrowser");
  my %gdb_labels1 = ("NCBI" => "NCBI assemblies",
                     "UniProt" => "UniProt proteomes",
                     "IMG" => "JGI/IMG genomes", "FitnessBrowser" => "Fitness Browser genomes",
                     "local" => "Uploaded proteome");
  my %gdb_labels = map { $_ => exists $gdb_labels1{$_} ? $gdb_labels1{$_} : "$_ genomes"} @gdbs;

  $orgsSpec = param('orgs');
  $orgsSpec = param('gdb') . "__" . param('gid')
    if !defined $orgsSpec && param('gdb') && param('gid');

  if (defined param('file')) {
    # Process the upload file and set $orgsSpec
    my $upFile = param('file');
    my $error;
    my %up;
    if (ref $upFile) {
      my $fhUp = $upFile->handle || die "Not a file handle";
      %up = ProcessUpload($fhUp);
      $error = $up{error};
    } else {
      $error = "No upload file specified";
    }
    if ($error) {
      start_page('title' => 'Upload Error in GapMind',
                 'banner' => $banner,
                 'bannerURL' => "gapView.cgi?set=${set}");
      print p(HTML::Entities::encode($error));
      Finish();
    }
    # else
    $orgsSpec = $up{gdb} . "__" . $up{gid};
  }

  if (!defined $orgsSpec && param('gquery')) {
    # mode: Find genome
    my $gquery = param('gquery');
    my $gdb = param('gdb') || die "Must specify gdb with gquery";
    die "Unknown genome database: $gdb\n"
      if !exists $gdb_labels{$gdb};
    start_page('title' => "GapMind for $setDesc",
               'banner' => $banner,
               'bannerURL' => "gapView.cgi?set=$set");
    print p("Searching", $gdb_labels{$gdb}, "for", "'" . HTML::Entities::encode($gquery) . "'"), "\n";
    my @rows = GetMatchingAssemblies($gdb, $gquery);
    my $limit = GetMaxNAssemblies();
    if (@rows > 0) {
      my $desc = "Found " . scalar(@rows) . " assemblies";
      if (@rows == 1) {
        $desc = "Found 1 assembly";
      } elsif (@rows >= $limit) {
        $desc = "Found the first " . scalar(@rows) . " matching assemblies";
      }
      $desc .= ", please choose one" if @rows > 1;
      print p($desc . ":"), "\n";
      print start_form(-method => 'get', -action => 'gapView.cgi'),
        hidden(-name => 'gdb', -value => $gdb, -override => 1),
        hidden(-name => 'set', -value => $set, -override => 1);
      foreach my $row (@rows) {
        my $checked = @rows == 1 ? "CHECKED" : "";
        print qq{<INPUT type="radio" NAME="gid" VALUE="$row->{gid}" $checked />},
          a({ -href => $row->{URL} }, $row->{genomeName} ),
            " ", small("(" . $row->{gid} . ")"), br();
      }
      print p(submit(-name => "Analyze")),
            end_form;
      Finish();
    } else {
      print p("Sorry, no matching genomes were found.");
    }
    print p("Try", a({-href => "gapView.cgi?set=$set&gdb=$gdb"}, "another genome"));
    Finish();
  } elsif (!defined $orgsSpec) {
    # mode: Front page
    start_page('title' => "GapMind: Automated annotation of $setDesc",
               'banner' => $banner,
               'bannerURL' => "gapView.cgi?set=${set}");
    my $exampleURL = "gapView.cgi?set=${set}&orgs=orgsFit";
    $exampleURL .= "&gaps=1" if $set eq "aa";
    my $exampleTitle = $set eq "aa" ? "View 'gaps' in" : "View results for";
    print
      GetMotd(),
      p($exampleTitle,
        a({ -href => $exampleURL }, "35 bacteria that grow in minimal media,"),
        "or choose a genome to analyze:"),
      start_form(-method => 'get', -action => "gapView.cgi", -autocomplete => 'on'),
      hidden(-name => 'set', -value => $set, -override => 1),
      p("Genome database to search:",
        popup_menu(-name => 'gdb', -values => \@gdbs, -labels => \%gdb_labels, -default => $gdbs[0])),
      p(textfield(-name => 'gquery', -value => '', -size => 50, -maxlength => 200)),
      p(small("Example:", a({-href => "gapView.cgi?set=${set}&gdb=NCBI&gquery=Desulfovibrio vulgaris"}, "Desulfovibrio vulgaris"))),
      p(submit(-name => "findgenome", -value => 'Find Genome')),
      end_form,
      start_form(-method => 'post', -action => "gapView.cgi", -autocomplete => 'on'),
      hidden(-name => 'set', -value => $set, -override => 1),
      p("Or upload a proteome in fasta format:",
        filefield(-name=>'file', -size=>40), br(), submit('Upload')),
      end_form;
    Finish();
  }

  $orgsSpec =~ m/^[a-zA-Z0-9._-]+$/ || die "Invalid orgs $orgsSpec";
  my $orgPre = "../tmp/$orgsSpec/orgs";
  my $sumPre = "../tmp/$orgsSpec/$set.sum";

  my $alreadyBuilt = NewerThan("$sumPre.db", "$stepsDir/steps.db");

  my $orgSetsFile = "../tmp/path.$set/orgSets.tsv";
  if (! $alreadyBuilt && -e $orgSetsFile
      && ! param("orgId")
      && param("gdb") && param("gid")) {
    # try to find this organism in a standard set
    my @orgsInSets = ReadTable($orgSetsFile, qw{orgId gdb gid orgSet});
    my $orgId = join("__", param("gdb"), param("gid"));
    @orgsInSets = grep { $_->{orgId} eq $orgId } @orgsInSets;
    if (@orgsInSets > 0) {
      my $orgSet = $orgsInSets[0]{orgSet};
      my $orgPre2 = "../tmp/$orgSet/orgs";
      my $sumPre2 = "../tmp/$orgSet/$set.sum";
      if (-e "$orgPre2.org" && NewerThan("$sumPre2.db", "$stepsDir/steps.db")) {
        # Use this set instead
        $orgsSpec = $orgSet;
        $orgPre = $orgPre2;
        $sumPre = $sumPre2;
        $alreadyBuilt = 1;
        param("orgId", $orgId); # orgId is initialized below
      }
    }
  }

  if ($alreadyBuilt) {
    # Verify that the gaps database is based on the correct version
    $dbhG = DBI->connect("dbi:SQLite:dbname=${sumPre}.db","","",{ RaiseError => 1 }) || die $DBI::errstr;
    my ($gapsVersion) = $dbhG->selectrow_array("SELECT stepsVersion FROM Version");
    if ($gapsVersion ne $stepsVersion) {
      $alreadyBuilt = 0;
      $dbhG->disconnect();
    }
  }

  # Wait up to 5 minutes for a previously running job to finish
  if (! $alreadyBuilt
      && -e "$sumPre.begin"
      && stat("$sumPre.begin")->mtime >= time() - 5*60
      && !param('force')) {
    # mode: The analysis is already running
    start_page('title' => 'Analysis in progress',
               'banner' => $banner,
               'bannerURL' => "gapView.cgi?set=${set}");
    print
      p("Analysis of $setDesc is already underway. Please check",
        a({-href => "gapView.cgi?orgs=$orgsSpec&set=$set"}, "here"),
        "in a few seconds");
      Finish();
  }

  unless ($alreadyBuilt) {
    # mode: Run the analysis
    start_page('title' => "Analyzing $setDesc",
               'banner' => $banner,
               'bannerURL' => "gapView.cgi?set=${set}");
    print "\n";
    unless (-e "$orgPre.org") {
      # Try to load the organism
      $orgsSpec =~ m/^([^_]+)__(.*)$/ || die "Invalid organism specifier $orgsSpec";
      my ($gdb,$gid) = ($1,$2);
      print p("Fetching assembly $gid from $gdb"), "\n";
      CacheAssembly($gdb, $gid, $dataDir) || die;
      mkdir("$dataDir/$orgsSpec");
      # buildorgs.pl creates $orgPre.org and $orgPre.faa and runs formatdb
      my @cmd = ("../bin/buildorgs.pl", "-out", $orgPre, "-orgs", $gdb.":".$gid);
      system(@cmd) == 0
        || die "command @cmd\nfailed with error: $!";
      FaaToDb("$orgPre.faa", "$orgPre.faa.db");
    }
    my @orgs = ReadOrgTable("$orgPre.org");
    die "No organisms for $orgPre.org" unless @orgs > 0;
    die "No such file: $orgPre.faa" unless -e "$orgPre.faa";
    system("touch", "$sumPre.begin");
    my $timeStart = time();
    my $nSec = scalar(@orgs) * ($set eq "aa" ? 12 : 30);
    print p("Analyzing $setDesc in", scalar(@orgs), "genomes.",
            "This should take around", $nSec, "seconds.")."\n";
    my @cmds = ();
    push @cmds, ["../bin/gapsearch.pl", "-orgs", $orgPre, "-set", $set,
                 "-out", "$dataDir/$orgsSpec/$set.hits",
                 "-nCPU", $nCPU ];
    push @cmds, ["../bin/gaprevsearch.pl", "-orgs", $orgPre,
                 "-hits", "$dataDir/$orgsSpec/$set.hits",
                 "-curated", "$stepsDir/curated.faa.udb",
                 "-out", "$dataDir/$orgsSpec/$set.revhits",
                 "-nCPU", $nCPU ];
    push @cmds, ["../bin/gapsummary.pl", "-orgs", $orgPre, "-set", $set,
                 "-hits", "$dataDir/$orgsSpec/$set.hits",
                 "-rev", "$dataDir/$orgsSpec/$set.revhits",
                 "-out", "$dataDir/$orgsSpec/$set.sum" ];
    my $gapReqFile = "$sumPre.warn";
    push @cmds, ["../bin/checkGapRequirements.pl", "-org", "../tmp/$orgsSpec", "-set", $set,
                 "-out", $gapReqFile];
    my $markerFaa = "../gaps/$set/$set.known.gaps.markers.faa";
    my $knownSimFile = "$sumPre.knownsim";
    if (-e $markerFaa) {
      push @cmds, ["../bin/orgsVsMarkers.pl", "-orgs", $orgPre,
                   "-vs", "../gaps/$set/$set.known.gaps.markers.faa",
                   "-out", $knownSimFile ]
    } else {
      unlink($knownSimFile);
    }
    my $buildCmd = ["../bin/buildGapsDb.pl", "-gaps", $sumPre,
                    "-requirements", $gapReqFile,
                    "-steps", "$stepsDir/steps.db",
                    "-out", "$sumPre.db"];
    push @$buildCmd, ("-markersim", $knownSimFile) if -e $markerFaa;
    push @cmds, $buildCmd;
    my %label = ('gapsearch.pl' => 'Searching for candidates for each step',
                 'gaprevsearch.pl' => 'Comparing candidates to other curated proteins',
                 'gapsummary.pl' => 'Scoring each candidate and pathway',
                 'orgsVsMarkers.pl' => 'Comparing to marker genes for known gaps');
    foreach my $cmd (@cmds) {
      my $show = $cmd->[0]; $show =~ s!.*/!!;
      print p($label{$show})."\n" if exists $label{$show};
      system(@$cmd) == 0 || die "Command failed\n@$cmd\nError code: $!";
    }
    unlink("$sumPre.begin");
    my $timeDiff = time() - $timeStart;
    print "</pre>\n",
      p("Analysis succeeded in $timeDiff seconds, please",
      a({-href => "gapView.cgi?orgs=$orgsSpec&set=$set"}, "view results")),
      end_html;
    exit(0);
  }

  #else -- viewing modes
  my @orgs = ReadOrgTable("$orgPre.org");
  die "No organisms for $orgPre.org" unless @orgs > 0;
  %orgs = map { $_->{orgId} => $_ } @orgs;
  my $faafile = "$orgPre.faa";
  die "No such file: $faafile" unless -e $faafile;

  my $pathSpec = param("path");
  $pathSpec = "" if !defined $pathSpec;
  die "Invalid path parameter $pathSpec"
    unless $pathSpec eq "" || $pathSpec =~ m/^[a-zA-Z0-9._-]+$/;
  my @path = (); # all pathways, or the one specified

  if ($pathSpec eq "") {
    @path = map { $_->{pathwayId} } @$pathInfo;
  } else {
    die "Unknown pathway $pathSpec" unless exists $pathDesc{$pathSpec};
    die "Invalid path $pathSpec" unless $pathSpec =~ m/^[a-zA-Z0-9._'-]+$/;
  }

  my $orgId = param("orgId"); # the organism specified, or, ""
  $orgId = "" if !defined $orgId;
  if ($orgId ne "") {
    die "Unknown orgId $orgId" unless exists $orgs{$orgId};
  }
  $orgId = $orgs[0]{orgId} if @orgs == 1;

  my $stepSpec = param("step");
  $stepSpec = "" if !defined $stepSpec;
  if ($stepSpec ne "") {
    die "step specified without path specified" unless $pathSpec ne "";
    $stepSpec =~ m/^[a-zA-Z0-9._'"-]+$/ || die "Invalid step $stepSpec";
    my $st = $dbhS->selectall_arrayref("SELECT * FROM Step WHERE pathwayId = ? AND stepId = ?",
                                       {}, $pathSpec, $stepSpec);
    die "Non-existent step $stepSpec" unless @$st == 1;
  }

  my $locusSpec = param("locusId");
  $locusSpec = "" if !defined $locusSpec;
  $locusSpec =~ m/^[$charsInId]*$/ || die "Invalid locus $locusSpec";

  my $findGene = param('findgene');
  $findGene = "" if !defined $findGene;
  $findGene =~ s/^[ \t]+//;
  $findGene =~ s/[ \t]+$//;
  my $findGeneShow = HTML::Entities::encode($findGene);

  my $title = $setDesc;
  $title = "Potential Gaps in $setDesc" if param('gaps');
  if ($locusSpec ne "") {
    $title = $stepSpec ne "" ? "Aligments for a candidate for $stepSpec" : "Protein $locusSpec";
  } elsif ($pathSpec ne "") {
    $title = $pathDesc{$pathSpec};
  }
  $title = "Finding step $stepSpec for $pathDesc{$pathSpec}"
    if $stepSpec ne "" && $orgId ne "" && $pathSpec ne "" && $locusSpec eq "";
  $title = "Searching for proteins"
    if $orgId ne "" && $findGene ne "";
  $title .= " in $orgs{$orgId}{genomeName}"
    if $orgId ne "";
  $title = "Definition of $pathDesc{$pathSpec}" if $pathSpec ne "" && param("showdef");
  my $nOrgs = scalar(@orgs);
  start_page('title' => $title,
             'banner' => $banner,
             'bannerURL' => "gapView.cgi?set=${set}");
  print "\n";

  # Load from the database:
  # step descriptions
  my $stepDesc = $dbhS->selectall_arrayref("SELECT pathwayId,stepId,desc FROM Step");
  foreach my $row (@$stepDesc) {
    my ($pathwayId, $stepId, $desc) = @$row;
    $stepDesc{$pathwayId}{$stepId} = $desc;
  }

  # marker comparisons
  my $markerSim = $dbhG->selectall_arrayref("SELECT * FROM MarkerSimilarity", { Slice => {} });
  foreach my $row (@$markerSim) {
    push @{ $markerSim{ $row->{orgId} } }, $row;
  }

  # known gaps
  my $knownGaps = $dbhS->selectall_arrayref("SELECT * from KnownGap", { Slice => {} });
  foreach my $row (@$knownGaps) {
    my $orgId = $row->{gdb} . "__" . $row->{gid};
    $row->{orgId} = $orgId; # so we don't need to compute this again
    $knownGaps{$orgId}{ $row->{pathwayId} }{ $row->{stepId} } = $row;
  }

  my @orgsSorted = sort { $a->{genomeName} cmp $b->{genomeName} } @orgs;
  my @ruleScoreLabels = ("has a gap", "may have a gap", "all steps were found");

  my @links = ();     # a list of items to put inside li at the bottom
  if ($pathSpec ne "" && param("showdef") && param("showdef") eq "literal") {
    # mode: Show the definition of this pathway, as literal raw contents
    print p("As text, or see",
            a({ -href => "gapView.cgi?orgs=$orgsSpec&set=$set&path=$pathSpec&showdef=1" },
              "rules and steps"));
    my $stepPath = "../gaps/$set";
    my $stfile = "$stepPath/$pathSpec.steps";
    open (my $fh, "<", $stfile) || die "No such file: $stfile\n";
    my @lines = <$fh>;
    close($fh) || die "Error reading $fh";
    print pre(join("",@lines)), "\n";
    push @links, a({-href => "$stepsDir/$pathSpec.query"}, "Table of queries for $pathSpec")
      . " (tab-delimited)";
  } elsif ($pathSpec ne "" && param("showdef")) {
    # mode: Show the definition of this pathway, formatted
    # Some information is in the steps file only, not in the database, including
    # comments on steps and the preferred ordering of steps is only in the steps file.
    my $stepsObj = Steps::ReadSteps("../gaps/$set/$pathSpec.steps");
    print p("As rules and steps, or see",
            a({ -href => "gapView.cgi?set=$set&orgs=$orgsSpec&path=$pathSpec&showdef=literal" },
            "full text"));
    print h3("Rules");
    print p("Overview:", LinkifyComment($stepsObj->{topComment})) if $stepsObj->{topComment} ne "";
    print RulesToHTML($stepsObj, $pathSpec, ""), # no orgId
      "\n";
    print h3("Steps");
    my @stepsInOrder = sort { $a->{i} <=> $b->{i} } values %{ $stepsObj->{steps} };
    foreach my $stepObj (@stepsInOrder) {
      print p(a({ -name => $stepObj->{name} }, b(i($stepObj->{name}) . ": " . $stepObj->{desc})));
      print start_ul;
      my $data = DataForStepParts($pathSpec, $stepObj->{name});
      my $stepParts = $dbhS->selectall_arrayref("SELECT * from StepPart WHERE pathwayId = ? AND stepId = ?",
                                                { Slice => {} }, $pathSpec, $stepObj->{name});
      foreach my $stepPart (@$stepParts) {
        print li(FormatStepPart($data, $stepPart, ""));
      }
      print li("Comment:", LinkifyComment($stepObj->{comment}))
        if $stepObj->{comment} ne "";
      print end_ul, "\n";
    }
    push @links, a({-href => "$stepsDir/$pathSpec.query"}, "Table of queries for $pathSpec")
      . " (tab-delimited)";
  } elsif (param('gaps')) {
    # mode: Overview of gaps, either for 1 organism or all organisms
    my $gaps;
    if ($orgId eq "") {
      $gaps = $dbhG->selectall_arrayref("SELECT * from StepScore WHERE score != 2 AND onBestPath=1",
                                       { Slice => {} });
    } else {
      $gaps = $dbhG->selectall_arrayref("SELECT * from StepScore WHERE score != 2 AND onBestPath=1 AND orgId = ?",
                                        { Slice => {} }, $orgId);
    }
    my @sorted = sort { lc($a->{pathwayId}) cmp lc($b->{pathwayId})
                          || lc($a->{stepId}) cmp lc($b->{stepId})
                          || $orgs{ $a->{orgId} }{genomeName} cmp $orgs{ $b->{orgId} }{genomeName} 
                      } @$gaps;
    if (@sorted == 0) {
      print p("Each pathway in",
              a({ -href => "gapView.cgi?orgs=$orgsSpec&set=$set&orgId=$orgId" },
                $orgs{$orgId}{genomeName}),
              "has a high-confidence path");
    } else {
      my $nLo = scalar(grep {$_->{score} eq "0" || $_->{score} eq ""} @sorted);
      my $nMed = scalar(grep {$_->{score} eq "1"} @sorted);
      my $nTot = scalar(@sorted);
      die unless $nTot == $nLo + $nMed;
      my $totals = "Found $nLo low-confidence and $nMed medium-confidence steps on the best paths for "
        . scalar(@$pathInfo) . " pathways";
      $totals .= " x " . scalar(@orgs) . " genomes" if $orgId eq "" && @orgs > 1;
      my $nCurated = 0;
      my $nKnown = 0;
      my $nThis = 0;
      foreach  my $gap (@sorted) {
        my $knownGap = StepScoreToKnownGap($gap);
        if (defined $knownGap) {
          my $knownOrgId = $knownGap->{orgId};
          if ($knownOrgId eq $gap->{orgId} && $knownGap->{gapClass} ne "") {
            $nCurated++;
          } elsif ($knownOrgId eq $gap->{orgId}) {
            $nThis++;
          } else {
            $nKnown++;
          }
        }
      }
      $totals .= ".";
      $totals .= " $nCurated of $nTot gaps have been manually classified."
        if $nCurated > 0;
      my $stringThis = $nThis > 0 ? "this or" : "";
      $totals .= " " . ($nThis + $nKnown) . " of $nLo low-confidence steps are known gaps in $stringThis related organisms."
        if $nKnown > 0;
      print p($totals);

      my @th = qw{Pathway Step Organism Best-candidate 2nd-candidate};
      my $showOrg = $orgId eq "" && scalar(@orgs) > 1;
      @th = grep { $_ ne "Organism" } @th unless $showOrg;
      map s/-/ /g, @th;
      push @th, "Class of gap" if $nCurated > 0;
      push @th, "Known?" if $nKnown > 0;
      my @tr = ();
      push @tr, Tr(th({-valign => "top"}, \@th));
      foreach my $row (@sorted) {
        my ($show1, $show2) = ShowCandidatesForStep($row);
        my $p = $row->{pathwayId};
        my $s = $row->{stepId};
        my $o = $row->{orgId};
        my @td = ( a({ -href => "gapView.cgi?orgs=$orgsSpec&set=$set&orgId=$o&path=$p",
                       -title => $pathDesc{$p},
                       -style => ScoreToStyle($row->{score}) }, $p),
                   a({ -href => "gapView.cgi?orgs=$orgsSpec&set=$set&orgId=$o&path=$p&step=$s",
                       -title => ScoreToLabel($row->{score}),
                       -style => ScoreToStyle($row->{score}) },
                     "$s: $stepDesc{$p}{$s}") );
        push @td, a({ -href => "gapView.cgi?orgs=$orgsSpec&set=$set&orgId=$o" }, $orgs{$o}{genomeName})
          if $showOrg;
        push @td, ( $show1, $show2 );
        my $known = StepScoreToKnownGap($row);
        if ($nCurated > 0) {
          my $curatedComment = "&nbsp;";
          $curatedComment = a({ -title => $known->{comment}, -style => "color: darkgreen;" },
                              $known->{gapClass})
            if defined $known && $known->{gapClass} ne "";
          push @td, $curatedComment;
        }
        if ($nKnown) {
          my $knownDesc = "&nbsp;";
          if (defined $known) {
            my $title = "Despite the apparent lack of $s,"
              . " $known->{genomeName} performs $pathDesc{$p}.";
            my $text = "known gap";
            if ($known->{orgId} ne $o) {
              my $idShow = int(0.5 + $known->{identity}) . "%";
              $title .= " Across $known->{nMarkers} ribosomal proteins, it is"
                            . " $idShow identical to $orgs{$o}{genomeName}.";
              $text .= " ($idShow id.)";
            }
            $knownDesc = a({ -title => $title,
                             -href => "gapView.cgi?gid=$known->{gid}&gdb=$known->{gdb}&set=$set",
                             -style => "color: black;" },
                           $text);
          }
          push @td, $knownDesc;
        }
        push @tr, Tr(td({-valign => "top" }, \@td));
      }
      print table({-cellpadding=>2, -cellspacing=>0, -border=>1}, @tr), "\n";
      print LegendForColorCoding();
      ShowWarnings($orgId, undef); # all pathways
    } # end else #gaps > 0
  } elsif ($orgId eq "" && $pathSpec eq "") {
    # mode: List pathways & genomes
    print p(scalar(@path), "pathways");
    print start_ul;
    foreach my $path (@path) {
      print li(a({ -href => "gapView.cgi?orgs=$orgsSpec&set=$set&path=$path" },
                 $pathDesc{$path}));
    }
    print end_ul, "\n";
    print p(scalar(@orgsSorted), "genomes"), "\n";
    print start_ul;
    foreach my $org (@orgsSorted) {
      my $orgId = $org->{orgId};
      my $URL = "gapView.cgi?orgs=$orgsSpec&set=$set&orgId=$orgId";
      print li(a({ -href => $URL }, $org->{genomeName} ));
    }
    print end_ul, "\n";
  } elsif ($orgId eq "" && $pathSpec ne "") {
    # mode: Overview of this pathway across organisms
    print p("Analysis of pathway $pathSpec in", scalar(@orgs), "genomes"), "\n";
    my $ruleScores = $dbhG->selectall_arrayref("SELECT * from RuleScore WHERE pathwayId = ? AND ruleId = 'all'",
                                               { Slice => {} }, $pathSpec);
    my %orgScore = map { $_->{orgId} => $_ } @$ruleScores;
    my $stepScores = $dbhG->selectall_arrayref("SELECT * from StepScore WHERE pathwayId = ?",
                                               { Slice => {} }, $pathSpec);
    my %orgStepScore = (); # orgId => step => row
    foreach my $row (@$stepScores) {
      $orgStepScore{$row->{orgId}}{$row->{stepId}} = $row;
    }
    my @tr = ();
    my @th = qw{Genome Best-path};
    map s/-/ /, @th;
    push @tr, Tr(th({-valign => "top"}, \@th));
    foreach my $org (@orgsSorted) {
      my $orgId = $org->{orgId};
      my $all = $orgScore{$orgId} || die "No all entry for $orgId and $pathSpec\n";
      die unless $all->{expandedPath};
      my @expandedPath = split / /, $all->{expandedPath};
      my $score = RuleToMinScore($all);
      push @tr, Tr(td({-valign => "top"},
                      [ a({ -href => "gapView.cgi?orgs=$orgsSpec&set=$set&path=$pathSpec&orgId=$orgId",
                            -style => ScoreToStyle($score),
                            -title => "$pathSpec $ruleScoreLabels[$score]" },
                          $orgs{$orgId}{genomeName}),
                        PathToHTML(\@expandedPath, $orgStepScore{$orgId}) ]));
    }
    print table({-cellpadding=>2, -cellspacing=>0, -border=>1}, @tr), "\n";
    print LegendForColorCoding();
    ShowWarnings(undef, $pathSpec); # all orgs, this pathway
  } elsif ($orgId ne "" && $findGene ne "") {
    # mode: Search for a gene
    print p(qq{Searching for "$findGeneShow" in}, $orgs{$orgId}{genomeName}), "\n";
    my $assembly = OrgToAssembly($orgId);
    my $regexp = quotemeta($findGene);
    $regexp =~ s/\\%/.*/g;
    my @hits = (); # each is a list of [locusId, desc]
    if (exists $assembly->{features}) {
      foreach my $row (@{ $assembly->{features} }) {
        next unless $row->{class} eq "with_protein";
        next unless $row->{name} =~ m/$regexp/i
          || $row->{locus_tag} =~ m/^$regexp/i
          || $row->{product_accession} =~ m/^$regexp/i
          || $row->{"non-redundant_refseq"} =~ m/^$regexp/i
          || (exists $assembly->{oldid}{$row->{locus_tag}} && $assembly->{oldid}{$row->{locus_tag}} =~ m/^$regexp/i);
        my $id = $row->{product_accession} || $row->{"non-redundant_refseq"};
        die "Invalid identifier $id in feature file\n"
          unless defined $id && $id ne "";
        my $desc = $row->{name};
        my @moreids = ();
        if ($row->{locus_tag}) {
          push @moreids, $row->{locus_tag};
          push @moreids, $assembly->{oldid}{$row->{locus_tag}}
            if exists $assembly->{oldid}{$row->{locus_tag}};
        }
        push @hits, [ $id, join(" ", @moreids, $desc) ];
      }
    } else {
      die "No faa file for this assembly\n" unless exists $assembly->{faafile};
      open(my $fhA, "<", $assembly->{faafile}) || die "Cannot read $assembly->{faafile}";
      my $state = {};
      while (my ($header, undef) = ReadFastaEntry($fhA, $state)) {
        if ($header =~ m/$regexp/i) {
          my @pieces = split / /, $header;
          my $id = shift @pieces;
          die "Blank header in fasta file" unless defined $id && $id ne "";
          push @hits, [ $id, join(" ", @pieces) ];
        }
      }
      close($fhA) || die "Error reading $assembly->{faafile}";
    }
    if (@hits == 0) {
      print p("No matching proteins were found");
    } else {
      my $cand = $dbhG->selectall_arrayref("SELECT * from Candidate WHERE orgId = ?",
                                           { Slice => {} }, $orgId);
      my %locusNCand = ();
      foreach my $row (@$cand) {
        $locusNCand{$row->{locusId}}++;
        $locusNCand{$row->{locusId2}}++;
      }
      foreach my $row (@hits) {
        my ($id, $desc) = @$row;
        $desc .= small(" (candidate for", $locusNCand{$id}, "steps)")
          if exists $locusNCand{$id};
        print p(a({ -href => "gapView.cgi?orgs=$orgsSpec&set=$set&orgId=$orgId&locusId=$id"}, $id),
                $desc)."\n";
      }
    }
  } elsif ($orgId ne "" && $pathSpec eq "" && $locusSpec eq "") {
    # mode: Overview of pathways for this organism
    my @hr = ("Pathway", span({-title=>"Best path"}, "Steps"));
    my @tr = ();
    push @tr, th({-valign => "top"}, \@hr);
    my $allScores = $dbhG->selectall_arrayref("SELECT * from RuleScore WHERE ruleId = 'all' AND orgId = ?",
                                               { Slice => {} }, $orgId);
    my %pathwayScore = map { $_->{pathwayId} => $_ } @$allScores;
    my $stepScores = $dbhG->selectall_arrayref("SELECT * from StepScore WHERE orgId = ?",
                                               { Slice => {} }, $orgId);
    my %stepScores = (); # pathwayId => stepId => row
    foreach my $row (@$stepScores) {
      $stepScores{ $row->{pathwayId} }{ $row->{stepId} } = $row;
    }
    my @allSorted = ();
    foreach my $pathwayId (@path) {
      my $all = $pathwayScore{$pathwayId}
        || die "Missing result for rule = all, $pathwayId, orgId = $orgId";
      $all->{minScore} = RuleToMinScore($all); # 0, 1, or 2
      $all->{n} = $all->{nHi} + $all->{nMed} + $all->{nLo};
      push @allSorted, $all;
    }
    my $baseURL = "gapView.cgi?orgs=$orgsSpec&set=$set&orgId=$orgId";
    if (param('pathByScore')) {
      # The same sorting as used to choose between paths in gapsummary.pl
      @allSorted = sort { $b->{minScore} <=> $a->{minScore}
                              || $b->{score} <=> $a->{score}
                              || $b->{n} <=> $a->{n} } @allSorted;
      print p(small("Pathways are sorted by completeness.",
                    a({-href => $baseURL}, "Sort by name instead.")));
    } else {
      print p(small("Pathways are sorted by name.",
                    a({-href => "${baseURL}&pathByScore=1"}, "Sort by completeness instead.")));
    }
    foreach my $all (@allSorted) {
      my $pathwayId = $all->{pathwayId};
      my @show = ();
      my @expandedPath = split / /, $all->{expandedPath};
      my $score = $all->{minScore};
      push @tr, Tr({-valign => "top"},
                   td([ a({ -href => "gapView.cgi?orgs=$orgsSpec&set=$set&orgId=$orgId&path=$pathwayId",
                            -style => ScoreToStyle($score),
                            -title => $pathDesc{$pathwayId} . " - " . $ruleScoreLabels[$score] }, $pathwayId),
                        PathToHTML(\@expandedPath, $stepScores{$pathwayId}) ]));
    }
    print table({-cellpadding=>2, -cellspacing=>0, -border=>1}, @tr), "\n";
    print LegendForColorCoding();
    ShowWarnings($orgId, undef); # all pathways, this organism
  } elsif ($orgId ne "" && $pathSpec ne "" && $stepSpec eq "" && $locusSpec eq "") {
    # mode: Show this pathway in this organism:
    # pathway=level gap comment (if any),
    # a summary of the best path, the rules (hierarchically),
    # a table of all of the steps, and the dependency warnings.
    if (exists $knownGaps{$orgId}{$pathSpec}{""}) {
      my $kg = $knownGaps{$orgId}{$pathSpec}{""};
      print p(i("Manual classification of gap:"), $kg->{gapClass},
                br(),
              i("Rationale:"), $kg->{comment});
    }
    my $ruleScores = $dbhG->selectall_hashref("SELECT * from RuleScore WHERE orgId = ? AND pathwayId = ?",
                                             "ruleId",
                                             { Slice => {} }, $orgId, $pathSpec);
    my $stepScores = $dbhG->selectall_hashref("SELECT * from StepScore WHERE orgId = ? AND pathwayId = ?",
                                              "stepId",
                                              { Slice => {} }, $orgId, $pathSpec);
    # Comments are in the steps file only, not in the database.
    my $stepsObj = Steps::ReadSteps("../gaps/$set/$pathSpec.steps");
    my @expandedPath = split / /, $ruleScores->{all}{expandedPath};
    die unless @expandedPath > 0;
    print h3(scalar("Best path")),
      p(i(PathToHTML(\@expandedPath, $stepScores))),
      "\n";

    my @bestcand = (); # the best candidate per step (or the top two, if the same score)
    foreach my $stepId (@expandedPath) {
      my $stepS = $stepScores->{$stepId} || die;
      push @bestcand, $stepS->{locusId} if $stepS->{locusId};
      push @bestcand, $stepS->{locusId2} if $stepS->{locusId2} && ($stepS->{score}+0) == ($stepS->{score2}+0);
    }
    if ($orgs{$orgId}{gdb} eq "FitnessBrowser" && @bestcand > 0) {
      # Link to fitness data for the best candidates
      my @loci = UniqueLoci(SplitLoci(@bestcand));
      my $URL = "http://fit.genomics.lbl.gov/cgi-bin/genesFit.cgi?orgId=$orgs{$orgId}{gid}&"
        . join("&", map { "locusId=$_" } @loci);
      print p("Also see", a({ -href => $URL }, "fitness data"), "for the top candidates");
    }

    print h3("Rules");
    print p("Overview:", LinkifyComment($stepsObj->{topComment})) if $stepsObj->{topComment} ne "";
    print RulesToHTML($stepsObj, $pathSpec, $orgId),
      "\n";

    # Show the steps on the best path, and then the other steps, sorted by their name (case insensitive)
    # Extra columns (not always shown) for curated gaps (same organism and gapClass set) and known gaps
    # (Since pathway-wide curated gaps are shown at the top, they are not shown here)
    my $nWithCand = 0;
    my %stepToCuratedGap;
    my %stepToKnownGap;
    foreach my $stepScore (values %$stepScores) {
      $nWithCand++ if $stepScore->{locusId} ne "";
      my $known = StepScoreToKnownGap($stepScore);
      if (defined $known && $known->{stepId} ne "") {
        if ($known->{orgId} eq $orgId && $known->{gapClass} ne "") {
          $stepToCuratedGap{ $stepScore->{stepId} } = $known;
        } else {
          $stepToKnownGap{ $stepScore->{stepId} } = $known;
        }
      }
    }
    my $hasCurated = keys(%stepToCuratedGap) > 0;
    my $hasKnown = keys(%stepToKnownGap) > 0;
    print h3(scalar(keys %$stepScores), "steps", "(${nWithCand}", "with candidates)"), "\n";
    print p("Or see", a({ -href => "gapView.cgi?set=$set&orgs=$orgsSpec&path=$pathSpec&showdef=1"}, "definitions of steps"));
    my @header = qw{Step Description Best-candidate 2nd-candidate};
    foreach (@header) { s/-/ /; }
    push @header, "Class of gap" if $hasCurated;
    push @header, "Known gap?" if $hasKnown;
    my @tr = ();
    push @tr, Tr(th({ -valign => "top" }, \@header));

    my @bestPathScores = map $stepScores->{$_}, @expandedPath;
    my @otherScores = sort { lc($a->{stepId}) cmp lc($b->{stepId}) } grep ! $_->{onBestPath}, (values %$stepScores);
    my $alternativeShown = 0;
    foreach my $stepScore (@bestPathScores, @otherScores) {
      my $stepId = $stepScore->{stepId};
      if (! $stepScore->{onBestPath} && ! $alternativeShown) {
        push @tr, Tr(td({ -colspan => scalar(@header) }, "Alternative steps:"));
        $alternativeShown = 1;
      }
      my ($show1, $show2) = ShowCandidatesForStep($stepScore);
      my @td = ( StepToShortHTML($stepId, $stepScore),
                 $stepDesc{$pathSpec}{$stepId},
                 $show1, $show2 );
      if ($hasCurated) {
        my $curatedHTML = "&nbsp;";
        $curatedHTML = a({ -title => $stepToCuratedGap{$stepId}{comment}, -style => "color: darkgreen;" },
                         $stepToCuratedGap{$stepId}{gapClass})
          if exists $stepToCuratedGap{$stepId};
        push @td, $curatedHTML;
      }
      if ($hasKnown) {
        my $knownHTML = "&nbsp;";
        if (exists $stepToKnownGap{$stepId}) {
          my $knownGap = $stepToKnownGap{$stepId};
          my $title = "Despite the apparent lack of $stepId, $knownGap->{genomeName}"
            . " performs $pathDesc{$pathSpec}.";
          if (exists $knownGap->{identity}) {
            my $idShow = int(0.5 + $knownGap->{identity});
            $title .= " Across $knownGap->{nMarkers} ribosomal proteins, it is"
              . " ${idShow}% identical to $orgs{$orgId}{genomeName}.";
            $knownHTML = a({ -title => $title }, "known gap (${idShow}% id.)");
          } else {
            $knownHTML = a({ -title => $title }, "known gap");
          }
        }
        push @td, $knownHTML;
      }
      push @tr, Tr(td({ -valign => "top" }, \@td));
    } # end loop over steps to show in table
    print table({-cellpadding=>2, -cellspacing=>0, -border=>1}, @tr), "\n";
    print LegendForColorCoding();
    ShowWarnings($orgId, $pathSpec);
  } elsif ($orgId ne "" && $pathSpec ne "" && $stepSpec ne "" && $locusSpec eq "") {
    # mode: Show step in an organism: if it is a known gap,
    # candidates, fitness link for candidates (if relevant),
    # the step's definition, and known gaps in related organisms (if any)
    my $stepScore = $dbhG->selectrow_hashref("SELECT * FROM StepScore WHERE orgId = ? AND pathwayId = ? AND stepId = ?",
                                             {}, $orgId, $pathSpec, $stepSpec);
    die "No score for $pathSpec $stepSpec and org $orgId" unless $stepScore;
    if ($stepScore->{score} ne "2") {
      my $knownGap = StepScoreToKnownGap($stepScore);
      if ($knownGap && $knownGap->{stepId} ne "") {
        if ($knownGap->{orgId} eq $orgId && $knownGap->{gapClass} ne "") {
          print p(i("Manual classification of gap:"), $knownGap->{gapClass},
                  br(),
                  i("Rationale:"), $knownGap->{comment});
        } elsif ($knownGap->{orgId} eq $orgId) {
          print p(i("Known gap:"),
                  "Despite the apparent lack of $stepSpec,",
                  $knownGap->{genomeName}, "performs", $pathDesc{$pathSpec}.".");
        } else {
          die unless $knownGap->{identity};
          my $idShow = int(0.5 + $knownGap->{identity}) . "%";
          print p(i("Known gap:"), "The related organism",
                  a({ -href => "gapView.cgi?set=$set&gdb=$knownGap->{gdb}&gid=$knownGap->{gid}" },
                    $knownGap->{genomeName}),
                  "performs", $pathDesc{$pathSpec},
                  "and has an apparent gap at $stepSpec.",
                  "Across $knownGap->{nMarkers} ribosomal proteins, the two organisms share",
                  "$idShow amino acid identity.");
        }
      }
    }

    # Candidates table
    my $cand = $dbhG->selectall_arrayref(qq{ SELECT * from Candidate
                                             WHERE orgId = ? AND pathwayId = ? AND stepId = ? },
                                         { Slice => {} }, $orgId, $pathSpec, $stepSpec);
    if (@$cand == 0) {
      print h3("No candidates for $stepSpec: $stepDesc{$pathSpec}{$stepSpec}"), "\n";
    } else {
      my @sorted = sort { ($b->{score} || 0) <=> ($a->{score} || 0)
                            || max($b->{blastBits} || 0, $b->{hmmBits} || 0)
                              <=> max($a->{blastBits} || 0, $a->{hmmBits} || 0 ) } @$cand;
      print h3(scalar(@sorted), "candidates for $stepSpec:", $stepDesc{$pathSpec}{$stepSpec}), "\n";
      my @header = qw{Score Gene Description Similar-to Id. Cov. Bits  Other-hit Other-id. Other-bits};
      $header[-1] = span({-title => "A characterized protein that is similar to the gene but is not associated with step $stepSpec"},
                         $header[-1]);
      foreach (@header) {
        s/-/ /;
      }
      my @tr = Tr(th({-valign => "bottom"}, \@header));

      foreach my $cand (@sorted) {
        # potentially make two rows, one for BLAST and one for HMM
        my $id = a({-href => "gapView.cgi?orgs=$orgsSpec&set=$set&orgId=$orgId&locusId=".$cand->{locusId},
                   -title => $cand->{locusId} },
                   $cand->{sysName} || $cand->{locusId} );
        my $desc = HTML::Entities::encode( $cand->{desc} );
        # HMM hits are based on the 1st ORF only so ignore the split when showing the HMM part
        my $id1 = $id;
        my $desc1 = $desc;
        if ($cand->{locusId2}) { # (this should only happen for BLAST hits)
          $id .= "; " . a({-href => "gapView.cgi?orgs=$orgsSpec&set=$set&orgId=$orgId&locusId=".$cand->{locusId2} },
                          $cand->{sysName2} || $cand->{locusId2} );
          $desc .= "; " . HTML::Entities::encode( $cand->{desc2} );
        }
        my ($linkOther, $otherIdentity, $otherBits) = CandToOtherHTML($cand);
        if ($cand->{blastScore} ne "") {
          push @tr, Tr(td({-valign => "top"},
                          [ ShowScoreShort($cand->{blastScore}),
                            $id, $desc,
                            CuratedToLink($cand->{curatedIds}, $pathSpec, $stepSpec),
                            int(0.5 + $cand->{identity})."%",
                            a({ -href => "gapView.cgi?orgs=$orgsSpec&set=$set&orgId=$orgId&path=$pathSpec&step=$stepSpec&locusId=$cand->{locusId}",
                                -title => "View alignments" },
                              int(0.5 + 100 * $cand->{blastCoverage})."%"),
                            $cand->{blastBits},
                            $linkOther, $otherIdentity, $otherBits ]));
        }

        if ($cand->{hmmScore} ne "") {
          my $hmmURL = HMMToURL($cand->{hmmId});
          push @tr, Tr(td({ -valign => "top" },
                          [ ShowScoreShort($cand->{hmmScore}), $id1, $desc1,
                            $cand->{hmmDesc} . " (" .
                            a({-href => $hmmURL, -title => "curated family (HMM)"}, $cand->{hmmId}) . ")",
                            "",
                            a({ -href => "gapView.cgi?orgs=$orgsSpec&set=$set&orgId=$orgId&path=$pathSpec&step=$stepSpec&locusId=$cand->{locusId}",
                                -title => "View alignments" },
                              int(0.5 + 100 * $cand->{hmmCoverage})."%"),
                            $cand->{hmmBits},
                            $linkOther, $otherIdentity, $otherBits ]));
        }
      }
      print table({-cellpadding=>2, -cellspacing=>0, -border=>1}, @tr), "\n";
      print LegendForColorCoding();
      if ($orgs{$orgId}{gdb} eq "FitnessBrowser") {
        # link to candidates in the fitness browser
        my @loci = UniqueLoci(grep { $_ ne "" } map { $_->{locusId}, $_->{locusId2} } @sorted);
        @loci = splice(@loci, 0, 6) if @loci > 6;
        my $gid = $orgs{$orgId}{gid};
        my $URL = "http://fit.genomics.lbl.gov/cgi-bin/genesFit.cgi?orgId=${gid}&"
          . join("&", map { "locusId=$_" } @loci);
        print p("Also see", a({ -href => $URL }, "fitness data"),"for the candidates");
      }
    } # end if there are candidates

    # Show step definition
    print h3("Definition of step", i($stepSpec)), "\n";

    my $stepsObj = Steps::ReadSteps("../gaps/$set/$pathSpec.steps"); # for the comment
    my $stepObj = $stepsObj->{steps}{$stepSpec} || die;
    my $stepPartData = DataForStepParts($pathSpec, $stepSpec);
    my $stepParts = $dbhS->selectall_arrayref("SELECT * from StepPart WHERE pathwayId = ? AND stepId = ?",
                                              { Slice => {} }, $pathSpec, $stepSpec);
    print start_ul();
    foreach my $row (@$stepParts) {
      print li(FormatStepPart($stepPartData, $row, $orgId));
    }
    print li("Comment:", LinkifyComment($stepObj->{comment}))
      if $stepObj->{comment} ne "";
    print end_ul(), "\n";
    print p("Or cluster all characterized", a({-href => "curatedClusters.cgi?set=$set&path=$pathSpec&step=$stepSpec"},
                         "$stepSpec proteins"));

    # Known gaps in related organisms
    if (exists $markerSim{$orgId}) {
      my @knownGaps = ();
      foreach my $markerSim (@{ $markerSim{$orgId} }) {
        my $orgId2 = $markerSim->{hitOrgId};
        if (exists $knownGaps{$orgId2}{$pathSpec}{$stepSpec}) {
          my %kg = %{ $knownGaps{$orgId2}{$pathSpec}{$stepSpec} };
          next if $kg{gapClass} eq "spurious";
          $kg{identity} = $markerSim->{identity};
          $kg{nMarkers} = $markerSim->{nMarkers};
          push @knownGaps, \%kg;
        }
      }
      if (@knownGaps > 0) {
        print h3("Known gaps in related organisms");
        print start_ul();
        foreach my $kg (@knownGaps) {
          my $idShow = int(0.5 + $kg->{identity}) . "%";
          my @orgParts = split /__/, $kg->{orgId};
          my $gdb = shift @orgParts;
          my $gid = join("__", @orgParts);
          print li(a({-href => "gapView.cgi?set=$set&gdb=$gdb&gid=$gid",
                      -title => "Despite the apparent lack of $stepSpec,"
                      . " $kg->{genomeName}"
                      . " performs $pathDesc{$pathSpec}" },
                     $kg->{genomeName}),
                   "(" . a({-title => "$idShow identity across $kg->{nMarkers} ribosomal proteins"}, $idShow) . ")");
        }
        print end_ul(), "\n";
      } # end if @knownGaps > 0
    } # end if exists $markerSim{$orgId}

  } elsif ($orgId ne "" && $stepSpec eq "" && $locusSpec ne "") {
    # mode: Show a gene, including metadata, candidates, and tools
    # First fetch its header and sequence
    my $tmp = "/tmp/gapView.$locusSpec.$$";
    my $faaCand = "$tmp.genome.faa";
    GetOrgSequence($faafile, $orgId, $locusSpec, $faaCand);
    my %fasta = ReadFastaDesc($faaCand);
    unlink($faaCand);
    my $descs = $fasta{desc};
    my $seqs = $fasta{seq};
    die unless scalar(keys %$seqs) == 1 && scalar(keys %$descs) == 1;
    my ($desc) = values %$descs;
    my ($seq) = values %$seqs;
    print p("Annotation:", $desc);
    print p("Length:", length($seq), "amino acids");
    if ($orgs{$orgId}{gdb} eq "local") {
      print p("Source: uploaded fasta file", small("(hash $orgs{$orgId}{gid})"));
    } else {
      print p("Source:", $orgs{$orgId}{gid}, "in", $orgs{$orgId}{gdb});
    }

    # Candidates
    my $cand = $dbhG->selectall_arrayref("SELECT * from Candidate WHERE orgId = ? AND (locusId = ? OR locusId2 = ?)",
                                         { Slice => {} }, $orgId, $locusSpec, $locusSpec);
    if (@$cand == 0) {
      print p("Not a candidate for any step in $setDesc"), "\n";
    } else {
      print h3("Candidate for", scalar(@$cand), "steps in $setDesc");
      my @header = qw{Pathway Step Score Similar-to Id. Cov. Bits Other-hit Other-id. Other-bits};
      foreach (@header) { s/-/ /; }
      my @tr = Tr(th({-valign => "bottom"}, \@header));
      my @sorted = sort  { ($b->{score} || 0) <=> ($a->{score} || 0)
                            || max($b->{blastBits} || 0, $b->{hmmBits} || 0)
                              <=> max($a->{blastBits} || 0, $a->{hmmBits} || 0 ) } @$cand;
      foreach my $cand (@sorted) {
        my $pathwayId = $cand->{pathwayId};
        my $stepId = $cand->{stepId};
        die "Non-existent step $stepId" unless exists $stepDesc{$pathwayId}{$stepId};
        my ($linkOther, $otherIdentity, $otherBits) = CandToOtherHTML($cand);
        my $pathLink = a({ -href => "gapView.cgi?orgs=$orgsSpec&set=$set&orgId=$orgId&path=$pathwayId" },
                         $pathDesc{$pathwayId} );
        my $stepLink = a({ -href => "gapView.cgi?orgs=$orgsSpec&set=$set&orgId=$orgId&path=$cand->{pathwayId}&step=$cand->{stepId}",
                           -title => $stepDesc{$pathwayId}{$stepId} },
                         $cand->{stepId} );
        if ($cand->{blastScore} ne "") {
          my $asterisk = "";
          $asterisk = a({ -title => "Split hit"}, "*") if $cand->{locusId2};
          push @tr, Tr(td({-valign => "top"},
                          [ $pathLink, $stepLink,
                            ShowScoreShort($cand->{blastScore}),
                            CuratedToLink($cand->{curatedIds}, $pathwayId, $stepId),
                            int(0.5 + $cand->{identity})."%",
                            a({ -href => "gapView.cgi?orgs=$orgsSpec&set=$set&orgId=$orgId&path=$pathwayId&step=$stepId&locusId=$locusSpec",
                                -title => "View alignments" },
                              int(0.5 + 100 * $cand->{blastCoverage})."%") . $asterisk,
                            $cand->{blastBits},
                            $linkOther, $otherIdentity, $otherBits ]));
        }
        if ($cand->{hmmScore} ne "") {
          my $hmmURL = HMMToURL($cand->{hmmId});
          push @tr, Tr(td{-valign => "top"},
                       [ $pathLink, $stepLink,
                         ShowScoreShort($cand->{hmmScore}),
                         a({-href => $hmmURL, }, $cand->{hmmDesc}, "(". $cand->{hmmId} . ")"),
                         "",
                         a({ -href => "gapView.cgi?orgs=$orgsSpec&set=$set&orgId=$orgId"
                             . "&path=$pathwayId&step=$stepId&locusId=$locusSpec",
                             -title => "View alignments" },
                           int(0.5 + 100 * $cand->{hmmCoverage})."%"),
                         $cand->{hmmBits},
                         $linkOther, $otherIdentity, $otherBits ]);
        }
      }
      print table({-cellpadding=>2, -cellspacing=>0, -border=>1}, @tr), "\n";
    }

    # Show sequence analysis tools
    my @seqparts = $seq =~ /.{1,60}/g;
    my $newline = "%0A";
    print
      h3("Sequence Analysis Tools");
    my $URL = GeneURL($orgId,$locusSpec);
    print p("View",
            a({-href => $URL}, "$locusSpec"),
            "at", $orgs{$orgId}{gdb})
      if $URL ne "";
    print
      p(a({-href => "http://papers.genomics.lbl.gov/cgi-bin/litSearch.cgi?query=>${locusSpec}$newline$seq"},
          "PaperBLAST"),
        "(search for papers about homologs of this protein)"),
      p(a({-href => "http://www.ncbi.nlm.nih.gov/Structure/cdd/wrpsb.cgi?seqinput=>${locusSpec}$newline$seq"},
          "Search CDD"),
        "(the Conserved Domains Database, which includes COG and superfam)"),
      p(a({-href => "http://pfam.xfam.org/search/sequence?seqOpts=&ga=0&evalue=1.0&seq=$seq"},
        "Search PFam"),
        "(including for weak hits, up to E = 1)"),
      p("Predict protein localization: ",
        a({-href => "http://www.psort.org/psortb/results.pl?"
           . join("&",
                  "organism=bacteria",
                  "gram=negative",
                  "format=html",
                  "sendresults=display",
                  "email=",
                  "seqs=>${locusSpec}$newline$seq")},
          "PSORTb"),
        "(Gram negative bacteria)"),
      p("Predict transmembrane helices:",
        a({-href => "http://www.cbs.dtu.dk/cgi-bin/webface2.fcgi?"
           . join("&",
                  "configfile=/usr/opt/www/pub/CBS/services/TMHMM-2.0/TMHMM2.cf",
                  "outform=-noshort",
                  "SEQ=>${locusSpec}$newline$seq")},
          "TMHMM")),
      p("Check the SEED with",
      a({-href => "http://pubseed.theseed.org/FIG/seedviewer.cgi?page=FigFamViewer&fasta_seq=>${locusSpec}$newline$seq"},
        "FIGfam search")),
      p(a({-title => "Fitness BLAST compares a sequence to bacterial proteins that have mutant phenotypes"},
          "Fitness BLAST:"),
        span({-id => "fitblast_short"}, small("loading..."))),
      qq{<SCRIPT src="https://fit.genomics.lbl.gov/d3js/d3.min.js"></SCRIPT>
         <SCRIPT src="https://fit.genomics.lbl.gov/images/fitblast.js"></SCRIPT>
         <SCRIPT>
         var server_root = "https://fit.genomics.lbl.gov/";
         var seq = "$seq";
         fitblast_load_short("fitblast_short", server_root, seq);
         </SCRIPT>},
      h3("Sequence"),
      join("\n", "<pre>", @seqparts, "</pre>"), "\n";
  } elsif ($locusSpec ne "" && $stepSpec ne "") {
    #XX die unless $pathSpec ne "" & $orgId ne "";
    # mode: Show alignments for a gene
    push @links, a({-href => "gapView.cgi?orgs=$orgsSpec&set=$set&orgId=$orgId&path=$pathSpec&step=$stepSpec"},
                    "All candidates for step $stepSpec in", $orgs{$orgId}{genomeName});
    my $cand = $dbhG->selectall_arrayref(qq{ SELECT * from Candidate
                                             WHERE pathwayId = ? AND stepId = ? AND orgId = ?
                                             AND (locusId = ? OR locusId2 = ?) },
                                         { Slice => {} }, $pathSpec, $stepSpec, $orgId, $locusSpec, $locusSpec);
    foreach my $cand (@$cand) {
      if ($cand->{blastBits} ne "" && $cand->{blastBits} > 0) {
        my $curatedIds = $cand->{curatedIds};
        my @loci = (); # locusId, sysName, desc
        push @loci, [ $cand->{locusId}, $cand->{sysName}, $cand->{desc} ];
        push @loci, [ $cand->{locusId2}, $cand->{sysName2}, $cand->{desc2} ] if $cand->{locusId2} ne "";
        foreach my $row (@loci) {
          my ($locusId, $sysName, $desc) = @$row;
          # Should move the descShowCurated/idShowHit code above to a subroutine for showing what it hits
          print p(b("Align", CuratedToLink($curatedIds, $stepSpec, $stepSpec),
                    br(),
                    "to candidate", $locusId, $sysName, $desc));
          # Figure out the sequence
          my $curatedSeq = CuratedToSeq($curatedIds);;
          my $faaCurated = "$tmp.curated.faa";
          open(my $fhC, ">", $faaCurated) || die "Cannot write to $faaCurated";
          my @hitIds = split /,/, $cand->{curatedIds};
          print $fhC ">$hitIds[0]\n$curatedSeq\n";
          close($fhC) || die "Error writing to $faaCurated";
          my $faaCand = "$tmp.genome.faa";
          GetOrgSequence($faafile, $orgId, $locusId, $faaCand);
          my $bl2seq = "../bin/blast/bl2seq";
          die "No such executable: $bl2seq\n" unless -x $bl2seq;
          print "<pre>";
          system($bl2seq, "-p", "blastp", "-i", $faaCurated, "-j", $faaCand, "-e", 0.01, "-F", "m S") == 0
            || die "bl2seq failed: $!";
          unlink($faaCurated);
          unlink($faaCand);
          print "</pre>\n";
        }
        # Arguably, should show alignments to "other" as well
      }
      if ($cand->{hmmBits} ne "" && $cand->{hmmBits} > 0) {
        print p(b("Align candidate",
                  HTML::Entities::encode($cand->{locusId}),
                  HTML::Entities::encode($cand->{sysName}),
                  HTML::Entities::encode("($cand->{desc})"),
                  br(), "to HMM $cand->{hmmId} ($cand->{hmmDesc})"));
        my $query = $dbhS->selectrow_hashref(qq{ SELECT * from StepQuery
                                                 WHERE pathwayId = ? AND stepId = ? AND hmmId = ? },
                                             { Slice => {} }, $pathSpec, $stepSpec, $cand->{hmmId});
        die "Unknown hmm $cand->{hmmId} for $pathSpec $stepSpec" unless defined $query;
        my $hmmFile = "$stepsDir/$query->{hmmFileName}";
        die "No file for $cand->{hmmId}: $hmmFile is missing\n" unless -e $hmmFile;
        my $hmmsearch = "../bin/hmmsearch";
        die "No such executable: $hmmsearch\n" unless -x $hmmsearch;
        my $faaCand = "$tmp.genome.faa";
        FetchSeqs("../bin/blast", $faafile, [$orgId.":".$cand->{locusId}], $faaCand);
        print "<pre>";
        system($hmmsearch, $hmmFile, $faaCand) == 0
          || die "hmmsearch failed: $!";
        print "</pre>\n";
        unlink($faaCand);
      }
    } # end loop over candidates
  } else {
    die "Unknown mode\n";
  }

  my $format = quotemeta("+%b %d %Y");
  my $dateQuery = `date -r $stepsDir/steps.db $format`;
  my $dateAnalysis = `date -r $sumPre.db $format`;
  chomp $dateAnalysis;
  chomp $dateQuery;
  print p("This GapMind analysis is from $dateAnalysis. The underlying query database was built on $dateQuery.")
    if $findGene eq "" && ! param("showdef");

  if ($orgId ne "") {
    push @links, a({ -href => "gapView.cgi?orgs=$orgsSpec&set=$set&orgId=$orgId" },
                   "$setDesc in", $orgs{$orgId}{genomeName})
      if $pathSpec ne "" || $locusSpec ne "" || $findGene ne "" || param('gaps');
    push @links, a({ -href => "gapView.cgi?orgs=$orgsSpec&set=$set&orgId=$orgId&path=$pathSpec" },
                   "$pathDesc{$pathSpec} in $orgs{$orgId}{genomeName}")
      if $pathSpec ne "" && $stepSpec ne "";
    my @form1 = (start_form(-method => 'get', -action => "genomeSearch.cgi"),
                 hidden(-name => 'gid', -value => $orgs{$orgId}{gid}, -override => 1),
                 hidden(-name => 'gdb', -value => $orgs{$orgId}{gdb}, -override => 1),
                 a({ -title => "Find characterized proteins whose descriptions match"
                     . " and have homologs in this genome",
                     -href => "genomeSearch.cgi?gid=$orgs{$orgId}{gid}&gdb=$orgs{$orgId}{gdb}" },
                   "Curated BLAST:"),
                 textfield(-name => 'query', -value => '', -size => 30, -maxlength => 200),
                 submit("Go"),
                 end_form);
    push @links, join("\n", @form1);
    my @form2 = (start_form(-method => 'get', -action => 'gapView.cgi'),
                 hidden(-name => 'orgs', -value => $orgsSpec, -override => 1),
                 hidden(-name => 'set'),
                 hidden(-name => 'orgId'),
                 a({-title => "Search through the gene descriptions."
                    . " You can use % as a wild-card character"}, "Search annotations: "),
                 textfield(-name => 'findgene', -value => '', -size => 30, -maxlength => 200, -override => 1),
                 submit("Go"),
                 end_form);
    push @links, join("", @form2);
    push @links, join(" ",
                      a({ -href => OrgIdToURL($orgId) }, $orgs{$orgId}{genomeName}),
                      small('(' . $orgs{$orgId}{gid} . ')'),
                      "at", $orgs{$orgId}{gdb})
      unless $orgs{$orgId}{gdb} eq "local";
  }

  push @links, a({-href => "gapView.cgi?orgs=$orgsSpec&set=$set&path=$pathSpec&showdef=1" },
                 "Definition of $pathDesc{$pathSpec}")
    if $pathSpec ne "" && !param("showdef");
  push @links, a({-href => "gapView.cgi?orgs=$orgsSpec&set=$set&path=$pathSpec"},
                 "$pathDesc{$pathSpec} across $nOrgs  genomes")
    if $pathSpec ne "" && (param("showdef") || $orgId ne "") && @orgs > 1;
  if (!param('gaps') && $pathSpec eq "") {
    if ($orgId eq "") {
      push @links, a({ -href => "gapView.cgi?orgs=$orgsSpec&set=$set&gaps=1" },
                     "Potential gaps across $nOrgs genomes and all pathways");
    } else {
      push @links, a({ -href => "gapView.cgi?orgs=$orgsSpec&set=$set&orgId=$orgId&gaps=1" },
                     "Potential gaps in", $orgs{$orgId}{genomeName});
    }
  }
  push @links, a({ -href => "gapView.cgi?orgs=$orgsSpec&set=$set"}, "All $nOrgs genomes and all pathways")
    unless ($orgId eq "" && $pathSpec eq "" && !param('gaps')) || @orgs == 1;
  my $otherSet = $set eq "aa" ? "carbon" : "aa";
  my $otherSetDesc = $otherSet eq "aa" ? "amino acid biosynthesis" : "catabolism of small carbon sources";
  my $otherSetURL = "gapView.cgi?orgs=$orgsSpec&set=$otherSet";
  $otherSetURL .= "&orgId=$orgId" if $orgId ne "";
  push @links, a({ -href => $otherSetURL  }, "GapMind for $otherSetDesc");

  print h3("Links"), start_ul(), li(\@links), end_ul
    if @links > 0;
  print h3("Downloads"),
    start_ul(),
    li(a({-href => "$sumPre.cand"}, "Candidates"), "(tab-delimited)"),
    li(a({-href => "$sumPre.steps"}, "Steps"), "(tab-delimited)"),
    li(a({-href => "$sumPre.rules"}, "Rules"), "(tab-delimited)"),
    li(a({-href => "$orgPre.faa"}, "Protein sequences"), "(fasta format)"),
    li(a({-href => "$orgPre.org"}, "Organisms"), "(tab-delimited)"),
    li("SQLite3 databases"),
    start_ul(),
    li(a({-href => "$stepsDir/curated.db" }, "Curated proteins")),
    li(a({-href => "$stepsDir/steps.db" }, "Rules, steps, and queries")),
    li(a({-href => "$sumPre.db" }, "Analysis results")),
    end_ul,
    end_ul;
  Finish();
}

sub GetOrgSequence($$$$) {
  my ($faaIn, $orgId, $locusSpec, $faaOut) = @_;
  my $db = "$faaIn.db";
  my $id = $orgId.":".$locusSpec;
  if (-e $db) {
    # use DB_File -- faster and more reliable than fastacmd
    my %seqs;
    tie %seqs, "DB_File", $db, O_RDONLY, 0666, $DB_HASH
      || die "Cannot open file $db -- $!";
    my $seq = $seqs{$id};
    untie %seqs;
    die "No sequence for $id in $db" unless defined $seq;
    open(my $fh, ">", $faaOut) || die "Cannot write to $faaOut";
    print $fh ">$id\n$seq\n";
    close($fh) || die "Error writing to $faaOut";
  } else {
    # older builds -- this uses fastacmd
    FetchSeqs("../bin/blast", $faaIn, [$id], $faaOut);
  }
}


sub FaaToDb($$) {
  my ($faaIn, $db) = @_;
  my %seqs;
  tie %seqs, "DB_File", $db, O_CREAT|O_RDWR, 0666, $DB_HASH
    || die "Cannot write to file $db -- $!";
  open (my $fh, "<", $faaIn) || die "Cannot read $faaIn";
  my $state = {};
  while (my ($header, $seq) = ReadFastaEntry($fh, $state)) {
    $header =~ s/ .*//;
    $seqs{$header} = $seq;
  }
  close($fh) || die "Error reading $faaIn";
  untie %seqs;
}

sub ProcessUpload($) {
  my ($fhUp) = @_;
  my $state = {};
  my %seq = (); # header => sequence
  my %ids = (); # ids seen already
  while (my ($header, $seq) = ReadFastaEntry($fhUp, $state, 1)) {
    return ('error' => "Invalid characters in sequence for $header -- only A-Z or * are allowed")
      unless $seq =~ m/^[A-Z*]+$/;
    my @pieces = split / /, $header;
    my $id = shift @pieces;
    return ('error' => "No identifier at beginning of header $header")
      unless defined $id && $id ne "";
    return ('error' => "Duplicate sequence for identifier $id")
      if exists $ids{$id};
    return ('error' => qq{Sorry, '$id' is an invalid protein identifier. The only characters allowed in the first word of each header line are any of: $charsInId})
      unless $id =~ m/^[$charsInId]+$/;
    $ids{$id} = 1;
    $seq{$header} = $seq;
  }
  return ('error' => "Sorry, uploaded file is not a valid fasta file: $state->{error}")
    if exists $state->{error};
  return ('error' => "No input sequences in upload")
    if scalar(keys %seq) == 0;
  return ('error' => "Sorry, the input has too many sequences (the limit is $maxNSeqsK,000)")
    if scalar(keys %seq) > $maxNSeqs;

  my $totLen = 0;
  my $nNucChar = 0;
  foreach my $seq (values %seq) {
    $totLen += length($seq);
    $nNucChar += ($seq =~ tr/ACGTUN//);
  }
  return ('error' => "The uploaded sequences seem to be nucleotide sequences rather than protein sequences")
    if $nNucChar/ $totLen >= 0.9;

  return('error' => "Uploaded protein sequences are limited to 10MB")
    if $totLen > 10e6;
  my $assembly = AASeqToAssembly(\%seq, $dataDir) || die;
  return ('gdb' => $assembly->{gdb}, 'gid' => $assembly->{gid});
}

sub ScoreToStyle($) {
  my ($score) = @_;
  $score = 0 if $score eq "";
  my $color = $score > 1 ? "#007000" : ($score < 1 ? "#CC4444" : "#000000");
  return "color: $color; font-weight: bold;" if $score > 1;
  return "color: $color;";
}

sub ScoreToLabel($) {
  my ($score) = @_;
  $score = 0 if $score eq "";
  return $score > 1 ? "high confidence" : ($score < 1 ? "low confidence" : "medium confidence");
}

sub ShowScoreShort($) {
  my ($score) = @_;
  $score = 0 if $score eq "";
  return span({ -style => ScoreToStyle($score), -title => ScoreToLabel($score) },
              $score > 1 ? "hi" : ($score < 1 ? "lo" : "med"));
}

sub RuleToMinScore($) {
  my ($ruleScore) = @_;
  die "Invalid input to RuleToMinScore"
    unless defined $ruleScore && ref $ruleScore eq "HASH" && exists $ruleScore->{nLo};
  return 0 if $ruleScore->{nLo} > 0;
  return 1 if $ruleScore->{nMed} > 0;
  return 2;
}

sub HMMToURL($) {
  my ($hmmId) = @_;
  if ($hmmId =~ m/^TIGR/) {
    return "https://www.ncbi.nlm.nih.gov/Structure/cdd/".$hmmId;
  } elsif ($hmmId =~ m/^PF/) {
    my $hmmIdShort = $hmmId; $hmmIdShort =~ s/[.]\d+$//;
    return "http://pfam.xfam.org/family/".$hmmIdShort;
  }
  return "";
}

sub GeneURL($$) {
  my ($orgId,$locusId) = @_;
  die unless exists $orgs{$orgId};
  my $gdb = $orgs{$orgId}{gdb};
  my $gid = $orgs{$orgId}{gid};
  if ($gdb eq "FitnessBrowser") {
    return "http://fit.genomics.lbl.gov/cgi-bin/singleFit.cgi?orgId=${gid}&locusId=${locusId}";
  } elsif ($gdb eq "MicrobesOnline") {
    return "http://www.microbesonline.org/cgi-bin/fetchLocus.cgi?locus=$locusId";
  } elsif ($gdb eq "UniProt") {
    return "http://www.uniprot.org/uniprot/$locusId";
  } elsif ($gdb eq "NCBI") {
    my $assembly = OrgToAssembly($orgId);
    if (exists $assembly->{prot}{$locusId}) {
      my $g = $assembly->{prot}{$locusId};
      return "https://www.ncbi.nlm.nih.gov/gene/?term=" . $g->{GeneID}
        if $g->{GeneID};
      if ($g->{genomic_accession} && $g->{start} && $g->{end}) {
        my $center = int(($g->{start} + $g->{end})/2);
        my ($left,$right) = ($center-5000,$center+5000);
        # The NCBI sequence viewer is smart enough to clip to valid regions
        return "https://www.ncbi.nlm.nih.gov/nuccore/$g->{genomic_accession}/scaffold?report=graph&v=$left:$right";
      }
      # No longer build PaperBLAST links
      #return "http://papers.genomics.lbl.gov/cgi-bin/litSearch.cgi?query=".$g->{"non-redundant_refseq"}
      #  if (exists $g->{"non-redundant_refseq"} && $g->{"non-redundant_refseq"}) {
    }
    return "";
  } elsif ($gdb eq "IMG") {
    return "https://img.jgi.doe.gov/cgi-bin/w/main.cgi?section=GeneDetail&page=geneDetail&gene_oid=$locusId";
  } elsif ($gdb eq "local") {
    return "";
  }
  die "Unknown genome database $gdb\n";
}

my $assembly = undef; # cached
sub OrgToAssembly($) {
  my ($orgId) = @_;
  die "Invalid orgId $orgId" unless exists $orgs{$orgId};
  my $gdb = $orgs{$orgId}{gdb};
  my $gid = $orgs{$orgId}{gid};

  if (!defined $assembly || $assembly->{gdb} ne $gdb || $assembly->{gid} ne $gid) {
    $assembly = CacheAssembly($gdb, $gid, $dataDir)
      || die "Cannot fetch assembly $gid from database $gdb\n";
  }
  return $assembly;
}

sub OrgIdToURL($) {
  my ($orgId) = @_;
  die "Invalid orgId" unless exists $orgs{$orgId};
  my $gdb = $orgs{$orgId}{gdb};
  my $gid = $orgs{$orgId}{gid};
  die unless $gdb && $gid;
  # Avoid running any queries or fetching the genome (via CacheAssembly)
  # if all we need is a link to the genome page
  return "https://www.ncbi.nlm.nih.gov/assembly/$gid"
    if $gdb eq "NCBI";
  return "http://fit.genomics.lbl.gov/cgi-bin/org.cgi?orgId=$gid"
    if $gdb eq "FitnessBrowser";
  return "http://www.microbesonline.org/cgi-bin/genomeInfo.cgi?tId=$gid"
    if $gdb eq "MicrobesOnline";
  # the other databases are trickier, so just fetch everything
  my $assembly = OrgToAssembly($orgId);
  return $assembly->{URL};
}

sub Finish() {
  my $set = param("set") || "aa";
  my $otherSet = $set eq "aa" ? "carbon" : "aa";
  my $otherDesc = $otherSet eq "aa" ? "GapMind for amino acid biosynthesis" : "GapMind for catabolism of small carbon sources";
  print
    h3("Related tools"),
    start_ul();
  print li(a{ -href => "gapView.cgi?set=$otherSet" }, $otherDesc) if $orgsSpec eq "";
  print
    li(a({ -href => "litSearch.cgi" }, "PaperBLAST: Find papers about a protein or its homologs")),
    li(a({ -href => "genomeSearch.cgi", -title => "Search a genome for proteins that are related to a query term" },
         "Curated BLAST for Genomes")),
    end_ul();

  my $email = 'funwithwords26@gmail.com';
  print <<END
<h3>About GapMind</h3>
<P>Each pathway is defined by a set of rules based on individual steps or genes. Candidates for each step are identified by using ublast against a database of manually-curated proteins (most of which are experimentally characterized) or by using HMMer. Ublast hits may be split across two different proteins.

<P>A candidate for a step is "high confidence" if either:
<UL>
<LI>ublast finds a hit to a characterized protein at above 40% identity and 80% coverage, and bits >= other bits+10.
<UL><LI>(Hits to curated proteins without experimental data as to their function are never considered high confidence.)</UL>
<LI>HMMer finds a hit with 80% coverage of the model, and either other identity < 40 or other coverage < 0.75.
</UL>
where "other" refers to the best ublast hit to a sequence that is not annotated as performing this step (and is not "ignored").

<P>Otherwise, a candidate is "medium confidence" if either:
<UL>
<LI>ublast finds a hit at above 40% identity and 70% coverage (ignoring otherBits).
<LI>ublast finds a hit at above 30% identity and 80% coverage, and bits >= other bits.
<LI>HMMer finds a hit (regardless of coverage or other bits).
</UL>

<P>Other blast hits with at least 50% coverage are "low confidence."
<P>Steps with no high- or medium-confidence candidates may be considered "gaps."
For the typical bacterium that can make all 20 amino acids, there are 1-2 gaps in amino acid biosynthesis pathways.
Gaps may be due to:
<UL>
<LI>our ignorance of proteins' functions,
<LI>omissions in the gene models,
<LI>frame-shift errors in the genome sequence, or
<LI>the organism lacks the pathway.
</UL>

<P>GapMind relies on the predicted proteins in the genome and does not search the six-frame translation. In most cases, you can search the six-frame translation by clicking on links to Curated BLAST for each step definition (in the per-step page).
END
    ;

  my @info = ("For more information, see the",
              a({ -href => "https://msystems.asm.org/content/5/3/e00291-20",
                  -title => "GapMind: Automated annotation of amino acid biosynthesis" },
                "paper from 2019"),
              "on GapMind for amino acid biosynthesis, or view the",
              a({ -href => "https://github.com/morgannprice/PaperBLAST" }, "source code"));
  my $changesFile = "../gaps/$set/changes";
  if (defined $setDesc && -e $changesFile) {
    $info[-1] .= ",";
    push @info, ("or see", a({ -href => $changesFile }, "changes"), "to", i($setDesc), "since the publication.");
  } else {
    $info[-1] .= ".";
  }
  print p(@info);
  print p("If you notice any errors or omissions in the step descriptions, or any questionable results, please",
          a( { -href => "mailto:$email" }, "let us know"));
  print p({ -align => 'center' }, "by",
          a( { -href => "http://morgannprice.org/" }, "Morgan Price").",",
          a( { -href => "http://genomics.lbl.gov/" }, "Arkin group").",",
          "Lawrence Berkeley National Laboratory");
  print end_html;
  exit(0);

}

sub StepScoreToKnownGap($) {
  my ($row) = @_;
  my $orgId = $row->{orgId};
  my $pathwayId = $row->{pathwayId};
  my $stepId = $row->{stepId};
  my $score = $row->{score};
  die unless defined $orgId && defined $pathwayId && defined $stepId && defined $score;
  return undef if $score eq "2";
  if (exists $knownGaps{$orgId}{$pathwayId}{$stepId}) {
    return $knownGaps{$orgId}{$pathwayId}{$stepId}
  }
  #else
  # For curated organisms, use cases with step=""
  if (exists $knownGaps{$orgId}{$pathwayId}{""}) {
    return $knownGaps{$orgId}{$pathwayId}{""};
  }
  #else
  # (Ignores similarity to other org if score is 1)
  if (exists $markerSim{$orgId} && $score ne "1") {
    foreach my $markerSim (@{ $markerSim{$orgId} }) {
      my $orgId2 = $markerSim->{hitOrgId};
      if (exists $knownGaps{$orgId2}{$pathwayId}{$stepId}) {
        my %out = %{ $knownGaps{$orgId2}{$pathwayId}{$stepId} };
        next if $out{gapClass} eq "spurious"; # Do not map spurious gaps
        $out{identity} = $markerSim->{identity};
        $out{nMarkers} = $markerSim->{nMarkers};
        return \%out;
      }
    }
  }
  return undef;
}

# Given a row from the step score table (or undef),
# returns formatted HTML for locusId/sysName/score
# (including a link to the page for the gene)
# for the top two candidates
sub ShowCandidatesForStep($) {
  my ($stepScoreRow) = @_;
  return ("","") unless defined $stepScoreRow && $stepScoreRow->{orgId};
  my $orgId = $stepScoreRow->{orgId};
  my @work = ();
  push @work, [ $stepScoreRow->{locusId}, $stepScoreRow->{sysName},  $stepScoreRow->{score} ]
    if $stepScoreRow->{locusId} ne "";
  push @work, [ $stepScoreRow->{locusId2}, $stepScoreRow->{sysName2}, $stepScoreRow->{score2} ]
    if $stepScoreRow->{locusId2} ne "";

  my @show = ();
  foreach my $work (@work) {
    my ($locusId,$sysName,$score) = @$work;
    # Create two links if this is a split hit
    my @sysNameParts = split /,/, $sysName;
    my @locusParts = split /,/, $locusId;
    die unless @locusParts >= 1 && @locusParts <= 2;
    my $locus1 = $locusParts[0];
    my $locus2 = "";
    $locus2 = $locusParts[1] if @locusParts  > 1;
    my $cand = $dbhG->selectrow_hashref("SELECT * FROM Candidate WHERE orgId=? AND pathwayId=? AND stepId=? AND locusId=? AND locusId2=?",
                                    {},
                                    $orgId, $stepScoreRow->{pathwayId}, $stepScoreRow->{stepId},
                                    $locus1, $locus2);
    die "No cand for " . join(" ", %$stepScoreRow) unless $cand;
    my @descParts = ($cand->{desc}, $cand->{desc2});
    my @parts = ();
    foreach my $i (0..(scalar(@locusParts)-1)) {
      my $locusPart = $locusParts[$i];
      my $desc = $cand->{$i == 0 ? "desc" : "desc2"};
      my $title = ScoreToLabel($score);
      $title .= ", annotated as " . HTML::Entities::encode($desc)
        if $desc =~ m/[a-zA-Z0-9]/; # ignore empty descriptions
    push @parts, a({ -style => ScoreToStyle($score),
                     -title => $title,
                     -href => "gapView.cgi?orgs=$orgsSpec&set=$set&orgId=$orgId&locusId=$locusPart" },
                   $sysNameParts[$i] || $locusPart );
    }
    if (@locusParts > 1) {
      push @show, join(" " . a({-title => "split protein"}, "with") . " ", @parts);
    } else {
      push @show, $parts[0];
    }
  }
  while (@show < 2) {
    push @show, "";
  }
  return @show;
}

sub LegendForColorCoding() {
  my @titles = ("Low confidence candidates are highly diverged, have low coverage of the characterized homolog, or are similar to proteins that have other functions.",
                "Medium confidence candidates are less than 40% identical to a characterized protein; or the alignment (to either a characterized protein or an HMM) had under 80% coverage; or the candidate was found by similarity to a uncharacterized (but well-curated) protein.",
                "High confidence candidates match an HMM or are over 40% similar to a characterized protein; and the alignment covers 80% of the characterized protein or the HMM; and the candidate is less similar to characterized proteins that have other functions.");
  my @showScores = map span({ -style => ScoreToStyle($_), -title => $titles[$_] },
                            ScoreToLabel($_)), (2,1,0);
  my @parts = ("Confidence:", @showScores, br(),
               "?", "&ndash;", "known gap:",
               "despite the lack of a good candidate for this step,",
               "this organism (or a related organism) performs the pathway"
              );
  push @parts, (br(), span({-style => $transporterStyle}, "transporter"),
                "&ndash;", "transporters and PTS systems are shaded because",
                "predicting their specificity is particularly challenging.")
    if $set eq "carbon" && ! param('gaps');
  return p(@parts)."\n";
}

sub ShowWarnings($$) {
  my ($orgIdFilter, $pathSpec) = @_;
  my @clauses;
  push @clauses, qq{orgId = "$orgIdFilter"} if defined $orgIdFilter;
  push @clauses, qq{pathwayId = "$pathSpec"} if $pathSpec;
  my $where = "";
  $where = "WHERE " . join(" AND ", @clauses) if @clauses > 0;
  my $reqNotMet = $dbhG->selectall_arrayref("SELECT * FROM RequirementNotMet $where",
                                            { Slice => {} });
  return if @$reqNotMet == 0;
  print h3("Dependencies"), start_ul;
  foreach my $warn (@$reqNotMet) {
    my $pathShow = a({-href => "gapView.cgi?orgs=$orgsSpec&set=$set&path=" . $warn->{pathwayId}
                      . "&orgId=" . $warn->{orgId} }, $warn->{pathwayId});
    my $gn = $orgs{ $warn->{orgId} }{genomeName};
    my $out = scalar(keys %orgs) > 1 && $orgIdFilter eq "" ? "In $gn, $pathShow" : $pathShow;
    $out .= " (using rule $warn->{ruleId})" if $warn->{ruleId} ne "all";
    my $reqShow = $warn->{requiredPathwayId};
    my $reqPart = $warn->{requiredRuleId} || $warn->{requiredStepId} || "";
    $reqShow .= ":" . $reqPart unless $reqPart eq "all";
    my $reqURL = "gapView.cgi?orgs=$orgsSpec&set=$set&path=" . $warn->{requiredPathwayId}
      . "&orgId=" . $warn->{orgId};
    if ($warn->{isNot}) {
      $out .= " is not allowed with " . a({-href => $reqURL}, $reqShow);
    } else {
      $out .= " also requires " . a({ -href => $reqURL }, $reqShow)
        . ", which is not high-confidence";
    }
    $out .= ".";
    $out .= start_ul . li($warn->{comment}) . end_ul if $warn->{comment} ne "";
    print li($out);
  }
  print end_ul, "\n";
}

sub PathToHTML($$) {
  my ($stepList, $stepScores) = @_;
  my @out = ();
  die unless $stepList && @$stepList > 0;
  foreach my $stepId (@$stepList) {
    my $stepScore = $stepScores->{$stepId} || die "No scoring for step $stepId";
    push @out, StepToShortHTML($stepId, $stepScore);
  }
  return join(", ", @out);
}

sub StepToShortHTML($$) {
  my ($stepId, $stepScore) = @_;
  die unless $stepId eq $stepScore->{stepId};
  my $orgId = $stepScore->{orgId} || die;
  my $pathwayId = $stepScore->{pathwayId} || die;
  my $score = $stepScore->{score} || 0;
  my $id = $stepScore->{sysName} || $stepScore->{locusId} || "";
  # For split ORFs, if sysNames are empty then the joint one will be ","
  # which is not useful to show
  $id = $stepScore->{locusId} if $id eq ",";
  my $title = $stepDesc{$pathwayId}{$stepId} || die;
  $title .= " $id" if $id ne "";
  my $knownGap = StepScoreToKnownGap($stepScore);
  if (defined $knownGap && $knownGap->{gapClass} ne "") {
    $title .= " (" . $knownGap->{gapClass} . " gap)";
  } elsif ($knownGap && $knownGap->{orgId} eq $orgId) {
    $title .= " (a known gap)";
  } elsif ($knownGap) {
    my $idShow = int(0.5 + $knownGap->{identity})."%";
    $title .= " (this is a known gap in $knownGap->{genomeName}, which is $idShow identical across $knownGap->{nMarkers} ribosomal proteins)";
  } else {
    $title .= " (" . ScoreToLabel($score) . ")";
  }
  my $showMaybe = defined $knownGap? a({-title => $title }, small("?")) : "";
  my $showSplit = "";
  if ($stepScore->{locusId} =~ m/,/) {
    my $ids = $id;
    $ids =~ s/,/ and /;
    # Do not superscript the asterisk, the result is too subtle
    $showSplit = a({ -title => "In this organism, $stepId may be split across two proteins: $ids",
                     -style => "font-weight: bold;" }, "*")
  }
  my $style = ScoreToStyle($score);
  my ($isTransport) = $dbhS->selectrow_array("SELECT isTransport FROM Step WHERE pathwayId = ? AND stepId = ?",
                                             {}, $pathwayId, $stepId);
  $style .= $transporterStyle if $isTransport;
  return a({ -href => "gapView.cgi?orgs=$orgsSpec&set=$set&orgId=$orgId&path=$pathwayId&step=$stepId",
             -style => $style,
             -title => $title },
           $stepId) . $showMaybe . $showSplit;
}

sub UniqueLoci(@) {
  my %seen = ();
  my @out = ();
  foreach my $in (@_) {
    push @out, $in if !exists $seen{$in};
    $seen{$in} = 1;
  }
  return @out;
}

sub SplitLoci(@) {
  my @out = ();
  foreach my $in (@_) {
    push @out, split /,/, $in;
  }
  return @out;
}

sub CuratedToLink($$$) {
  my ($curatedIds, $pathwayId, $stepId) = @_;
  die "Undefined curatedIds" unless defined $curatedIds;
  my $stepDefURL = "gapView.cgi?orgs=$orgsSpec&set=$set&showdef=1&path=$pathwayId#$stepId";
  my $curatedDesc;
  if ($curatedIds =~ m/^curated2:(.*)$/) {
    my $protId = $1;
    ($curatedDesc) = $dbhC->selectrow_array("SELECT desc FROM Curated2 WHERE protId = ?",
                                            {}, $protId);
    die $protId unless defined $curatedDesc;
  } elsif ($curatedIds =~ m/^uniprot:(.*)$/) {
    my $uniprotId = $1;
    ($curatedDesc) = $dbhS->selectrow_array("SELECT desc FROM StepQuery WHERE uniprotId = ?",
                                            {}, $uniprotId);
    die $uniprotId unless defined $curatedDesc;
  } else {
    ($curatedDesc) = $dbhC->selectrow_array("SELECT descs FROM CuratedInfo WHERE curatedIds = ?",
                                            {}, $curatedIds);
    die $curatedIds unless defined $curatedDesc;
  }
  $curatedDesc =~ s/;;.*//;
  my ($first) = split /,/, $curatedIds;
  my $charLabel = "characterized";
  my $idShowHit = $first;
  $idShowHit =~ s/^.*://;
  my $charTitle = "$idShowHit has been studied experimentally.";
  if ($first =~ m/^curated2:/) {
    $first =~ s/^curated2://;
    $charLabel = "uncharacterized";
    $charTitle = qq{$idShowHit has not been studied experimentally.
                    It is included in GapMind's database because its annotation was manually curated.};
    $charTitle =~ s/\s+/ /g;
  } elsif ($first =~ m/^uniprot:/) {
    $first =~ s/^uniprot://;
    $curatedDesc =~ s/^(Sub|Rec)Name: Full=//;
    $curatedDesc =~ s/[{][A-Za-z0-9:|_. ;,-]+[}]//g;
    $curatedDesc =~ s/AltName:.*//;
    $curatedDesc =~ s/EC=/EC /g;
    $curatedDesc =~ s/ +;/;/g;
    $curatedDesc =~ s/;+ *$//;
    $charTitle = qq{$idShowHit has been studied experimentally.
                    It was manually added to GapMind's database and may not be in PaperBLAST.};
    $charTitle =~ s/\s+/ /g;
    $charLabel .= ", see " . a({-href => $stepDefURL, -title => $charTitle}, "rationale");
  }
  my $URL = "http://papers.genomics.lbl.gov/cgi-bin/litSearch.cgi?query=" . $first;
  my $link = a({-href => $URL, -title => "View $idShowHit in PaperBLAST"}, $curatedDesc);
  $link .= " (" . a({-title => $charTitle}, $charLabel) . ")";
  return $link;
}

sub CuratedToSeq($) {
  my ($curatedIds) = @_;
  if ($curatedIds =~ m/^curated2:(.*)$/) {
    my $protId = $1;
    my ($seq) = $dbhC->selectrow_array("SELECT seq FROM Curated2 WHERE protId = ?",
                                       {}, $protId);
    die "Unknown protId $protId in Curated2" unless $seq;
    return $seq;
  }
  #else
  if ($curatedIds =~ m/^uniprot:(.*)$/) {
    my $uniprotId = $1;
    my ($seq) = $dbhS->selectrow_array(qq{SELECT seq FROM StepQuery WHERE queryType = "uniprot" AND uniprotId = ?},
                                       {}, $uniprotId);
    die "Unknown uniprot id $uniprotId" unless $seq;
    return $seq;
  }
  #else
  my ($seq) = $dbhC->selectrow_array("SELECT seq FROM CuratedSeq WHERE curatedIds = ?",
                                     {}, $curatedIds);
  die "Unknown curatedIds $curatedIds" unless $seq;
  return $seq;
}

sub CandToOtherHTML($) {
  my ($cand) = @_;
  return ("", "", "") if $cand->{otherIds} eq "";

  my ($desc) = $dbhC->selectrow_array("SELECT descs from CuratedInfo WHERE curatedIds = ?",
                                      {}, $cand->{otherIds});
  die "Unknown curatedIds $cand->{otherIds}" unless defined $desc;
  $desc =~ s/;;.*//;
  my $first = $cand->{otherIds}; $first =~ s/,.*//;
  my $idShowOther = $first; $idShowOther =~ s/^.*://;
  my $linkOther = small(a({ -href => "http://papers.genomics.lbl.gov/cgi-bin/litSearch.cgi?query=$first",
                         -title => "view $idShowOther in PaperBLAST"},
                       $desc));
  my $otherIdentity = small(int(0.5 + $cand->{otherIdentity})."%");
  my $otherBits = small(sprintf("%.1f",$cand->{otherBits}));
  return ($linkOther, $otherIdentity, $otherBits);
}

sub LinkifyComment($) {
  my ($comment) = @_;
  my @words = split /\s/, $comment;
  my @out = ();
  foreach my $word (@words) {
    # Pull of leading brackets or parentheses
    my $pre = "";
    $pre = $1 if $word =~ m/^([\[({]+)/;
    $word =~ s/^([\[({]+)//;
    if ($word =~ m/^pfam:(PF\d+)/i) {
      my $pfam = $1;
      $word =~ s/^pfam:PF\d+//i;
      push @out, $pre . a({ -href => "hmmSearch.cgi?hmmId=$pfam" }, $pfam) . $word;
    } elsif ($word =~ m/^uniprot:([A-Z][A-Z0-9_]+)/i) {
      my $uniprotId = $1;
      $word =~ s/^uniprot:[A-Z][A-Z0-9_]+//i;
      push @out, $pre . a({ -href => "https://www.uniprot.org/uniprot/$uniprotId" }, $uniprotId) . $word;
    } elsif ($word =~ m/^pmid:(\d+)/i) {
      my $pmId = $1;
      $word =~ s/^pmid:\d+//i;
      push @out, $pre . a({ -href => "https://pubmed.ncbi.nlm.nih.gov/$pmId/" }, "PMID:$pmId") . $word;
    } elsif ($word =~ m/^(PMC\d+)/) {
      my $pmcId = $1;
      $word =~ s/^PMC\d+//;
      push @out, $pre . a({ -href => "http://www.ncbi.nlm.nih.gov/pmc/articles/" . lc($pmcId) . "/" },
                   $pmcId) . $word;
    } elsif ($word =~ m/^metacyc:([0-9A-Z][A-Z0-9-]+)/i) {
      my $metacycId = $1;
      $word =~ s/^metacyc:([0-9A-Z][A-Z0-9-]+)//i;
      push @out, $pre. a({ -href => "https://metacyc.org/META/NEW-IMAGE?object=$metacycId" },
                         "link") . $word;
    } elsif ($word =~ m!^URL:(http[A-Za-z0-9_,:./?&-]+)!i) {
      my $URL = $1;
      $word =~ s!^URL:(http[A-Za-z0-9_,:./?&-]+)!!;
      push @out, $pre. a({ -href => $URL }, "link") . $word;
    } elsif ($word =~ m/^EC:([0-9][.][0-9.]+[0-9])/i) {
      my $ec = $1;
      $word =~ s/^EC:([0-9][.][0-9.]+[0-9])//i;
      push @out, $pre . "EC " . a({ -href => "https://enzyme.expasy.org/EC/$ec" }, $ec) . $word;
    } else {
      push @out, $pre . $word;
    }
  }
  return join(" ", @out);
}

sub DataForStepParts($$) {
  my ($pathwayId, $stepId) = @_;
  my $stepParts = $dbhS->selectall_arrayref("SELECT * from StepPart WHERE pathwayId = ? AND stepId = ?",
                                            { Slice => {} }, $pathwayId, $stepId);
  my $stepQueries = $dbhS->selectall_arrayref("SELECT * from StepQuery WHERE pathwayId = ? AND stepId = ?",
                                              { Slice => {} }, $pathwayId, $stepId);
  my %curatedQuery = (); # curatedId (not ids -- a single compoennt) to stepquery row
  # (Includes entries of type curated or ignore)
  my %uniprotQuery = (); # uniprotId to stepquery row
  foreach my $sq (@$stepQueries) {
    if ($sq->{queryType} eq "curated" || $sq->{queryType} eq "ignore") {
      $sq->{desc} =~ s/;;/. /g;
      foreach my $id (split /,/, $sq->{curatedIds}) {
        $curatedQuery{$id} = $sq;
      }
    } elsif ($sq->{queryType} eq "uniprot") {
      $uniprotQuery{$sq->{uniprotId}} = $sq;
    }
  }
  return { 'uniprotQuery' => \%uniprotQuery,
           'curatedQuery' => \%curatedQuery };
}

sub FormatStepPart($$$) {
  my ($data, $stepPart, $orgId) = @_;
  die unless defined $data->{curatedQuery} && defined $data->{uniprotQuery};
  my $type = $stepPart->{partType};
  my $value = $stepPart->{value};

  if ($type eq "EC") {
    # Use local URLs for Curated BLAST links, instead of using the official papers.genomics.lbl.gov
    # site, because the genome may not exist at the public site
    my $out = "Curated proteins or TIGRFams with EC "
      . a({-href => "https://enzyme.expasy.org/EC/$value" }, $value);
    $out .= " (" . a({ -href => "genomeSearch.cgi?gdb="
                       . $orgs{$orgId}{gdb}
                       . "&gid=" . $orgs{$orgId}{gid}
                       . "&query=$value",
                       -title => "Run Curated BLAST" },
                     "search")
      . ")" if $orgId ne "";
    return $out;
  } elsif ($type eq "hmm") {
    return "HMM " . a({-href => HMMToURL($value) }, $value);
  } elsif ($type eq "term") {
    # nead to uri_escape the query because it may contain %
    my $URL = "curatedClusters.cgi?set=$set&word=1&query=" . uri_escape($value);
    $URL = "genomeSearch.cgi?gdb=" . $orgs{$orgId}{gdb}
          . "&gid=" . $orgs{$orgId}{gid}
          . "&word=1"
          . "&query=" . uri_escape($value)
            if $orgId ne "";
    return "Curated proteins matching "
      . a({ -href => $URL,
            -title => $orgId eq "" ? "" : "Run Curated BLAST" }, $value);
  } elsif ($type eq "curated") {
    my $URL = "http://papers.genomics.lbl.gov/cgi-bin/litSearch.cgi?query=".$value;
    my $show_id = $value; $show_id =~ s/^.*://;
    # Find the relevant step query
    return "Curated sequence " . a({-href => $URL, -title => "View in PaperBLAST"}, $show_id)
      . ": " . $data->{curatedQuery}{$value}{desc};
  } elsif ($type eq "uniprot") {
    my $URL = "https://www.uniprot.org/uniprot/".$value;
    return "UniProt sequence " . a({-href => $URL, -title => "View in UniProt"}, $value)
      . ": " . $data->{uniprotQuery}{$value}{desc};
  } elsif ($type eq "ignore_other") {
    my $URL = "http://papers.genomics.lbl.gov/cgi-bin/curatedSearch.cgi?word=1"
      . "&query=" . uri_escape($value); # value could contain %
    return "Ignore hits to items matching "
      . a({-href => $URL}, $value)
      . " when looking for 'other' hits";
  } elsif ($type eq "ignore") {
    my $URL = "http://papers.genomics.lbl.gov/cgi-bin/litSearch.cgi?query=$value";
    my $showId = $value; $showId =~ s/^.*://;
    return  "Ignore hits to "
      . a({-href => $URL, -title => "View in PaperBLAST"}, $showId)
      . " when looking for 'other' hits"
      . " (" . $data->{curatedQuery}{$value}{desc} . ")";
  } elsif ($type eq "ignore_hmm") {
    return "Do not include HMM " . a({ -href => HMMToURL($value) }, $value)
      . " (when considering EC numbers)";
  }
  die "Unknown StepPart type $type";
}

# orgId is optional
# If present, steps link to the candidates page
# If not present, steps link to #step (presumed to be present in show definition mode)
sub RulesToHTML($$$) {
  my ($stepsObj, $pathwayId, $orgId) = @_;
  my $stepScores = $dbhG->selectall_hashref("SELECT * from StepScore WHERE orgId = ? AND pathwayId = ?",
                                         "stepId",
                                         { Slice => {} }, $orgId, $pathwayId)
    if $orgId ne "";
  my $ruleScores = $dbhG->selectall_hashref("SELECT * from RuleScore WHERE orgId = ? AND pathwayId = ?",
                                            "ruleId",
                                            { Slice => {} }, $orgId, $pathwayId)
    if $orgId ne "";
  my $ruleOrdering = $dbhS->selectall_arrayref(qq{ SELECT ruleId, min(instanceId) AS minInstanceId FROM RuleInstance
                                                     WHERE pathwayId = ?
                                                     GROUP BY ruleId
                                                     ORDER BY minInstanceId },
                                               { Slice => {} }, $pathwayId);
  my @rulesInOrder = map $_->{ruleId}, @$ruleOrdering;
  my $out = start_ul;
  foreach my $ruleId (reverse @rulesInOrder) {
    my @instanceHTML = ();
    my $instances = $dbhS->selectcol_arrayref("SELECT instanceId from RuleInstance WHERE pathwayId = ? AND ruleId = ?",
                                              {}, $pathwayId, $ruleId);
    foreach my $instanceId (@$instances) {
      my $components = $dbhS->selectall_arrayref("SELECT * from InstanceComponent WHERE pathwayId = ? AND ruleId = ? AND instanceId = ?",
                                                 { Slice => {} },
                                                 $pathwayId, $ruleId, $instanceId);
      die $instanceId unless @$components > 0;
      my @parts;
      foreach my $component (@$components) {
        if ($component->{stepId} ne "") {
          if ($orgId eq "") {
            push @parts, a({ -title => $stepDesc{$pathwayId}{ $component->{stepId} },
                             -href => "#" . $component->{stepId} },
                              i($component->{stepId}));
          } else {
            push @parts, i(StepToShortHTML($component->{stepId}, $stepScores->{ $component->{stepId} }));
          }
        } else {
          my $subRuleId = $component->{subRuleId};
          die if $subRuleId eq "";
          my $param = { -title => "see rules for $subRuleId below" };
          $param->{style} = ScoreToStyle(RuleToMinScore($ruleScores->{$subRuleId}))
            if $orgId ne "";
          push @parts, span($param, $subRuleId);
        }
      }
      if (@parts == 1) {
        push @instanceHTML, $parts[0];
      } else {
        my $last = pop @parts;
        push @instanceHTML, join(", ", @parts) . " and $last";
      }
    }
    die $ruleId unless @instanceHTML > 0;
    my $ruleHTML = b("${ruleId}:");
    my $comment = "";
    $comment = "Comment: " . LinkifyComment($stepsObj->{ruleComments}{$ruleId})
      if exists $stepsObj->{ruleComments}{$ruleId}
        && $stepsObj->{ruleComments}{$ruleId} =~ m/\S/;
    if (@instanceHTML > 1) {
      foreach my $i (1..(scalar(@instanceHTML)-1)) {
        $instanceHTML[$i] = "or " . $instanceHTML[$i];
      }
      $out .= li($ruleHTML);
      $out .= start_ul;
      $out .= join("", map li($_), @instanceHTML);
      $out .= li($comment) if $comment ne "";
      $out .= end_ul;
    } else {
      $out .= li($ruleHTML, $instanceHTML[0]);
      $out .= start_ul . li($comment) . end_ul if $comment ne "";
    }
  }
  $out .= end_ul;
  return $out;
}

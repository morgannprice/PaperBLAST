#!/usr/bin/perl -w
# The second phase of the pipeline: identifying query terms
use strict;
use Getopt::Long;
use FindBin qw($Bin);

my @allsteps = qw{explodeam explodeeco parsepm pm pmc am mo refseq sprot byorg comb};
my $allsteps = join(",",@allsteps);

my $usage = <<END
run_terms.pl -in downloads -work work [ -test ] [ -steps $allsteps ]

This script build the list of query terms to search against
EuropePMC. It also explodes tar balls in the download directory, as
needed, and parses the pubmed and ecocyc downloads.
END
    ;

sub read_list( $ );
sub maybe_run( $ );

my ($test, $indir, $workdir);
my $dosteps = $allsteps;

die $usage
    unless GetOptions('in=s' => \$indir,
                      'work=s' => \$workdir,
                      'steps=s' => \$dosteps,
                      'test' => \$test)
    && @ARGV == 0
    && defined $indir && defined $workdir;
die "No such directory: $indir\n" unless -d $indir;
foreach my $subdir (qw{am oa pubmed refseq}) {
  die "No such directory: $indir/$subdir\n" unless -d "$indir/$subdir";
}
die "No such directory: $workdir\n" unless -d $workdir;

my @dosteps = split /,/, $dosteps;
my %dosteps = map { $_ => 1 } @dosteps;
my %allsteps = map { $_ => 1 } @allsteps;
foreach my $step (keys %dosteps) {
  die "Unrecognized step: $step\n" unless exists $allsteps{$step};
}

print STDERR "Test mode\n" if defined $test;

if (exists $dosteps{"explodeam"}) {
  print STDERR "Step: explode\n";
  my @files = &read_list("$indir/am/files");
  foreach my $gz (@files) {
    $gz =~ m/[.]xml[.]tar[.]gz$/ || die "Do not know how to handle $gz";
    &maybe_run("(cd $indir/am; tar xzf $gz)");
  }
}

if (exists $dosteps{"explodeeco"}) {
  my $ecotar = "$indir/ecoli.tar.gz";
  die "No such file: $ecotar\n" unless -e $ecotar;
  # Figure out what directory it will create
  print "STDERR listing $ecotar to find the target directory -- ignore any write error warning\n";
  open(FILES, "zcat $ecotar | tar --list |") || die "Error running tar on $ecotar";
  my $first = <FILES>;
  close(FILES);
  chomp $first;
  $first =~ s!/$!!;
  die "Error listing $ecotar" unless $first;

  &maybe_run("(cd $indir; tar xzf $ecotar)");
  &maybe_run("(ln -s ecocyc $first)");
}

if (exists $dosteps{"parsepm"}) {
  my @files1 = read_list("$indir/pubmed/baseline/files");
  my @files2 = read_list("$indir/pubmed/updatefiles/files");
  my @out = ();
  my $cmdsfile = "$workdir/parsepm.cmds";
  my $listfile = "$workdir/parsepm.list";
  open(CMDS, ">", $cmdsfile) || die "Cannot write to $cmdsfile";
  open(LIST, ">", $listfile) || die "Cannot write to $listfile";
  foreach my $file (@files1) {
    my $base = $file; $base =~ s/[.]xml[.]gz$//;
    my $out = "$workdir/pubmed1_$base.tab";
    print CMDS "zcat $indir/baseline/$base.xml.gz | $Bin/pubmedparse.pl > $out\n";
    print LIST "$out\n";
  }
  foreach my $file (@files2) {
    my $base = $file; $base =~ s/[.]xml[.]gz$//;
    my $out = "$workdir/pubmed2_$base.tab";
    print CMDS "zcat $indir/updatefiles/$base.xml.gz | $Bin/pubmedparse.pl > $out\n";
    print LIST "$out\n";
  }
  close(CMDS) || die "Error writing to $cmdsfile";
  close(LIST) || die "Error writing to $listfile";
  &maybe_run("$Bin/submitter.pl $cmdsfile");
}

if (exists $dosteps{"pm"}) {
  # Find locus tag like terms in the pubmed abstracts
  my $pmlist = "$workdir/parsepm.list";
  die "No such file: $pmlist\n" unless -e $pmlist;
  &maybe_run("$Bin/pubmed_words.pl -list $pmlist > $workdir/words.pubmed");
}

if (exists $dosteps{"pmc"}) {
  # Find locus tag like terms in open acess pmc articles
  my @files = &read_list("$indir/oa/files");
  my $cmdsfile = "$workdir/pmcwords.cmds";
  my $listfile = "$workdir/pmcwords.list";
  open(CMDS, ">", $cmdsfile) || die "Cannot write to $cmdsfile";
  open(LIST, ">", $listfile) || die "Cannot write to $listfile";
  foreach my $file (@files) {
    $file =~ m/[.]xml[.]gz$/ || die "Cannot parse file name $file";
    $file =~ s/[.]xml[.]gz$//;
    my $out = "$workdir/pmc_$file.words";
    print CMDS "gunzip -c $indir/oa/$file.xml.gz | $Bin/words.pl > $out\n";
    print LIST "$out\n";
  }
  close(CMDS) || die "Error writing to $cmdsfile";
  close(LIST) || die "Error writing to $listfile";
  &maybe_run("$Bin/submitter.pl $cmdsfile");
}

if (exists $dosteps{"am"}) {
  # Find locus like tags in author manuscript articles
  # Because of how exploding works out, need to use find to build the list
  unless (defined $test && -e "$indir/am/xml.list") {
    system("(cd $indir/am; find . -name '*.xml' -print > xml.list)") == 0
      || die "Error running find";
  }
  my @files = &read_list("$indir/am/xml.list");
  if (@files == 0) {
    print STDERR "Warning: no am files exploded yet\n";
  } else {
    print STDERR "Found " . scalar(@files) . " exploded am files\n";
  }
  # There are so many xml files -- to reduce the number of output files to deal with, group them into 200 sets
  my @sets = ();
  foreach my $i (0..(scalar(@files)-1)) {
    push @{ $sets[$i % 200] }, $files[$i];
  }
  my $cmdsfile = "$workdir/amwords.cmds";
  my $listfile = "$workdir/amwords.list";
  open(CMDS, ">", $cmdsfile) || die "Cannot write to $cmdsfile";
  open(LIST, ">", $listfile) || die "Cannot write to $listfile";
  foreach my $i (0..(scalar(@sets)-1)) {
    my $out = "$workdir/am_$i.words";
    my @in = map "$indir/am/$_", @{ $sets[$i] };
    print CMDS "cat " . join(" ", @in) . " | $Bin/words.pl > $out\n";
    print LIST "$out\n";
  }
  close(CMDS) || die "Error writing to $cmdsfile";
  close(LIST) || die "Error writing to $listfile";
  &maybe_run("$Bin/submitter.pl $cmdsfile");
}

my @wordsfiles = "$workdir/words.pubmed";
if (exists $dosteps{"mo"} || exists $dosteps{"refseq"} || exists $dosteps{"sprot"}) {
  push @wordsfiles, &read_list("$workdir/pmcwords.list");
  push @wordsfiles, &read_list("$workdir/amwords.list");
  &maybe_run("cat " . join(" ", @wordsfiles) . " > $workdir/words");
}

# Given words from pm/pmc/am, find identifiers in MicrobesOnline
if (exists $dosteps{"mo"}) {
  &maybe_run("$Bin/moIds.pl < $workdir/words > $workdir/mo.words");
  &maybe_run("$Bin/oaquery.pl < $workdir/mo.words > $workdir/mo.query");
}

# Given words from pm/pmc/am, find identifiers in RefSeq
if (exists $dosteps{"refseq"}) {
  my @files = &read_list("$indir/refseq/files");
  @files = grep m/genomic.gbff.gz$/, @files;
  print STDERR "Found " . scalar(@files) . " genomic RefSeq files\n";

  my $cmdsfile = "$workdir/refseqwords.cmds";
  my $listfile = "$workdir/refseqwords.list";
  open(CMDS, ">", $cmdsfile) || die "Cannot write to $cmdsfile";
  open(LIST, ">", $listfile) || die "Cannot write to $listfile";
  foreach my $file (@files) {
    $file =~ m/[.]gbff[.]gz$/ || die "Cannot parse file name $file";
    $file =~ s/[.]gbff[.]gz$//;
    my $out = "$workdir/refseq_$file.words";
    print CMDS "gunzip -c $indir/refseq/$file.gbff.gz | $Bin/findRefSeqQueries.pl $workdir/words > $out\n";
    print LIST "$out\n";
  }
  close(CMDS) || die "Error writing to $cmdsfile";
  close(LIST) || die "Error writing to $listfile";
  &maybe_run("$Bin/submitter.pl $cmdsfile");
  my @out = &read_list($listfile);
  &maybe_run("cat " . join(" ", @out) . " | $Bin/convertRefSeqQueries.pl > refseq.query");
}

# Given words from pm/pmc/am, find identifiers matching UniProt
if (exists $dosteps{"sprot"}) {
  &maybe_run("$Bin/select_uniprot_words.pl < $workdir/words > $workdir/words.uniprot");
  my $in = "$indir/uniprot_sprot.fasta.gz";
  die "No such file: $in\n" unless -e $in;
  &maybe_run("zcat $in | $Bin/sprotToQuery.pl -words $workdir/words.uniprot > sprot.query");
}

if (exists $dosteps{"byorg"}) {
  &maybe_run("$Bin/queryProtByOrg.pl < $Bin/../static/top_taxa > $workdir/byorg.query");
}

if (exists $dosteps{"comb"}) {
  my @in = map "$workdir/$_.query", qw{mo refseq sprot byorg};
  &maybe_run("$Bin/removeDupQueries.pl " . join(" ", @in) . " > $workdir/comb.query");
}

sub read_list($) {
  my ($file) = @_;
  open(LIST, "<", $file) || die "Cannot read $file";
    my @lines = <LIST>;
    close(LIST) || die "Error reading $file";
    @lines = map { chomp; $_ } @lines;
    return(@lines);
}

sub maybe_run($) {
    my ($cmd) = @_;
    if (defined $test) {
        print STDERR "Would run\t$cmd\n";
    } else {
        print STDERR "Running $cmd\n";
        system($cmd) == 0
            || die "Error running $cmd: $!\n";
    }
}

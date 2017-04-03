#!/usr/bin/perl -w
# The second phase of the pipeline: identifying query terms
use strict;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use pbutils;

my @allsteps = qw{explodeam explodeeco parsepm pm pmc am mo generif refseq sprot byorg comb stats};
my $allsteps = join(",",@allsteps);

my $usage = <<END
run_terms.pl -in downloads [ -work work ] [ -test ] [ -steps ... ]

Given the downloads directory (populated by download.pl), build the
list of query terms to search against EuropePMC.

The steps are:
explodeam: explode the author manuscript tarballs
explodeeco: explode the EcoCyc tarball
parsepm: parse the abstracts from PubMed
pm: identify query terms in PubMed abstracts
pmc: identify query terms in open access papers from EuropePMC
am: identify query terms in author manuscript papers from EuropePMC
mo: given the query terms, find locus tags in MicrobesOnline
generif: convert gene ids to refseq protein ids
refseq: given the query terms, find locus tags or refseq ids in RefSeq
sprot: given the query terms, find UniProt accessions or entry_names
byorg: given the top genomes (see static/top_taxa), identify query terms
comb: combine all the queries

You can use the -steps argument with a comma-delimited subset of steps. The default is:
	-steps $allsteps
END
    ;

my ($test, $indir);
my $workdir = "work";
my $dosteps = $allsteps;

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
  print STDERR "Listing $ecotar to find the target directory -- ignore any write error warning\n";
  open(FILES, "zcat $ecotar | tar --list |") || die "Error running tar on $ecotar";
  my $first = <FILES>;
  close(FILES);
  chomp $first;
  $first =~ s!/$!!;
  die "Error listing $ecotar" unless $first;

  &maybe_run("(cd $indir; tar xzf ecoli.tar.gz)");
  &maybe_run("(cd $indir; rm ecocyc >& /dev/null; ln -s $first ecocyc)");
  if (! defined $test) {
    foreach my $file ("$indir/ecocyc/data/proteins.dat", "$indir/ecocyc/data/protseq.fsa") {
      die "Cannot find ecocyc file $file\n" unless -e $file;
    }
  }
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
    print CMDS "zcat $indir/pubmed/baseline/$base.xml.gz | $Bin/pubmedparse.pl > $out\n";
    print LIST "$out\n";
  }
  foreach my $file (@files2) {
    my $base = $file; $base =~ s/[.]xml[.]gz$//;
    my $out = "$workdir/pubmed2_$base.tab";
    print CMDS "zcat $indir/pubmed/updatefiles/$base.xml.gz | $Bin/pubmedparse.pl > $out\n";
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
  &maybe_run("cut -f 2 $workdir/words | sort -u > $workdir/words.u");
}

# Given words from pm/pmc/am, find identifiers in MicrobesOnline
if (exists $dosteps{"mo"}) {
  &maybe_run("$Bin/moIds.pl < $workdir/words > $workdir/mo.words");
  &maybe_run("$Bin/oaquery.pl < $workdir/mo.words > $workdir/mo.query");
}

if (exists $dosteps{"generif"}) {
  die "No such file: $indir/generifs_basic" unless -e "$indir/generifs_basic";
  # Note -- these protein ids are without the version number
  &maybe_run("cut -f 2 $indir/generifs_basic | sort -u | $Bin/geneIdToProtein.pl > $workdir/generifs_prot");
}

# Given words from pm/pmc/am/, find identifiers in RefSeq
if (exists $dosteps{"refseq"}) {
  # get just the protein ids from the pmc links
  &maybe_run("cut -d, -f 1 ind/RefSeq_PMC.csv | grep 'P_' | sed -e 's/[.][0-9]*//' > $workdir/words.refseq.pmclinks");
  &maybe_run("(cut -f 2 $workdir/generifs_prot; cat $workdir/words.u; cat $workdir/words.refseq.pmclinks) | sort -u > $workdir/words.refseq");

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
    print CMDS "gunzip -c $indir/refseq/$file.gbff.gz | $Bin/findRefSeqQueries.pl $workdir/words.refseq > $out\n";
    print LIST "$out\n";
  }
  close(CMDS) || die "Error writing to $cmdsfile";
  close(LIST) || die "Error writing to $listfile";
  &maybe_run("$Bin/submitter.pl $cmdsfile");
  my @out = &read_list($listfile);
  &maybe_run("(cat " . join(" ", @out) . " | $Bin/convertRefSeqQueries.pl > $workdir/refseq.query) >& $workdir/refseq.query.log");
}

# Given words from pm/pmc/am, find identifiers matching UniProt
if (exists $dosteps{"sprot"}) {
  &maybe_run("$Bin/select_uniprot_words.pl < $workdir/words > $workdir/words.uniprot1");
  &maybe_run("cut -d, -f 1 < $indir/UniProt_PMC.csv | grep -v UniProt | sort -u > $workdir/words.uniprot2");
  &maybe_run("sort -u $workdir/words.uniprot1 $workdir/words.uniprot2 > $workdir/words.uniprot");

  my $in = "$indir/uniprot_sprot.fasta.gz";
  die "No such file: $in\n" unless -e $in;
  my $in2 = "$indir/uniprot_trembl.fasta.gz";
  die "No such file: $in2\n" unless -e $in2;
  &maybe_run("zcat $in $in2 | $Bin/sprotToQuery.pl -words $workdir/words.uniprot > $workdir/sprot.query");
}

if (exists $dosteps{"byorg"}) {
  &maybe_run("$Bin/queryProtByOrg.pl < $Bin/../static/top_taxa > $workdir/byorg.query");
}

if (exists $dosteps{"comb"}) {
  my @in = map "$workdir/$_.query", qw{mo refseq sprot byorg};
  &maybe_run("($Bin/uniqueQueries.pl " . join(" ", @in) . " > $workdir/comb.query) >& $workdir/comb.query.log");
}

if (exists $dosteps{"stats"}) {
  if (! defined $test) {
    print STDERR join("\t", "file", "thousands\n");
    foreach my $file (qw{words.u byorg.query mo.query refseq.query sprot.query comb.query}) {
      my $nlines = `wc -l < $workdir/$file`;
      chomp $nlines;
      print STDERR join("\t", $file, sprintf("%.1f", $nlines/1000))."\n";
    }
  }
}

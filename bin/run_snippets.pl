#!/usr/bin/perl -w
# The fourth phase of the pipeline: identifying snippets
use strict;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use pbutils;

my @allsteps = qw{parse pmclinks oa am elsevier crossref pubmed comb generif stats};
my $allsteps = join(",", @allsteps);
my $usage = <<END
run_snippets.pl -in downloads [ -work work ] [ -test ]
	[ -cache cache ] [ -cache-only ]
	[ -steps $allsteps ]

Given the EuropePMC results in work/epmc, identify
additional links from EuropePMC, and compute snippets using:
  The open access manuscripts in downloads/
  The author manuscripts in downloads/
  Full text fetched with the Elsevier API
	[skipped if in cache only mode]
  Full text fetched with CrossRef API,
	and/or using the files in cache/
  Pubmed abstracts in work/pubmed*.tab

  (The elsevier key and the crossref token must be in the cache directory
  unless this is cache-only mode.)

Also parses GeneRIF.
END
;

my ($indir, $test, $cacheOnly);
my $workdir = "work";
my $cachedir = "cache";
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
                    'test' => \$test,
                    'cache=s' => \$cachedir,
                    'cache-only' => \$cacheOnly,
                   'steps=s' => \$dosteps)
  && @ARGV == 0
  && defined $indir && defined $workdir && defined $cachedir;
die "No such directory: $indir\n" unless -d $indir;
die "No such directory: $workdir\n" unless -d $workdir;
die "No such directory: $cachedir\n" unless -d $cachedir;

my @dosteps = split /,/, $dosteps;
my %dosteps = map { $_ => 1 } @dosteps;
my %allsteps = map { $_ => 1 } @allsteps;
foreach my $step (keys %dosteps) {
  die "Unrecognized step: $step\n" unless exists $allsteps{$step};
}

print STDERR "Test mode\n" if defined $test;
print STDERR "Cache only\n" if defined $cacheOnly;

foreach my $file (qw{comb.query epmc}) {
  die "No such file: $workdir/$file\n" unless -e "$workdir/$file";
}

my $elsevier_key = "$cachedir/elsevier_key";
die "Elsevier key not find: no such file: $elsevier_key\n"
  unless -e $elsevier_key || defined $cacheOnly;

my $crossref_token = "$cachedir/crossref_token";
die "Crossref token not fund: no such file: $crossref_token\n"
  unless -e $crossref_token || defined $cacheOnly;

if (exists $dosteps{"parse"}) {
  &maybe_run("$Bin/parseEuropePMCHits.pl -in $workdir/comb.query -hits $workdir/epmc -out $workdir/hits >& $workdir/epmc_parse.log");
}

if (exists $dosteps{"pmclinks"}) {
  &maybe_run("$Bin/addPMCLinks.pl -query $workdir/comb.query -papers $workdir/hits.papers -in $indir -out $workdir/pmclinks >& $workdir/addpmclinks.log");
}

if (exists $dosteps{"oa"}) {
  my @files = &read_list("$indir/oa/files");
  my @in = map "$indir/oa/$_", @files;
  &write_list(\@in, "$workdir/snippets.oa.list");
  &maybe_run("$Bin/buildSnippets.pl -list $workdir/snippets.oa.list -out $workdir/snippets_oa $workdir/hits.papers $workdir/pmclinks.papers >& $workdir/snippets_oa.log");
}

if (exists $dosteps{"am"}) {
  my @files = &read_list("$indir/am/xml.list");
  my @in = map "$indir/am/$_", @files;
  &write_list(\@in, "$workdir/snippets.am.list");
  &maybe_run("$Bin/buildSnippets.pl -list $workdir/snippets.am.list -out $workdir/snippets_am $workdir/hits.papers $workdir/pmclinks.papers >& $workdir/snippets_am.log");
}

if (exists $dosteps{"elsevier"} && !defined $cacheOnly) {
  &maybe_run("$Bin/elsevierFetch.pl -key $elsevier_key -dir $cachedir -papers $workdir/hits.papers -journals $Bin/../static/elsevier_journals >& $workdir/snippets_elsevier.log");
}

if (exists $dosteps{"crossref"}) {
  my $cache_only_arg = defined $cacheOnly ? "-cache-only" : "";
  &maybe_run("$Bin/crossrefSnippets.pl $cache_only_arg -token $crossref_token -papers $workdir/hits.papers -dir $cachedir -out $workdir/snippets_crossref >& $workdir/snippets_crossref.log");
}

if (exists $dosteps{"pubmed"}) {
  &maybe_run("$Bin/abstractSnippets.pl -list $workdir/parsepm.list -out $workdir/snippets_pubmed $workdir/hits.papers");
}

if (exists $dosteps{"comb"}) {
  &maybe_run("$Bin/combineSnippets.pl -out $workdir/snippets_comb $workdir/snippets_oa $workdir/snippets_am $workdir/snippets_crossref $workdir/snippets_pubmed  >& $workdir/snippets_comb.log");
}

if (exists $dosteps{"generif"}) {
  &maybe_run("$Bin/join.pl -header 0 -match 1.2=2.1 -keep 1.3 $indir/generifs_basic $workdir/generifs_prot | sed -e 's/,/\\n/g' | sort -u > $workdir/generifs_pmid");
  my @files1 = read_list("$indir/pubmed/baseline/files");
  my @files2 = read_list("$indir/pubmed/updatefiles/files");
  my $cmdsfile = "$workdir/parsepm_generif.cmds";
  my $listfile = "$workdir/parsepm_generif.list";
  open(CMDS, ">", $cmdsfile) || die "Cannot write to $cmdsfile";
  open(LIST, ">", $listfile) || die "Cannot write to $listfile";
  foreach my $file (@files1) {
    my $base = $file; $base =~ s/[.]xml[.]gz$//;
    my $out = "$workdir/pubmed1_generif_$base.tab";
    print CMDS "$Bin/pubmedFields.pl -gz $indir/pubmed/baseline/$base.xml.gz < $workdir/generifs_pmid > $out\n";
    print LIST "$out\n";
  }
  foreach my $file (@files2) {
    my $base = $file; $base =~ s/[.]xml[.]gz$//;
    my $out = "$workdir/pubmed2_generif_$base.tab";
    print CMDS "$Bin/pubmedFields.pl -gz $indir/pubmed/updatefiles/$base.xml.gz < $workdir/generifs_pmid > $out\n";
    print LIST "$out\n";
  }
  close(CMDS) || die "Error writing to $cmdsfile";
  close(LIST) || die "Error writing to $listfile";
  &maybe_run("$Bin/submitter.pl $cmdsfile");
  &maybe_run("xargs cat < $listfile > $workdir/generifs_pmid.tab");
  # Note -- we give generifTables *all* of the potential refseq queries to work with.
  # These should virtually all have hits, either from GeneRIF or from searching EuropePMC, but
  # there could be exceptions. This will lead to unnecessary entries in generif_tab.queries
  # buildLitDb.pl will remove these from the final FAA file, so they will not show up as BLAST hits,
  # and they will not be in the GenePaper table, but, they could be in the final Gene table.
  &maybe_run("$Bin/generifTables.pl -rif $indir/generifs_basic -prot $workdir/generifs_prot -query $workdir/refseq.query -known $workdir/hits.queries -papers $workdir/generifs_pmid.tab -out $workdir/generif_tab >& $workdir/generifTables.log");
}

if (exists $dosteps{"stats"} && !defined $test) {
  print STDERR join("\t", "file", "thousands\n");
  foreach my $file (qw{hits.queries hits.papers pmclinks.queries pmclinks.papers snippets_oa snippets_am snippets_crossref snippets_pubmed snippets_comb generif_tab.papers generif_tab.rif}) {
    my $nlines = `wc -l < $workdir/$file`;
    chomp $nlines;
    print STDERR join("\t", $file, sprintf("%.1f", $nlines/1000))."\n";
  }
  my $nFull = `grep -c full $workdir/snippets_comb.access`;
  chomp $nFull;
  my $nAbstract = `grep -c abstract $workdir/snippets_comb.access`;
  chomp $nAbstract;
  print STDERR sprintf("By paper: full-text %.1fK abstract-only %.1fK\n", $nFull/1000.0, $nAbstract/1000.0);
}

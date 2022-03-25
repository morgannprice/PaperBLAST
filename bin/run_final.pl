#!/usr/bin/perl -w
# The last phase of the pipeline: building the database
use strict;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use pbutils;
use DBI;

my @allsteps = qw{sprot ecocyc biolip db stats compare};
my $allsteps = join(",",@allsteps);
my $comparedir = "$Bin/../data";

my $usage = <<END
run_final.pl -in downloads [ -work work ] [ -test ]
	[ -steps $allsteps ]
	[ -compare $comparedir ]

Given the downloads and work directories, and the output of
run_snippets.pl, builds a database in work/data, and compares it to
the existing database (from $comparedir)

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
                      'test' => \$test,
                      'compare=s' => \$comparedir)
    && @ARGV == 0
    && defined $indir && defined $workdir;
die "No such directory: $indir\n" unless -d $indir;
die "No such directory: $workdir\n" unless -d $workdir;
die "No such directory: $comparedir\n" unless -d $comparedir;

my @dosteps = split /,/, $dosteps;
my %dosteps = map { $_ => 1 } @dosteps;
my %allsteps = map { $_ => 1 } @allsteps;
foreach my $step (keys %dosteps) {
    die "Unrecognized step: $step\n" unless exists $allsteps{$step};
}

print STDERR "Test mode\n" if defined $test;

my $outdir = "$workdir/data";
&mkdir_if_needed($outdir);

if (exists $dosteps{sprot}) {
  my $sprot = "$indir/uniprot_sprot.dat.gz";
  die "No such file: $sprot\n" unless -e $sprot;
  &maybe_run("zcat $sprot | bin/sprotCharacterized.pl > $workdir/sprot.curated_parsed");
  &maybe_run("zcat $sprot | bin/parse_SwissProt_features.pl > $workdir/sprotFt.tab");
}

if (exists $dosteps{ecocyc}) {
  my $datfile = "$indir/ecocyc/data/proteins.dat";
  die "No such file: $datfile\n" unless -e $datfile;
  my $fsafile = "$indir/ecocyc/data/protseq.fsa";
  die "No such file: $fsafile\n" unless -e $fsafile;
  &maybe_run("bin/parse_ecocyc.pl $indir/ecocyc/data > $workdir/ecocyc.curated_parsed");
}

if (exists $dosteps{biolip}) {
  # build work/PDBLigands.tab and static/biolip.curated_parsed
  &maybe_run("zcat $indir/components.cif.gz | $Bin/parsePdbCdd.pl > $workdir/PDBLigands.tab");
  my @biolip = ("$indir/BioLiP_2013-03-6_nr.txt", "$indir/BioLiP_UP_nr.txt");
  my $protNames = "$indir/protnames.lst";
  foreach my $file ($protNames, @biolip) {
    die "No such file: $file\n" unless -e $file || defined $test;
  }
  &maybe_run("$Bin/biolipCurated.pl -biolip @biolip -pdbnames $protNames -pdbligands $workdir/PDBLigands.tab > static/biolip.curated_parsed");
}

if (exists $dosteps{db}) {
  foreach my $file (qw{snippets_comb hits.papers pmclinks.papers generif_tab.rif
                       sprot.curated_parsed ecocyc.curated_parsed
                       PDBLigands.tab}) {
    die "No such file: $workdir/$file" unless -e "$workdir/$file" || defined $test;
  }
  my @cmd = ();
  @cmd = ("$Bin/buildSiteTables.pl",
          "-sprot", "$workdir/sprotFt.tab",
          "-BioLiP", "$indir/BioLiP_2013-03-6_nr.txt", "$indir/BioLiP_UP_nr.txt",
          "-seqres", "$indir/pdb_seqres.txt",
          "-sprotFasta", "$indir/uniprot_sprot.fasta.gz",
          "-pdbnames", "$indir/protnames.lst",
          "-outdir", $workdir);
  &maybe_run(join(" ", @cmd));
  @cmd = ("$Bin/buildLitDb.pl",
          "-dir", $outdir,
          "-snippets", "$workdir/snippets_comb",
          "-rif", "$workdir/generif_tab.rif",
          "-curated",
          "$workdir/sprot.curated_parsed", "$workdir/ecocyc.curated_parsed",
          "static/CAZy.curated_parsed", "static/CharProtDB.curated_parsed",
          "static/metacyc.curated_parsed", "static/reanno.curated_parsed",
          "static/REBASE.curated_parsed", "static/BRENDA.curated_parsed",
          "static/TCDB.curated_parsed", "static/biolip.curated_parsed",
          "-prefix",
          "$workdir/hits",
          "$workdir/pmclinks",
          "$workdir/generif_tab");
  &maybe_run(join(" ", @cmd));
  &maybe_run("$Bin/loadSiteTables.pl -db $outdir/litsearch.db -indir $workdir");
  &maybe_run("cp $workdir/hassites.faa $outdir/");
  &maybe_run("$Bin/blast/formatdb -p T -o T -i $outdir/hassites.faa");
}

if (exists $dosteps{stats}) {
  if (!defined $test) {
    my $sqldb = "$outdir/litsearch.db";
    my $dbh = DBI->connect("dbi:SQLite:dbname=$sqldb","","",{ RaiseError => 1 }) || die $DBI::errstr;
    my $pm1 = $dbh->selectcol_arrayref(qq{ SELECT DISTINCT pmId FROM GenePaper WHERE pmId <> "" });
    my $pm2 = $dbh->selectcol_arrayref(qq{ SELECT DISTINCT pmId FROM GeneRIF });
    my $pm3 = $dbh->selectcol_arrayref(qq{ SELECT DISTINCT pmId FROM CuratedPaper });
    my $pmc = $dbh->selectcol_arrayref(qq{ SELECT DISTINCT pmcId FROM GenePaper
                                           WHERE pmId = "" AND pmcId <> "" });
    my %pm = map { $_ =>  1 } @$pm1;
    foreach my $pm (@$pm2) { $pm{$pm} = 1; }
    foreach my $pm (@$pm3) { $pm{$pm} = 1; }
    my $nPaper = scalar(keys %pm) + scalar(@$pmc);

    my $nSeq = `grep -c ">" $outdir/uniq.faa`;
    chomp $nSeq;
    my $date = `date -r $workdir/epmc '+%B %-d %Y'`;
    chomp $date;
    my %stats = ( "nSeq" => $nSeq, "nPaper" => $nPaper, "date" => $date );
    open(STATS, ">", "$outdir/stats") || die "Error writing to $outdir/stats";
    while (my ($key,$value) = each %stats) {
      print STATS "$key\t$value\n";
    }
    close(STATS) || die "Error writing to $outdir/stats\n";
    $dbh->disconnect();
  }
  &maybe_run("cat $outdir/stats");
}

if (exists $dosteps{compare}) {
  my $newdb = "$outdir/litsearch.db";
  my $olddb = "$comparedir/litsearch.db";
  if (! -e $olddb) {
    print STDERR "No such file: $olddb -- skipping the comparison step\n";
  } else {
    &maybe_run("$Bin/compare_dbs.pl -old $olddb -new $newdb -out $workdir/dbcmp");
  }
}

#!/usr/bin/perl -w
# The first phase of the pipeline: downloading
use strict;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use pbutils;

my @allsteps = qw{ecocyc oa am pmclinks refseq generif pubmed uniprot pdb biolip};
my $dosteps = join(",", @allsteps);

my $usage = <<END
download.pl -dir downloads [ -steps $dosteps ] [ -test ]
	[ -ecocyc ecoli_tar_gz_or_URL ]

This script downloads all of the inputs for building a PaperBLAST
database. These are all put into the specified directory. These inputs
are:

From EuropePMC:
The open access manuscripts, into dir/oa/*.xml.gz (21 GB)
	These are listed in dir/oa/files
The author manuscripts, into dir/am/*.tar.gz (8.5 GB)
	These are then exploded to give dir/am/*/*.xml
	The directory names are listed in dir/am/files
Their links of PMC ids to RefSeq or UniProt, into
	dir/uniprot.csv and dir/refseq.csv

From RefSeq:
The compressed genbank format files, into dir/refseq/complete.*.gbff.gz (86 GB)
	These are listed in dir/refseq/files

From GeneRIF:
ftp://ftp.ncbi.nih.gov/gene/GeneRIF/generifs_basic.gz

From PubMed: metadata about articles (used primarily for finding snippets in abstracts)
	These are placed within dir/pubmed/updatefiles/*.xml.gz
        or dir/pubmed/baseline/*.xml.gz
	and are listed in dir/pubmed/updatefiles/files and dir/pubmed/baseline/files
	dir/ftp.ncbi.nlm.nih.gov/pubmed/updatefiles/*.xml.gz
	dir/ftp.ncbi.nlm.nih.gov/pubmed/updatefiles/*.xml.gz

From UniProt:
SwissProt (curated) annotations: dir/uniprot_sprot.dat.gz (0.5 GB)
SwissProt (curated) sequences: dir/uniprot_sprot.fasta.gz (0.1 GB)
TReMBL (non-curated) sequences: dir/uniprot_trembl.fasta.gz (16 GB)

From EcoCyc: ecoli.tar.gz (158 MB)
	This is then exploded to give the dir/ecocyc/version.number subdirectory,
	and then a symlink is created for dir/ecocyc_latest/
	Unfortunately, there is no standard public URL for downloading EcoCyc.
	You need to contact them to obtain a URL and password.
	You can specify a URL with the -ecocyc argument, with something like
	-ecocyc http://biocyc-flatfiles:data-12345\@bioinformatics.ai.sri.com/ecocyc/dist/flatfiles-12345678/ecoli.tar.gz
	or you can put the URL in the file ecocyc.URL
	or you can download the tar ball manually into ecoli.tar.gz

PDB metadata (300 MB):
  http://www.ebi.ac.uk/thornton-srv/databases/pdbsum/data/protnames.lst
  ftp://ftp.wwpdb.org/pub/pdb/data/monomers/components.cif.gz

BioLiP annotation files (364 MB):
  https://zhanglab.ccmb.med.umich.edu/BioLiP/download/BioLiP_nr.tar.bz2
  (extract BioLiP_2013-03-6_nr.txt)
  BioLiP_updated_set/*txt
  (see download_BioLip_updates.pl)

(Most sizes for downloads are as of January 2017; PDB/BioLip sizes are from December 2019)

You can set the environment variable PB_DOWNLOAD_PASS = 1 to avoid
re-downloading files that already exist (as non-empty files) due to a
previous partially- successful download.
END
    ;

my $test;
my $ecocyc_URL;

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

sub maybe_wget($$) {
    my ($url, $file) = @_;
    if (defined $test) {
        print STDERR "Would wget $url into $file\n";
    } else {
        &wget($url, $file);
    }
}

my $dir;
die $usage
    unless GetOptions('dir=s' => \$dir, 'steps=s' => \$dosteps, 'test' => \$test,
                     'ecocyc=s' => \$ecocyc_URL)
    && @ARGV == 0;
die $usage unless defined $dir;
die "No such directory: $dir\n" unless -d $dir;

my @dosteps = split /,/, $dosteps;
my %dosteps = map { $_ => 1 } @dosteps;
my %allsteps = map { $_ => 1 } @allsteps;
foreach my $step (keys %dosteps) {
    die "Unrecognized step: $step\n" unless exists $allsteps{$step};
}

print STDERR "Test mode\n" if defined $test;
print STDERR "Not replacing already-existing non-empty files\n"
  if $ENV{PB_DOWNLOAD_PASS};

my $listfile = "$dir/listing.$$";

if (exists $dosteps{"ecocyc"}) {
    print STDERR "Step ecocyc\n";
    if (!defined $ecocyc_URL) {
      if (-e "ecocyc.URL") {
        open(URL, "<", "ecocyc.URL") || die "Cannot read ecocyc.URL";
        $ecocyc_URL = <URL>;
        chomp $ecocyc_URL;
        close(URL) || die "Error reading ecocyc.URL";
        die "Invalid URL in $ecocyc_URL" unless $ecocyc_URL =~ m/^http|ftp/
          || -e $ecocyc_URL;
        print STDERR "URL for ecocyc: $ecocyc_URL\n";
      }
    }
    $ecocyc_URL = "ecoli.tar.gz" if !defined $ecocyc_URL;
    if (-e $ecocyc_URL) {
      print STDERR "Using manually downloaded tarball $ecocyc_URL for ecocyc\n";
      &maybe_run("cp $ecocyc_URL $dir/");
    } else {
      &maybe_wget($ecocyc_URL, "$dir/ecoli.tar.gz");
    }
}

if (exists $dosteps{"oa"}) {
    print STDERR "Step oa\n";
    &mkdir_if_needed("$dir/oa");
    unlink($listfile); # in case it exists -- wget will not overwrite if PB_DOWNLOAD_PASS is set
    &wget("http://europepmc.org/ftp/oa/", $listfile);
    my @files = &ftp_html_to_files($listfile);
    @files = grep m/[.]xml[.]gz$/, @files;
    die "No .xml.gz files in oa, see $listfile" if @files == 0;
    print STDERR "Found " . scalar(@files) . " oa gz files to fetch\n";
    &write_list(\@files, "$dir/oa/files");
    foreach my $file (@files) {
        &maybe_wget("http://europepmc.org/ftp/oa/$file", "$dir/oa/$file");
    }
}

if (exists $dosteps{"am"}) {
    print STDERR "Step am\n";
    &mkdir_if_needed("$dir/am");
    unlink($listfile);
    &wget("http://europepmc.org/ftp/manuscripts/", $listfile);
    my @files = &ftp_html_to_files($listfile);
    @files = grep m/[.]xml[.]tar[.]gz$/, @files;
    die "No xml.tar.gz files in am, see $listfile" if @files == 0;
    print STDERR "Found " . scalar(@files) . " am gz files to fetch\n";
    &write_list(\@files, "$dir/am/files");
    foreach my $file (@files) {
        &maybe_wget("http://europepmc.org/ftp/manuscripts/$file", "$dir/am/$file");
    }
}

if (exists $dosteps{"pmclinks"}) {
  print STDERR "Step pmclinks\n";
  foreach my $file (qw{uniprot.csv refseq.csv}) {
    &maybe_wget("ftp://ftp.ebi.ac.uk/pub/databases/pmc/TextMinedTerms/$file", "$dir/$file");
  }
}

if (exists $dosteps{"refseq"}) {
    print STDERR "Step refseq\n";
    &mkdir_if_needed("$dir/refseq");
    unlink($listfile);
    &wget("ftp://ftp.ncbi.nlm.nih.gov/refseq/release/complete/", $listfile);
    my @files = &ftp_html_to_files($listfile);
    @files = grep m/^complete.*gbff[.]gz$/, @files;
    die "No complete*.gbff.gz files in refseq, see $listfile" if @files == 0;
    print STDERR "Found " . scalar(@files) . " refseq gbff.gz files to fetch\n";
    &write_list(\@files, "$dir/refseq/files");
    foreach my $file (@files) {
        &maybe_wget("ftp://ftp.ncbi.nlm.nih.gov/refseq/release/complete/$file", "$dir/refseq/$file");
    }
}

if (exists $dosteps{"generif"}) {
  print STDERR "Step generif\n";
  &maybe_wget("ftp://ftp.ncbi.nih.gov/gene/GeneRIF/generifs_basic.gz", "$dir/generifs_basic.gz");
  &maybe_run("gunzip $dir/generifs_basic.gz");
}

if (exists $dosteps{"pubmed"}) {
    print STDERR "Step pubmed\n";
    &mkdir_if_needed("$dir/pubmed");
    foreach my $part (qw{baseline updatefiles}) {
        &mkdir_if_needed("$dir/pubmed/$part");
        unlink($listfile);
        &wget("ftp://ftp.ncbi.nlm.nih.gov/pubmed/$part/", $listfile);
        my @files = &ftp_html_to_files($listfile);
        @files = grep m/[.]xml[.]gz$/, @files;
        die "No *.xml.gz files in pubmed $part, see $listfile" if @files == 0;
        print STDERR "Found " . scalar(@files) . " pubmed $part xml.gz files to fetch\n";
        &write_list(\@files, "$dir/pubmed/$part/files");
        foreach my $file (@files) {
            &maybe_wget("ftp://ftp.ncbi.nlm.nih.gov/pubmed/$part/$file", "$dir/pubmed/$part/$file");
        }
    }
}

if (exists $dosteps{"uniprot"}) {
    print STDERR "Step uniprot\n";
    foreach my $file (qw{uniprot_sprot.dat.gz uniprot_sprot.fasta.gz uniprot_trembl.fasta.gz}) {
        &maybe_wget("ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/$file",
                    "$dir/$file");
    }
}

if (exists $dosteps{"pdb"}) {
  &maybe_wget("http://www.ebi.ac.uk/thornton-srv/databases/pdbsum/data/protnames.lst",
              "$dir/protnames.lst");
  &maybe_wget("ftp://ftp.wwpdb.org/pub/pdb/data/monomers/components.cif.gz",
              "$dir/components.cif.gz");
  &maybe_wget("ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz",
              "$dir/pdb_seqres.txt.gz");
  &maybe_run("gunzip $dir/pdb_seqres.txt.gz");
}

if (exists $dosteps{biolip}) {
  my $biolipCacheDir = "BioLiP_updated_set";
  mkdir($biolipCacheDir);
  die "Not a directory: $biolipCacheDir\n" unless -d $biolipCacheDir;
  my $mainFile = "$dir/BioLiP_2013-03-6_nr.txt";
  unless (-s $mainFile) {
    my $tarfile = "$dir/BioLiP_nr.tar";
    &maybe_wget("https://zhanglab.ccmb.med.umich.edu/BioLiP/download/BioLiP_nr.tar.bz2",
                "$tarfile.bz2");
    &maybe_run("bunzip2 -c $tarfile.bz2 | (cd $dir; tar xf -)");
    unless (defined $test) {
      die "Failed to create $mainFile from $tarfile.bz2 in the biolip step of download.pl\n"
        unless -s $mainFile;
      print STDERR "Successfully created $mainFile\n";
    }
  }
  # Fetch all the weekly updates and combine into BioLiP_UP_nr.txt
  &maybe_run("$Bin/download_BioLip_updates.pl");
  unless (defined $test) {
    die "Failed to assemble $biolipCacheDir/BioLiP_UP_nr.txt\n"
      unless -s "$biolipCacheDir/BioLiP_UP_nr.txt";
    print STDERR "Assembled $biolipCacheDir/BioLiP_UP_nr.txt\n";
  }
  # Save those in the data directory
  &maybe_run("cp $biolipCacheDir/BioLiP_UP_nr.txt $dir/");
}

unlink($listfile);
if (defined $test) {
    print STDERR "Finished test\n";
} else {
    print STDERR "Downloads complete\n";
}

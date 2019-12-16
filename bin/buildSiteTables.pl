#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw{$Bin};
use lib "$Bin/../lib";
use pbutils qw{ReadFastaEntry SQLiteLine};

my $outdir = ".";
my @outtypes = qw{binding modified functional mutagenesis};
my %outtypes = map { $_ => 1 } @outtypes;
my $usage = <<END
Usage: buildSiteTables.pl -sprot sprotFt.tab -biolip anno_nr.txt
       -seqres pdb_seqres.txt -sprotFasta uniprot_sprot.fasta.gz
  (-biolip may have multiple arguments)
Optional arguments:
  -outdir $outdir -- directory to put the output files in

Builds the Sites.tab, HasSites.tab, and hassites.faa files in the
output directory, and formats hassites.faa as a BLAST database

Sites.tab has the fields db, id, chain, ligandId, ligandChain, type,
posFrom, posTo, pdbFrom, pdbTo, comment, pmIds,
  where db is SwissProt or PDB and type is one of
  @outtypes

HasSites.tab has the fields db, id, chain (optional), id2 (optional), desc

hassites.faa has identifiers of the form db:id or db:id:chain
END
;

my $sprotFtFile;
my $sprotFastaFile;
my @biolipAnnoFiles;
my $pdbSeqFile;

die $usage
  unless GetOptions('sprot=s' => \$sprotFtFile,
                    'biolip=s{1,}' => \@biolipAnnoFiles,
                    'seqres=s' => \$pdbSeqFile,
                    'sprotFasta=s' => \$sprotFastaFile,
                    'outdir=s' => \$outdir)
  && @ARGV == 0;
die "Not a directory: $outdir\n" unless -d $outdir;
die "Must specify 1 or more biolip annotation files:\n$usage" unless @biolipAnnoFiles;
die "Must use -sprotFt:\n$usage" unless defined $sprotFtFile;
die "Must use -seqres:\n$usage" unless defined $pdbSeqFile;
die "Must use -sprotFasta:\n$usage" unless defined $sprotFastaFile;

my $formatdb = "$Bin/blast/formatdb";
die "No such executable: $formatdb\n" unless -x $formatdb;

foreach my $file ($pdbSeqFile,$sprotFastaFile,$sprotFtFile,@biolipAnnoFiles) {
  die "No such file: $file\n" unless -e $file;
}

# sprotId => [acc, desc]
my %sprotIds = ();

# Write out SwissProt and PDB sites
open(my $fhSites, ">", "$outdir/Sites.tab") || die "Cannot write to $outdir/Sites.tab";

open(my $fhFt, "<", $sprotFtFile) || die "Cannot read $sprotFtFile\n";
while(my $line = <$fhFt>) {
  chomp $line;
  my ($id, $acc, $desc, $type, $from, $to, $comment, $pmIds) = split /\t/, $line, -1;
  my ($from2, $comment2); # for the second line, for disulfides
  my $outtype;

  if ($type eq "ACT_SITE") {
    $outtype = "functional";
    $comment = "active site, $comment";
  } elsif ($type eq "BINDING") {
    $outtype = "binding";
  } elsif ($type eq "CA_BIND") {
    $outtype = "binding";
    $comment = "Ca2+"; # the comment field does not seem to be useful
  } elsif ($type eq "CARBOHYD") {
    $outtype = "modified";
    $comment = "carbohydrate, $comment";
  } elsif ($type eq "CHAIN") {
    $outtype = "modified";
    $comment = "mature protein, $comment"
  } elsif ($type eq "COILED") {
    next; # ignore
  } elsif ($type eq "CONFLICT") {
    $outtype = "modified";
    $comment = "sequence conflict, $comment";
  } elsif ($type eq "CROSSLNK") {
    $outtype = "modified";
    # two lines, one for each, unless from=to, which can happen if interchain
    if ($from ne $to) {
      my $commentOrig = $comment;
      $comment = "Crosslink with $to, $commentOrig";
      $from2 = $to;
      $comment2 = "Crosslink with $from, $commentOrig";
      $to = $from;
    }
  } elsif ($type eq "DISULFID") {
    $outtype = "modified";
    # two lines, one for each, unless from=to
    # If labeled as interchain, then usually from=to, unless it is
    # between two different copies of the same chain.
    # And there are also rare exceptions; so, ignore the interchain issue in the logic
    if ($from ne $to) {
      my $commentOrig = $comment;
      $comment = "Disulfide link with $to, $commentOrig";
      $from2 = $to;
      $comment2 = "Disulfide link with $from, $commentOrig";
      $to = $from;
    }
  } elsif ($type eq "DNA_BIND") {
    $outtype = "binding";
  } elsif ($type eq "DOMAIN") {
    next; # ignore, not specific to residues
  } elsif ($type eq "INIT_MET") {
    $outtype = "modified";
    $comment = "Initiator methionine, $comment";
  } elsif ($type eq "INTRAMEM") {
    next; # ignore
  } elsif ($type eq "LIPID") {
    $outtype = "modified";
  } elsif ($type eq "METAL") {
    $outtype = "binding";
  } elsif ($type eq "MOD_RES") {
    $outtype = "modified";
  } elsif ($type eq "MOTIF") {
    $outtype = "functional";
  } elsif ($type eq "MUTAGEN") {
    # Note this includes cases where mutating the residue has little effect,
    # so kept this separate from "functional"
    $outtype = "mutagenesis";
  } elsif ($type eq "NON_CONS") {
    $outtype = "modified";
    $comment = "Non-consecutive (sequence gap), $comment";
  } elsif ($type eq "NON_STD") {
    # worth highlighting, so, show here even though the Swiss-Prot sequence will have U or O
    $outtype = "modified";
    $comment = "nonstandard, $comment";
  } elsif ($type eq "NP_BIND") {
    $outtype = "binding";
  } elsif ($type eq "PEPTIDE") {
    $outtype = "modified";
    $comment = "released peptide, $comment";
  } elsif ($type eq "PROPEP") {
    $outtype = "modified";
    $comment = "propeptide, $comment";
  } elsif ($type eq "REGION") {
    $outtype = "functional";
  } elsif ($type eq "REPEAT") {
    next; # ignore
  } elsif ($type eq "SIGNAL") {
    $outtype = "modified";
    $comment = "signal peptide, $comment";
  } elsif ($type eq "SITE") {
    $outtype = "functional";
  } elsif ($type eq "TOPO_DOM") {
    next;
  } elsif ($type eq "TRANSIT") {
    $outtype = "modified";
    $comment = "transit peptide, $comment";
  } elsif ($type eq "TRANSMEM") {
    next;
  } elsif ($type eq "UNSURE") {
    $outtype = "modified";
    $comment = "sequence uncertain, $comment";
  } elsif ($type eq "VAR_SEQ") {
    $outtype = "modified";
    $comment = "Variant sequence, $comment";
  } elsif ($type eq "VARIANT") {
    # ignore if the comment contains only a gene name or strain name, i.e. "in OPD2" or "in strain:",
    # or dbSNP references, i.e. "dbSNP:rs121912635"
    $outtype = "modified";
    if ($comment =~ m/[(](.*)[)]/) {
      my $function = $1;
      my @pieces = split /; /, $function;
      @pieces = grep { !m/^in strain/ && !m/^dbSNP/ && !m/^in [A-Z0-9]+$/ } @pieces;
      next if @pieces == 0;
    }
    $comment = "natural variant: $comment";
  } elsif ($type eq "ZN_FING") {
    next; # ignore, like a domain annotation
  } else {
    print STDERR "Ignoring unknown feature type from Swiss-Prot: $type\n";
    next;
  }

  $comment =~ s/,? $//;
  $comment2 =~ s/,? $// if defined $comment2;
  die "Invalid out type $outtype from line\n$line"
    unless defined $outtype && exists $outtypes{$outtype};
  print $fhSites SQLiteLine("SwissProt", $id,
                            # no chain, ligandId, or ligandChain
                            "", "", "",
                            $outtype, $from, $to,
                            # no pdbfrom or pdbto
                            "", "",
                            $comment, $pmIds);
  if (defined $from2) {
    die "from2 set without comment2 from line\n$line" unless defined $comment2;
    print $fhSites SQLiteLine("SwissProt", $id,
                              "", "", "",
                              $outtype, $from2, $from2,
                              "", "",
                              $comment2, $pmIds);
  }
  # save this identifier
  $sprotIds{$id} = [$acc, $desc];
}
close($fhFt) || die "Error reading $sprotFtFile\n";

my %pdbSeq = (); # pdbId => chain => sequence
foreach my $biolipAnnoFile (@biolipAnnoFiles) {
  open(my $fhAnno, "<", $biolipAnnoFile) || die "Cannot read $biolipAnnoFile\n";
  my %seenActive = (); # pdbId => chain => position => 1
  while (my $line = <$fhAnno>) {
    chomp $line;
    my ($pdbId, $chain, $resolution, $siteId, $ligandId, $ligandChain, $ligandNo,
        $bindingResiduesPDB, $bindingResidues,
        $siteResiduesPDB, $siteResidues,
        $EC, $GO,
        $affinityPM, $affinityMOAD, $affinityPDBBind, $affinityBindingDB,
        $uniprotId, $pubmedIds, $seq) = split /\t/, $line;
    die "Invalid line $line" unless $seq =~ m/^[A-Z]+$/;
    if (exists $pdbSeq{$pdbId}{$chain}) {
      die "Inconsistent sequence for $pdbId:$chain"
        unless $pdbSeq{$pdbId}{$chain} = $seq;
    } else {
      $pdbSeq{$pdbId}{$chain} = $seq;
    }
    $pubmedIds =~ m/^[0-9,]*$/ || die "Invalid pubmed ids $pubmedIds in line\n$line";
    # "Different" sites are separated by ; but I do see residues repeat, so,
    # use a hash to ignore repeats
    my @bindingResiduesPDB = split /[; ]/, $bindingResiduesPDB;
    my @bindingResidues = split /[; ]/, $bindingResidues;
    my @siteResiduesPDB = split /[; ]/, $siteResiduesPDB;
    my @siteResidues = split /[; ]/, $siteResidues;
    foreach my $type ("functional","binding") {
      my $listPDB = $type eq "functional" ? \@siteResiduesPDB : \@bindingResiduesPDB;
      my $list = $type eq "functional" ? \@siteResidues : \@bindingResidues;
      my $lenAgree = scalar(@$list) == scalar(@$listPDB);
      print STDERR "Warning: Mismatched lengths of $type residues for $pdbId:$chain\n"
        unless $lenAgree;
      next if scalar(@$list) == 0;
      foreach my $i (0..(scalar(@$list)-1)) {
        my $pos = $list->[$i];
        $pos =~ s/^[A-Z]//;
        $pos =~ m/^\d+$/ || die "Invalid $type residue $pos in " . join(" ",@$list) . " from \n$line";
        my $posPDB = ""; # not reported if there are mismatches
        if ($lenAgree) {
          my $posPDB = $listPDB->[$i];
          $posPDB =~ s/^[A-Z]//;
          # Occasionally, PDB residue numbers have suffixes like A or D or a, I'm not sure why
          $posPDB =~ m/^-?\d+[A-Za-z]*$/ || die "Invalid PDB $type residue $posPDB in " . join(" ",@$listPDB) . " from \n$line";
        }
        if ($type eq "functional") {
          # prevent the same active site from being reported more than once (with different ligands)
          next if exists $seenActive{$pdbId}{$chain}{$pos};
          $seenActive{$pdbId}{$chain}{$pos} = 1;
        }
        print $fhSites join("\t", "PDB", $pdbId, $chain,
                            $type eq "functional" ? "" : $ligandId,
                            $type eq "functional" ? "" : $ligandChain,
                            $type,
                            $pos, $pos, $posPDB, $posPDB, "", $pubmedIds)."\n";
      }
    }
  }
  close($fhAnno) || die "Error reading $biolipAnnoFile\n";
}
print STDERR "Read sites for " . scalar(keys %sprotIds) . " SwissProt entries and " . scalar(keys %pdbSeq) . " PDB entries\n";
close($fhSites) || die "Error writing to $outdir/Sites.tab";
print STDERR "Wrote $outdir/Sites.tab\n";

# Write out hassites.faa by selecting sequences from $sprotFastaFile and writing out %pdbSeq
my $fhFaaIn;
if ($sprotFastaFile =~ m/[.]gz$/) {
  open($fhFaaIn, "-|", "zcat", $sprotFastaFile)
    || die "Cannot run zcat on $sprotFastaFile\n";
} else {
  open($fhFaaIn, "<", $sprotFastaFile) || die "Cannot read $sprotFastaFile\n";
}
open(my $fhFaaOut, ">", "$outdir/hassites.faa")
  || die "Cannot write to $outdir/hassites.faa\n";
my %sprotFound = (); # hash of ids whose sequence was found
my $state = {};
while (my ($header, $seq) = ReadFastaEntry($fhFaaIn, $state)) {
  $header =~ m/sp[|]([A-Z0-9]+)[|]/
    || die "Cannot extract SwissProt id from header line\n>$header\n";
  my $id = $1;
  if (exists $sprotIds{$id}) {
    die "Duplicate sequence for $id in $sprotFastaFile\n"
      if exists $sprotFound{$id};
    $sprotFound{$id} = 1;
    print $fhFaaOut ">SwissProt:$id\n$seq\n";
  }
}
close ($fhFaaIn) || die "Error reading $sprotFastaFile\n";
foreach my $pdbId (sort keys %pdbSeq) {
  my $hash = $pdbSeq{$pdbId};
  foreach my $chain (sort keys %$hash) {
    print $fhFaaOut ">PDB:$pdbId:$chain\n" . $hash->{$chain} . "\n";
  }
}

close($fhFaaOut) || die "Error writing to $outdir/hassites.faa\n";
print STDERR "Wrote $outdir/hassites.faa\n";
print STDERR "Warning: found sequences of " . scalar(keys %sprotFound)
  . "SwissProt identifiers, out of " . scalar(keys %sprotIds) . "\n"
  if scalar(keys %sprotFound) != scalar(keys %sprotIds);

system($formatdb,"-p","T","-o","T","-i","$outdir/hassites.faa") == 0
  || die "$formatdb failed -- $!";
print STDERR "Formatted $outdir/hassites.faa\n";

# Parse $pdbSeqFile
my %pdbDesc;
open(my $fhPdbSeq, "<", $pdbSeqFile) || die "Cannot read $pdbSeqFile\n";
$state = {};
while (my ($header,$seq) = ReadFastaEntry($fhPdbSeq, $state)) {
  $header =~ m/([0-9a-zA-Z]+)_([0-9a-zA-Z]+) mol:[a-z ]+:(\d+) +(.*)$/
    || die "Cannot parse header line from $pdbSeqFile\n$header";
  my ($pdbId, $chain, $seqlen, $desc) = ($1, $2, $3, $4);
  die "Incorrect sequence length for ${pdbId}_${chain}"
    unless $seqlen == length($seq);
  if (exists $pdbSeq{$pdbId}{$chain}) {
    # Do not check that sequences match -- having a few extra a.a at beginning or end is very common
    $pdbDesc{$pdbId}{$chain} = $desc;
  }
}
close($fhPdbSeq) || die "Error reading $pdbSeqFile\n";

# And write HasSites.tab
open(my $fhHasSites, ">", "$outdir/HasSites.tab") || die "Cannot write to $outdir/HasSites.tab";
foreach my $sprotId (sort keys %sprotIds) {
  my ($acc,$desc) = @{ $sprotIds{$sprotId} };
  print $fhHasSites SQLiteLine("SwissProt", $sprotId, "", $acc, $desc);
}

foreach my $pdbId (sort keys %pdbSeq) {
  my $hash = $pdbSeq{$pdbId};
  foreach my $chain (sort keys %$hash) {
    if (!exists $pdbDesc{$pdbId}{$chain}) {
      print STDERR "Warning, no description for ${pdbId}_${chain}\n";
    } else {
      print $fhHasSites SQLiteLine("PDB", $pdbId, $chain, "", $pdbDesc{$pdbId}{$chain});
    }
  }
}
close($fhHasSites) || die "Error writing to $outdir/HasSites.tab";
print STDERR "Wrote $outdir/HasSites.tab\n";

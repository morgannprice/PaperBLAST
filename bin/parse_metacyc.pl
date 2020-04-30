#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use pbutils; # for ParsePTools, ReadFastaEntry

my $usage = <<END
Usage: parse_metacyc.pl -data directory -faa uniprot_fasta trembl_fasta -out out

The data directory must contain proteins.dat, enzrxns.dat, and
reactions.dat (in attribute-value format).

The uniprot and trembl fasta filenames can be .gz files instead.

Writes to out.curated_parsed, out.reaction_links, and out.reaction_compounds

The curated_parsed file has the fields metacyc,
metacyc identifier, uniprot identifier, short name, description, organism,
sequence, blank (comment), pubmed ids (comma separated).

The reaction_links file has the fields reaction id (either RHEA:number
or metacyc:reaction-id), metacyc enzrxn id, and one or protein
identifiers (each of the form metacyc::protein-id)

The reaction_compounds file has the fields reaction id (either
RHEA:number or metacyc:reaction-id), reaction location (often empty),
and one or more participants, each of the form
  side:coefficient:compartment:compound-id:compound-description
where side is -1 for left side and +1 for right side. The leading CCO-
has been removed from the compartments.

Limitations -- EC numbers are obtained for proteins that link to
enzrxns.dat (via CATALYZES) and reactions.dat (via REACTION), or from
proteins that belong to complexes that have CATALYZES attributes. But
multi-level enzyme complexes (A is a COMPONENT-OF B which is a
COMPONENT-OF C which has CATALYZES) are not handled.

Optional arguments:
-debug -- verbose output and allow running with no fasta input (for testing)
END
;

my @faaFiles;
my ($metadir, $outpre, $debug);
die $usage
  unless GetOptions('data=s' => \$metadir,
                    'faa=s{,}' => \@faaFiles,
                    'out=s' => \$outpre,
                    'debug' => \$debug)
  && defined $metadir && defined $outpre;
die "No such directory: $metadir\n" unless -d $metadir;
die "Must specify faa files unless using -debug\n"
  if @faaFiles == 0 && !defined $debug;
foreach my $faa (@faaFiles) {
  die "No such file: $faa\n" unless -e $faa;
}

# A list of hashes, one per protein. Each entry must include UNIQUE-ID and may include:
# ABBREV-NAME, COMMON-NAME, CATALYZES (as a list of enzrxnIds),
# UNIPROT, CITATIONS (as a list of pmIds),
# COMPONENT-OF, COMPONENTS (as a list of identifiers), GENE
#
# Unlike in earlier versions of this code, @prot includes every protein
# object (even if it has no pmIds or no link to sequence)

my @prot = ();

my $pfile = "$metadir/proteins.dat";
open(my $fhProt, "<", $pfile) || die "Cannot read $pfile\n";
while (my $prot = ParsePTools($fhProt)) {
  # Note that these objects contain SPECIES attributes of the form
  # TAX-nnnn, where the number is an NCBI taxonomy id, but that is
  # ignored in favor of the OS field in the uniprot fasta file.

  my $id = $prot->{"UNIQUE-ID"}[0]{"value"};
  die unless $id;
  my $obj = { 'UNIQUE-ID' => $id };

  my @catalyzes = ();
  if (exists $prot->{CATALYZES}) {
    @catalyzes = map $_->{value}, @{ $prot->{CATALYZES} };
    $obj->{CATALYZES} = \@catalyzes;
  }

  # Try to link to sequence
  my @dblinks = map { $_->{"value"} } @{ $prot->{"DBLINKS"} };
  my @uniprot = grep m/UNIPROT/, @dblinks;
  if (@uniprot > 0) {
    my $uniprotLine = $uniprot[0];
    if ($uniprotLine =~ m/[(]UNIPROT "([0-9A-Z_]+)-?\d*"/) {
      $obj->{UNIPROT} = $1;
    } else {
      print STDERR "Cannot parse uniprot id from DBLINKS entry: $uniprotLine\n";
    }
  }

  # Links to papers may be in the CITATIONS attribute or in the CITATIONS annotation of GO-TERMS
  # or in the COMMENT field as entries like |CITS: [8621577]| or |CITS: [2693216][1425658]| or |CITS: [9755155] [9388228]|
  my @pmLines = ();
  foreach my $value ( @{ $prot->{"CITATIONS"} }) {
    push @pmLines, $value->{"value"};
  }
  foreach my $value ( @{ $prot->{"GO_TERMS"} }) {
    push @pmLines, $value->{CITATIONS}
      if exists  $value->{CITATIONS};
  }
  my %pmIds = ();
  foreach my $line (@pmLines) {
    my $pmId = $1
      if $line =~ m/^(\d+):/ || $line =~ m/^(\d+)$/ || $line =~ m/^\[(\d+)\]$/;
    # Some citations are like Wang02, not sure why
    $pmIds{$pmId} = 1 if $pmId;
  }
  my $comment = $prot->{"COMMENT"}[0]{"value"} || "";
  while ($comment =~ m/[|]CITS: ([\[\]0-9 ]+)[|]/g) {
    # this should extract a field like [8801422] or [10501935][11872485] or [6706930] [7142155]
    my $ids = $1;
    $ids =~ s/[\[\]]/ /g;
    my @ids = split / +/, $ids;
    foreach my $id (@ids) {
      $pmIds{$id} = 1 if $id;
    }
  }
  my @pmIds = sort { $a <=> $b } keys %pmIds;
  $obj->{CITATIONS} = \@pmIds;

  foreach my $key (qw/ABBREV-NAME COMMON-NAME COMPONENT-OF GENE/) {
    $obj->{$key} = $prot->{$key}[0]{"value"} || "";
  }
  my @components = ();
  foreach my $value ( @{ $prot->{"COMPONENTS"} }) {
    push @components, $value->{"value"};
  }
  $obj->{COMPONENTS} = \@components if @components > 0;
  push @prot, $obj;
}
close($fhProt) || die "Error reading $pfile\n";
my %prot = map { $_->{"UNIQUE-ID"} => $_ } @prot;

my $nWithPM = scalar(grep @{ $_->{CITATIONS} } > 0, @prot);
my $nWithUniProt = scalar(grep exists $_->{UNIPROT}, @prot);
print STDERR "Read " . scalar(@prot) . " protein entries, $nWithPM with PubMid ids, $nWithUniProt with UniProt ids\n";

my %enzrxn = (); # unique id => common name
my %enzrxnToRxn = (); # unique id => unique id

my %compounds = (); # compoundId => hash of compoundName, keggLigand, formula (as string)
my $compfile = "$metadir/compounds.dat";
open(my $fhc, "<", $compfile) || die "Cannot read $compfile\n";
while (my $cmp = ParsePTools($fhc)) {
  my $cmpId = $cmp->{"UNIQUE-ID"}[0]{"value"};
  die "Invalid compound in $compfile" unless $cmpId;
  my $name = $cmp->{"COMMON-NAME"}[0]{"value"} || "";
  my @formula = ();
  foreach my $l (@{ $cmp->{"CHEMICAL-FORMULA"} }) {
    my $string = $l->{"value"};
    die "Invalid component of chemical formula $string" unless $string =~ m/^[(]([A-Za-z-]+) (\d+)[)]$/;
    my ($element, $cnt) = ($1,$2);
    push @formula, ucfirst(lc($element)) . $cnt;
  }
  my $keggLigand = "";
  foreach my $l (@{ $cmp->{"DBLINKS"} }) {
    my $string = $l->{"value"};
    $keggLigand = $1 if $string =~ m/^[(]LIGAND-CPD "([A-Z0-9]+)"/;
  }
  die "Duplicate compound $cmpId" if exists $compounds{$cmpId};
  $compounds{$cmpId} = { "compoundName" => $name,
                         "keggLigand" => $keggLigand,
                         "formula" => join("", @formula) };
}
close($fhc) || die "Error reading $compfile";
print STDERR "Read " . scalar(keys %compounds) . " compounds\n";

# Classes are hierarchical, so a compound class is not directly labelled as such --
# one needs to infer that it is a compound by going up the hierarchy
# (i.e. All-Carbohydrates has TYPES = Compounds and Glycans has TYPES = Carbohydrates.)
# Instead, will just assume that it is a compound if a pathway or reaction refers to it
my %classes = (); # id => hash that includes id, name and types (also a hash)
my $classfile = "$metadir/classes.dat";
open(my $fhcc, "<", $classfile)
  || die "Cannot open $classfile";
while(my $cc = ParsePTools($fhcc)) {
  my $id = $cc->{"UNIQUE-ID"}[0]{"value"};
  die unless $id;
  die "Duplicate class $id" if exists $classes{$id};
  my %obj = ( "id" => $id, "types" => {}, "name" => "" );
  foreach my $l (@{ $cc->{"TYPES"} }) {
    $obj{"types"}{ $l->{"value"} } = 1;
  }
  $obj{"name"} = $cc->{"COMMON-NAME"}[0]{"value"}
    if exists $cc->{"COMMON-NAME"};
  $classes{$id} = \%obj;
}
close($fhcc) || die "Error reading $classfile";

# Parse the enzymatic reactions file to get descriptions for the catalyzed reactions
my $enzfile = "$metadir/enzrxns.dat";
open(my $fhEnz, "<", $enzfile) || die "Cannot read $enzfile";
while (my $enzrxn = ParsePTools($fhEnz)) {
  my $id = $enzrxn->{"UNIQUE-ID"}[0]{"value"};
  my $desc = $enzrxn->{"COMMON-NAME"}[0]{"value"};
  $enzrxn{$id} = $desc if $id && $desc;
  $enzrxnToRxn{$id} = $enzrxn->{"REACTION"}[0]{"value"}
    if exists $enzrxn->{"REACTION"};
}
close($fhEnz) || die "Error reading $fhEnz";
print STDERR "Parsed descriptions for " . scalar(keys %enzrxn) . " enzymatic reactions\n";

my %rxnToEc = ();
my %rxnToRhea = (); # rxnId => list of rhea id
 # rxnId => list of hashes with "compoundId", "coefficient" "side" (+1 for right side; -1 for left) and optionally "compartment"
my %rxnCompounds = ();
my %rxnLocation = ();
my $reactfile = "$metadir/reactions.dat";
open (my $fhReact, "<", $reactfile) || die "Cannot read $reactfile";
while (my $react = ParsePTools($fhReact)) {
  my $id = $react->{"UNIQUE-ID"}[0]{"value"};
  my $EC = $react->{"EC-NUMBER"}[0]{"value"} || "";
  # Many reactions have values like "|EC-1.17.98.c|" which seem not to be official EC numbers
  # and so are not saved.
  # (Although RXN-16226 is given |EC-1.8.3.b| which is supposedly official.)
  # Some other reactions have partially specified EC numbers -- again, skip those.
  # Minor bug: numbers like |EC-4.2.2.n1| are actually official and should be saved
  if ($EC =~ m/^EC-([0-9]+[.][0-9]+[.][0-9]+[.][0-9]+)$/) {
    die unless $id;
    $rxnToEc{$id} = $1;
  } else {
    # print STDERR "Ignoring non-standard EC number $EC for $id\n" if $EC ne "";
  }
  my @dblinks = map { $_->{"value"} } @{ $react->{"DBLINKS"} };
  foreach my $dblink (@dblinks) {
    if ($dblink =~ m/^[(]RHEA "(\d+)" /) {
      my $rheaId = "RHEA:" . $1;
      push @{ $rxnToRhea{$id} }, $rheaId;
    }
  }
  my @cmp = ();
  foreach my $l (@{ $react->{"LEFT"} }) {
    $l->{side} = -1;
    push @cmp, $l;
  }
  foreach my $l (@{ $react->{"RIGHT"} }) {
    $l->{side} = +1;
    push @cmp, $l
  }
  # set up coefficients and compartments
  my @cmp2 = ();
  foreach my $l (@cmp) {
    my $cmpId = $l->{"value"};
    my $coefficient = $l->{"COEFFICIENT"} || 1;
    my $compartment = $l->{"COMPARTMENT"} || "";
    $compartment =~ s/^CCO-//;
    $compartment =~ s/://g; # not sure this happens but it would break the output format
    push @cmp2, { "compoundId" => $cmpId,
                  "side" => $l->{"side"},
                  "coefficient" => $coefficient,
                  "compartment" => $compartment };
  }
  # ParseMetaCycPathways.pl from the FEBA code base tries to identify duplicate compounds,
  # and sums the coefficients. Do not do that here.
  $rxnCompounds{$id} = \@cmp2;
  $rxnLocation{$id} = $react->{"RXN-LOCATIONS"}[0]{value} || "";
}
close($fhReact) || die "Error reading $reactfile";
print STDERR "Read " . scalar(keys %rxnCompounds) . " reactions from $reactfile\n";

# Link proteins to EC numbers via enzrxns and via COMPONENT-OF
foreach my $prot (@prot) {
  # First, try to link it to enzymatic reactions, either directly or via component-of
  my $enzrxns = [];
  if (exists $prot->{CATALYZES}) {
    $enzrxns = $prot->{CATALYZES};
  } elsif (exists $prot->{"COMPONENT-OF"}) {
    my $id2 = $prot->{"COMPONENT-OF"};
    if (exists $prot{$id2} && exists $prot{$id2}{CATALYZES}) {
      $enzrxns = $prot{$id2}{CATALYZES};
    }
  }

  foreach my $enzrxnId (@$enzrxns) {
    $prot->{"COMMON-NAME"} = $enzrxn{$enzrxnId}
      if $prot->{"COMMON-NAME"} eq "";
    if (exists $enzrxnToRxn{$enzrxnId}) {
      my $rxnId = $enzrxnToRxn{$enzrxnId};
      if (exists $rxnToEc{$rxnId}) {
        push @{ $prot->{"EC"} }, $rxnToEc{$rxnId};
      }
    }
  }
}

my @protNamed = grep { $_->{"COMMON-NAME"} && exists $_->{UNIPROT} } @prot;
print STDERR "Reduced to " . scalar(@protNamed) . " proteins with sequences and names\n";

# Which sequences we are looking for; empty string for as-yet unknown sequence
my %uniprot = ();
foreach my $prot (@prot) {
  if (exists $prot->{UNIPROT}) {
    $uniprot{ $prot->{UNIPROT} } = "";
  }
}
my %uniprotOrg = ();

# And now fetch sequences from the fasta files
foreach my $faafile (@faaFiles) {
  my $fh;
  if ($faafile =~ m/[.]gz$/) {
    open($fh, "zcat $faafile |") || die "Cannot run zcat on $faafile";
  } else {
    open($fh, "<", $faafile);
  }
  my $state = {};
  while(my ($header,$seq) = ReadFastaEntry($fh,$state)) {
    # ignore non-full-length sequences because they are difficult to interpret for annotation
    next if $header =~ m/[(]Fragment[)]/;
    my ($ids) = split / /, $header;
    my @ids = split /[|]/, $ids;
    shift @ids; # ignore sp or tr at front
    die "Invalid header $header" if @ids == 0;
    foreach my $id (@ids) {
      if (exists $uniprot{$id}) {
        $uniprot{$id} = $seq;
        if ($header =~ m/ OS=([^=]+) /) {
          $uniprotOrg{$id} = $1;
        }
      }
    }
  }
  close($fh) || die "Error reading $faafile";
}

# Fetch gene names from genes.dat
my %geneName = ();
my $genesfile = "$metadir/genes.dat";
open(my $fhGenes, "<", $genesfile) || die "Cannot read $genesfile";
while(my $gene = ParsePTools($fhGenes)) {
  my $geneId = $gene->{"UNIQUE-ID"}[0]{value};
  die unless $geneId;
  my $geneName = $gene->{"COMMON-NAME"}[0]{value};
  $geneName{$geneId} = $geneName if $geneName;
}
close($fhGenes) || die "Error reading $genesfile";
print STDERR "Read names for " . scalar(keys %geneName) . " genes\n";

open(my $fhCP, ">", "$outpre.curated_parsed")
  || die "Cannot write to $outpre.curated_parsed\n";
my %protWritten = (); # UNIQUE-ID => 1 if monomer saved to $outpre.curated_parsed
foreach my $prot (@protNamed) {
  my $uniprotId = $prot->{"UNIPROT"};
  next unless exists $uniprot{$uniprotId};
  my $seq = $uniprot{$uniprotId};
  my $org = $uniprotOrg{$uniprotId} || "";
  if ($seq eq "" && !defined $debug) {
    print STDERR "Warning: no sequence found for $uniprotId\n";
  } else {
    # Does this entry link to a paper, either directly or via a complex?
    my $cit = $prot->{CITATIONS};
    if (@$cit == 0) {
      my $id2 = $prot->{"COMPONENT-OF"};
      $cit = $prot{$id2}{CITATIONS} if exists $prot{$id2};
    }
    if (@$cit > 0) { # only entries linked to papers
      my $desc = $prot->{"COMMON-NAME"};
      if (exists $prot->{"EC"}) {
        my @ec = @{ $prot->{"EC"} };
        my %seen = ();
        my @parts = ();
        foreach my $ec (@ec) {
          next if exists $seen{$ec};
          $seen{$ec} = 1;
          push @parts, "EC $ec";
        }
        $desc .= " (" . join("; ", @parts) . ")";
      }
      $protWritten{ $prot->{"UNIQUE-ID"} } = 1;
      my $id2 = $prot->{"ABBREV-NAME"} || "";
      $id2 = $geneName{$prot->{GENE}} if
        $id2 eq "" && exists $geneName{$prot->{GENE}};
      print $fhCP join("\t",
                       "metacyc", $prot->{"UNIQUE-ID"}, $uniprotId,
                       $id2, $desc, $org,
                       $seq, "", join(",", @$cit)
                      )."\n";
    }
  }
}
close($fhCP) || die "Error writing to $outpre.curated_parsed";
print STDERR "Wrote " . scalar(keys %protWritten) . " entries to $outpre.curated_parsed\n";

# And save reaction links
open (my $fhRxn, ">", "$outpre.reaction_links") || die "Cannot write to $outpre.reaction_links";
foreach my $prot (@prot) {
  my $enzrxnIds = $prot->{CATALYZES};
  next unless defined $enzrxnIds && @$enzrxnIds > 0;

  # Check that we have all of the relevant sequences
  my @parts = ();
  if (exists $prot->{COMPONENTS}) {
    @parts = map $prot{$_}, @{ $prot->{COMPONENTS} };
  } else {
    @parts = ($prot);
  }
  my @partsDef = grep defined $_ && exists $_->{"UNIQUE-ID"}, @parts;
  next unless scalar(@partsDef) == scalar(@parts);

  my @partsWritten = grep exists $protWritten{$_->{"UNIQUE-ID"}}, @parts;
  next unless scalar(@parts) == scalar(@partsWritten);

  my %rxnSeen = ();
  foreach my $enzrxnId (@$enzrxnIds) {
    my $rxnId = $enzrxnToRxn{$enzrxnId};
    next if exists $rxnSeen{$rxnId};
    $rxnSeen{$rxnId} = 1;
    my @ids = map "metacyc::" . $_->{"UNIQUE-ID"}, @parts;
    print $fhRxn join("\t", "metacyc:".$rxnId, $enzrxn{$enzrxnId} || "", @ids)."\n";
    if (exists $rxnToRhea{$rxnId}) {
      foreach my $rhea (@{ $rxnToRhea{$rxnId}}) {
        print $fhRxn join("\t", $rhea, $enzrxn{$enzrxnId} || "", @ids)."\n";
      }
    }
  }
}
close($fhRxn) || die "Error writing to $outpre.reaction_links";

open (my $fhrc, ">", "$outpre.reaction_compounds")
  || die "Cannot write to $outpre.reaction_compounds\n";
foreach my $rxnId (sort keys %rxnCompounds) {
  my $cmps = $rxnCompounds{$rxnId};
  next if @$cmps == 0; # not sure this actually happens
  my $showId = "metacyc::" . $rxnId;
  $showId = $rxnToRhea{$rxnId}[0] if exists $rxnToRhea{$rxnId};
  my @line = ( $showId, $rxnLocation{$rxnId} );
  foreach my $cmp (@$cmps) {
    my $id = $cmp->{compoundId};
    my $cmpName = $id;
    if (exists $compounds{$id}) {
      $cmpName = $compounds{$id}{compoundName};
    } elsif (exists $classes{$id}) {
      $cmpName = $classes{$id}{name} if exists $classes{$id}{name};
      # just use the class name as the id if no COMMON-NAME was provided
    } else {
      print STDERR "Warning: unknown compound $id in $rxnId\n" if $cmpName eq "";
    }
    push @line, join(":", $cmp->{side}, $cmp->{coefficient}, $cmp->{compartment}, $cmp->{compoundId},
                    $cmpName);
  }
  print $fhrc join("\t", @line)."\n";
}
close($fhrc) || die "Error writing to $outpre.reaction_compounds\n";

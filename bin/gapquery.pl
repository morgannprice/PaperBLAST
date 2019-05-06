#!/usr/bin/perl -w
# Find the relevant queries for a set of gaps
# Operates on the PaperBLAST database

use strict;
use Getopt::Long;
use DBI;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use Steps qw{ReadSteps FetchUniProtSequence};
use pbutils qw{CuratedMatch CuratedWordMatch ReadTable IdToUniqId UniqIdToSeq FetchCuratedInfo};

sub MyCuratedMatch($$); # paperblast database and term to word matches
# (returned as a reference to a list of rows in CuratedGene)

my $debug;

{
  my $pbdir = "papers.data";
  my $bindir = ".";
  my $tigrinfoFile = "tigrinfo";
  my $pfaminfoFile = "pfam.tab";
  my $tigrfile = "TIGRFAMs.hmm";
  my $pfamfile = "Pfam-A.hmm";
  my $dir;
  my $usage = <<END
Usage: gapquery.pl -steps stepsfile -dir .

stepsfile includes the following lines:
comment lines (that start with #)

step lines, which are tab-delimited with fields step name, step description,
  and attributes of the form type:value. Examples of attributes:
    EC:1.1.1.1	                   term:search string
    hmm:PF00001	hmm:TIGR0001
    uniprot:Q8TPT3_METAC	   curated:reanno::Miya:8499492
    ignore:CharProtDB::CH_123581   ignore_other:EC 2.2.1.6

  The step name and step description are arbitrary.
  EC numbers are expanded to matching TIGRFams as well.

  To obtain the list of queries for finding candidates by BLAST, all
  curated items matching the attributes are combined, and then those
  in the ignore list are removed. (The ignore list also affects
  scoring of reverse-hits later on.)

rules lines starting with name: followed by a list of step names or
rule names.

blank (whitespace only lines)

This script considers only the step lines.

Writes to dir/steps.query, a tab-delimited file with the fields
step, type, query, desc, file (optional), sequence (optional)
type is one of curated, hmm, uniprot, or ignore
query is the sequence identifier(s), HMM identifier, or the uniprot identifier
  (If multiple curated items have the same sequence, all identifiers are joined by ",")
sequence is not set for hmm items; those have file set instead (i.e., PF02965.17.hmm)
Also creates the HMM files in outdir.

Optional arguments:
-pbdir $pbdir -- data directory including litsearch.db and uniq.faa
-bin $bindir -- directory with executables for blast, usearch, hmmer
-tigrinfo $tigrinfoFile -- tab-delimited file describing TIGRFams
-pfaminfo $pfaminfoFile -- tab-delimited file describing PFams
-pfamfile $pfamfile -- the indexed HMM database for PFam
-tigrfile $tigrfile -- the indexed HMM database for TIGRFam
-debug
END
;

  my $stepsfile;
  my $outdir;
  die $usage
    unless GetOptions('pbdir=s' => \$pbdir,
                      'bin=s' => \$bindir,
                      'tigrinfo=s' => \$tigrinfoFile,
                      'pfaminfo=s' => \$pfaminfoFile,
                      'tigrfile=s' => \$tigrfile,
                      'pfamfile=s' => \$pfamfile,
                      'steps=s' => \$stepsfile,
                      'dir=s' => \$outdir,
                      'debug' => \$debug)
      && defined $stepsfile && defined $outdir;
  foreach my $dir ($pbdir, $bindir, $outdir) {
    die "No such directory: $dir\n" unless -d $dir;
  }
  foreach my $file ($stepsfile, $tigrinfoFile, $pfaminfoFile, $tigrfile, $pfamfile) {
    die "No such file: $file\n" unless -e $file;
  }
  my $hmmfetch = "$bindir/hmmer/hmmfetch";
  foreach my $x ($hmmfetch) {
    die "No such executable: $x\n" unless -x $x;
  }

  my @tigrinfo = ReadTable($tigrinfoFile, ["tigrId","ec","definition","type"]);
  my %tigrinfo = map { $_->{tigrId} => $_ } @tigrinfo;
  my %ecTIGR = ();
  foreach my $row (@tigrinfo) {
    my $ec = $row->{ec};
    push @{ $ecTIGR{$ec} }, $row if $ec ne "";
  }

  my @pfam = ReadTable($pfaminfoFile, ["acc","name"]);
  my %pfam = map { $_->{acc} => $_ } @pfam;
  my %pfToAcc = ();
  foreach my $acc (keys %pfam) {
    $acc =~ m/^(PF\d+)[.]\d+$/ || die "Invalid pfam accession $acc in $pfaminfoFile";
    my $pf = $1;
    $pfToAcc{$pf} = $acc;
  }

  my $pb = DBI->connect("dbi:SQLite:dbname=$pbdir/litsearch.db", "", "", { RaiseError => 1 }) || die $DBI::errstr;

  my $st = ReadSteps($stepsfile);
  my $steps = $st->{steps};
  my $rules = $st->{rules};

  # Check that each step or rule is used as a dependency, except for the rule named all
  my %dep = ();
  foreach my $lists (values %$rules) {
    foreach my $list (@$lists) {
      foreach my $member (@$list) {
        $dep{$member} = 1;
      }
    }
  }
  foreach my $step (keys %$steps) {
    print STDERR "Warning: step $step is not used by any rule\n"
      unless exists $dep{$step};
  }
  foreach my $rule (keys %$rules) {
    print STDERR "Warning: rule $rule is not used by any other rule\n"
      unless exists $dep{$rule} || $rule eq "all";
  }

  my %hmm = (); # hmm Id => relevant step => 1
  my %curated = (); # curated db => curated id => row in CuratedGene
  my %curatedToStep = (); # db => id => step => 1
  my %uniprotToStep = (); # uniprotId => step => 1
  my %uniprotSeq = ();
  my %uniprotDesc = ();
  my %ignore = (); # db => protId => step => desc
  foreach my $step (sort keys %$steps) {
    my $desc = $steps->{$step}{desc};
    my $l = $steps->{$step}{search};
    my @curatedThis = ();
    my @hmmThis = ();
    foreach my $row (@$l) {
      my ($type, $value) = @$row;
      my $curated = [];
      my @hmm = ();
      if ($type eq "EC") {
        $curated = MyCuratedMatch($pb, $value);
        push @hmm, map { $_->{tigrId} } @{ $ecTIGR{$value} } 
          if exists $ecTIGR{$value};
        print STDERR "$step\t$desc: no curated hits or TIGRFam for EC:$value\n"
          if @$curated == 0 && @hmm == 0;
      } elsif ($type eq "hmm") {
        push @hmm, $value;
      } elsif ($type eq "term") {
        $curated = MyCuratedMatch($pb, $value);
        print STDERR "$step\t$desc: no curated hits for $value\n"
          unless @$curated;
      } elsif ($type eq "curated") {
        $value =~ m/^([^:]+)::(.*)$/
          || die "Invalid curated protein specifier $value\n";
        my ($db, $protId) = ($1,$2);
        my $row = $pb->selectrow_hashref("SELECT * FROM CuratedGene WHERE db = ? AND protId = ?",
                                          {}, $db, $protId);
        die "Unknown curated protein $db $protId\n"
          unless defined $row && exists $row->{desc};
        push @$curated, $row;
      } elsif ($type eq "uniprot") {
        # Verify that this is a valid uniprot id
        die "Invalid uniprot identifier $value -- wrong format\n"
          unless $value =~ m/^[0-9A-Z_]+$/;
        unless (exists $uniprotSeq{$value}) {
          my ($seq, $desc) = FetchUniProtSequence($value);
          die "Cannot fetch uniprot sequence for identifier $value, is it invalid?\n" unless $seq;
          $uniprotSeq{$value} = $seq;
          $uniprotDesc{$value} = $desc;
        }
        $uniprotToStep{$value}{$step} = 1;
      } elsif ($type eq "ignore") {
        # Verify that this is a valid identifier
        die "Invalid ignore identifier $value -- should be of the form db::protId\n"
          unless $value =~ m/^[0-9A-Za-z_.]+::[0-9A-Za-z_.:-]+$/;
        my ($db,$protId) = split /::/, $value;
        my ($desc) = $pb->selectrow_array("SELECT desc FROM CuratedGene WHERE db = ? AND protId = ?",
                                          {}, $db, $protId);
        if (defined $desc) {
          $ignore{$db}{$protId}{$step} = $desc;
        } else {
          print STDERR "Warning: no curated rows for ignored identifier $value\n";
        }
      } elsif ($type eq "ignore_other") {
        # ignore curated items that match the term unless they have already matched
        my $list = MyCuratedMatch($pb, $value);
        print STDERR "$step\t$value: no curated hits for $value\n"
          unless @$list;
        foreach my $row (@$list) {
          $ignore{ $row->{db} }{ $row->{protId} }{$step} = $row->{desc}
            unless exists  $curatedToStep{ $row->{db} }{ $row->{protId} }{$step};
        }
      } else {
        die "Unknown attribute type $type in step $step\n";
      }
      if (@$curated) {
        push @curatedThis, @$curated;
        foreach my $row (@$curated) {
          $curated{ $row->{db} }{ $row->{protId} } = $row;
          $curatedToStep{ $row->{db} }{ $row->{protId} }{$step} = 1;
        }
      }
      foreach my $hmm (@hmm) {
        $hmm{$hmm}{$step} = 1;
      }
      push @hmmThis, @hmm;
    } # end loop over search terms
    my @curated2 = grep { !exists $ignore{$_->{db}}{$_->{protId}}{$step} } @curatedThis;
    print STDERR "Warning: no hmm or curated items for $step\n"
      if @curated2 == 0 && @hmmThis == 0;
    print STDERR join("\t", $step, "nCurated", scalar(@curatedThis),
                      "nNotIgnore", scalar(@curated2),
                      "HMM", scalar(@hmmThis))."\n" if defined $debug;
  }

  my $outfile = $stepsfile;
  $outfile =~ s!^.*[/]!!;
  $outfile =~ s/[.].*$//;
  $outfile = "$outdir/$outfile.query";
  open(my $fhO, ">", $outfile) || die "Cannot write to $outfile\n";

  print $fhO join("\t", qw{step type query desc file sequence})."\n";
  foreach my $db (sort keys %curated) {
    my $hash = $curated{$db};
    foreach my $protId (sort keys %$hash) {
      foreach my $step (sort keys %{ $curatedToStep{$db}{$protId} }) {
        if (exists $ignore{$db}{$protId}{$step}) {
          print STDERR "Ignoring entry $db $protId for $step (it is explicitly ignored, and also matches)\n";
        } else {
          my $uniqId = IdToUniqId($pb, $db."::".$protId);
          die "$db $protId" unless $uniqId;
          my $curated = FetchCuratedInfo($pb, $uniqId);
          die "$db $protId" unless @$curated > 0;
          my @ids = map { $_->[0] . "::" . $_->[1] } @$curated;
          my @descs = map { $_->[2] } @$curated;
          my $seq = UniqIdToSeq($pbdir, "$Bin/blast", $uniqId);
          print $fhO join("\t", $step, "curated",
                          join(",", @ids), join(";; ", @descs), "", $seq)."\n";
        }
      }
    }
  }
  foreach my $hmmId (sort keys %hmm) {
    my ($hmmdb, $acc, $desc);
    if ($hmmId =~ m/^PF\d+$/) {
      die "Unknown PFam $hmmId\n" unless exists $pfToAcc{$hmmId};
      $acc = $pfToAcc{$hmmId};
      $desc = $pfam{$pfToAcc{$hmmId}}{'name'};
      $hmmdb = $pfamfile;
    } elsif ($hmmId =~ m/^TIGR\d+$/) {
      die "Unknown HMM $hmmId\n"
        unless exists $tigrinfo{$hmmId};
      my $ti = $tigrinfo{$hmmId};
      print STDERR "Warning: $hmmId for steps " .  join(" ", sort keys %{ $hmm{$hmmId} })
        . " is of type $ti->{type}\n"
        unless $ti->{type} eq "equivalog" || $ti->{type} eq "equivalog_domain";
      $acc = $hmmId;
      $hmmdb = $tigrfile;
      $desc = $ti->{definition};
      $desc = $ti->{geneSymbol} . ": " . $desc if $ti->{geneSymbol};
      $desc .= " (EC $ti->{ec})" if $ti->{ec};
    } else {
      die "Unknown HMM $hmmId for steps " . join(" ", sort keys %{ $hmm{$hmmId} });
    }

    my $hmmfile = "$outdir/$acc.hmm";
    if (-e $hmmfile) {
      print STDERR "Using existing file $hmmfile\n";
    } else {
      my $cmd = "$hmmfetch $hmmdb $acc > $hmmfile.tmp";
      system($cmd)==0 || die "Failed to fetch hmm $acc from $hmmdb: $!";
      rename("$hmmfile.tmp", $hmmfile) || die "Failed to rename to $hmmfile";
      print STDERR "Created $hmmfile\n";
    }
    foreach my $step (sort keys %{ $hmm{$hmmId} }) {
      print $fhO join("\t", $step, "hmm", $hmmId, $desc, "$acc.hmm", "")."\n";
    }
  }
  foreach my $uniprotId (sort keys %uniprotToStep) {
    my $seq = $uniprotSeq{$uniprotId} || die;
    foreach my $step (sort keys %{ $uniprotToStep{$uniprotId} }) {
      print $fhO join("\t", $step, "uniprot", $uniprotId, $uniprotDesc{$uniprotId}, "", $seq)."\n";
    }
  }
  foreach my $db (sort keys %ignore) {
    foreach my $protId (sort keys %{ $ignore{$db} }) {
      my $hash = $ignore{$db}{$protId};
      my $uniqId = IdToUniqId($pb, $db."::".$protId);
      die "$db $protId" unless $uniqId;
      my $curated = FetchCuratedInfo($pb, $uniqId);
      die "$db $protId" unless @$curated > 0;
      my @ids = map { $_->[0] . "::" . $_->[1] } @$curated;
      my @descs = map { $_->[2] } @$curated;
      my $seq = UniqIdToSeq($pbdir, "$Bin/blast", $uniqId);

      foreach my $step (sort keys %$hash) {
        print $fhO join("\t", $step, "ignore",
                        join(",", @ids), join(";; ", @descs), "", $seq)."\n";
      }
    }
  }
  close($fhO) || die "Error writing to $outfile\n";
  print STDERR "Wrote $outfile\n";
}

sub MyCuratedMatch($$) {
  my ($dbh, $term) = @_;
  die "Searching for empty term"
    unless defined $term && $term ne "";
  my $limit = 10000;
  my $chits = CuratedMatch($dbh, $term, $limit);
  die "Limit reached for curated descriptions matching $term\n"
    if @$chits >= $limit;
  return CuratedWordMatch($chits, $term); # force word matches
}

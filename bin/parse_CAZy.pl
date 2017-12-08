#!/usr/bin/perl -w
use strict;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use pbutils;
use List::MoreUtils qw(uniq);

my $usage = <<END
parse_CAZy.pl faafile ecinfo

  Given the subset of entries that have EC numbers in CAZy, i.e. as downloaded from
  http://csbl.bmb.uga.edu/dbCAN/download/CAZyDB.07202017.fa
  produces a curated_parsed file, with the fields
  database id, protein identifier, blank (secondary identifier), blank (short name), description, blank (organism),
  sequence, blank (comment), blank (pmids)

  Also requires metadata on the EC numbers, as in
  http://csbl.bmb.uga.edu/dbCAN/download/CAZyDB-ec-info.txt.07-20-2017
  with the fields sequence identifier, family(ignored), ec number, and name.
  This additional information is necessary because many of the EC numbers are non-specific,
  such as 1.-.-.- for lytic monooxygenases with various substrates
END
;

die "Usage: $usage\n" unless @ARGV == 2;
my ($faafile, $ecfile) = @ARGV;
die "No such file: $faafile\n" unless -e $faafile;
die "No such file: $ecfile\n" unless -e $ecfile;

my %ecinfo = (); # protein id => list of [ EC number, description ]
open(EC, "<", $ecfile) || die "Cannot read $ecfile";
while(<EC>) {
  chomp;
  s/\r//g; # not sure why some lines have these, even within them
  my ($id, undef, $ec, $desc) = split /\t/, $_;
  die "Not enough fields or description is empty: $_\n" unless $desc;
  # Allow duplicate entries for some entries; in practice these often have the same ec info but different families?
  push @{ $ecinfo{$id} }, [ $ec, $desc ];
}
close(EC) || die "Error reading $ecfile";

my $seqs = pbutils::ReadFasta($faafile);
my $nOut = 0;
print STDERR "Read " . scalar(keys %$seqs) . " sequences\n";
my @unknown;
my @skip;
my %known = ();
while (my ($name,$seq) = each %$seqs) {
  my ($id,$family,$ec) = split /[|]/, $name;
  # an identifier may have more than one ec number -- but we will get that information from %ecinfo
  die "Invalid identifier $name in $faafile" unless $ec;
  if (!exists $ecinfo{$id}) {
    push @unknown, $id;
    next;
  }
  next if exists $known{$id};
  my @ecinfo = @{ $ecinfo{$id} };
  my @ec = uniq( sort map { $_->[0] } @ecinfo );
  my @desc = uniq( sort map { $_->[1] } @ecinfo );
  # db, id, id2, shortname, desc, organism, sequence, comment, pmids
  my $comb_desc = join("; ", @desc);
  my $comb_ec = join("; ", @ec);
  if ($comb_desc =~ m/frameshift|fragment/) {
    push @skip, $id;
  } else {
    print join("\t", "CAZy", $id, "", "",
               "$comb_desc (EC $comb_ec)",
               "",
               uc($seq),
               "", "")."\n";
    $nOut++;
  }
  $known{$id} = 1;
}

print STDERR "Skipped " . scalar(@unknown) . " proteins in the fasta file but not the ecinfo file\n";
print STDERR "Skipped " . scalar(@skip) . " proteins that may be fragments or frameshifts\n";
print STDERR "Wrote $nOut different proteins\n";

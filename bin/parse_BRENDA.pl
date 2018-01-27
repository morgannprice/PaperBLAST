#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use pbutils; # for ReadFastaEntry

my $usage = <<END
Usage: parse_BRENDA.pl uniprot_fasta trembl_fasta < brenda_download.txt

The input should be the brenda_download.txt file in a UniProt-like
format.

The output is in curated_parsed format with fields
database id ("BRENDA")
protein identifier (the uniprot id)
secondary identifier (blank)
short name (blank)
description (the recommended name, along with subunit number information if an oligomer)
organism
sequence
comment (blank)
pubmedids (comma delimited)
END
;

die $usage unless @ARGV > 0;
my @faa = @ARGV;
foreach my $file (@faa) {
  die "No such file: $file\n" unless -e $file;
}

# Read an entire record at a time
my $nEntry = 0;
local $/ = "\n///\n";

my %acc = (); # uniprot or swissprot accession to list of [ organism, description, ec#, comma-delimited pubmedids ]
# (potentially more than one EC # so more than one entry)

while(my $record =<STDIN>) {
  # continuation lines start with a tab
  $record =~ s/\n\t/ /g;
  my @lines = split /\n/, $record;
  my $id = "";
  my $enzname;
  my %rf; # number to string
  my %pr; # number to a hash of organism, protdb, protids (a list), and refnos (a list)
  foreach my $line (@lines) {
    next unless $line =~ m/^([A-Z]+)\t(.*)$/;
    my ($field, $value) = ($1,$2);
    if ($field eq "ID") {
      die "Duplicate identifier for $id" if $id;
      # ignore transferred records -- these should have no proteins anyway
      next if $value =~ m/tr[ae]nsferred/ || $value =~ m/deleted/;
      $id = $value;
      # The one time I checked, a "preliminary" EC number was valid, so allow those
      $id =~ s/ *[(]preliminary.*[)]$//;
      $id =~ s/ *[(]reinstated.*[)]$//;

      # Not sure why there is a trailing () in some cases
      $id =~ s/ *[(][)] *$//;
      die "Invalid identifier $value" unless $id =~ m/^[0-9][0-9.]+[0-9]$/;
    } elsif ($field eq "PR") {
      # a PR line includes a number, i.e. #78#, an organism name,
      # a swiss or uniprot identifier formatted like "Q9ZT38 UniProt" or "Q8U259 SwissProt",
      # an optional gene name like (#120# At5g56760, SERAT1.1 <211>) [which is ignored for now]
      # and references like <127,236,239>

      # The sequence identifiers may be mulitples, i.e. P29703 AND P22007 -- this is not handled yet

      # First, remove the number from the beginning
      $value =~ m/^#([0-9]+)# (.*)$/ || die "Cannot parse PR line $value";
      my ($prnum, $line) = ($1,$2);
      # Not sure why there are duplicate PR entries. In general they give redundant entries, i.e. pointing
      # at both SwissProt and UniProt with the same id.
      next if exists $pr{$prnum};
      # ignore entries with no reference and no SwissProt or UniProt identifier
      next unless $line =~ m/SwissProt|UniProt/ && $line =~ m/<[0-9, ]+>/;
      $line =~ m/^(.*)<([0-9, ]+)>$/ || die "Cannot parse references in PR line $value for id $id (prnum $prnum)";
      my ($rest,$refnos) = ($1,$2);
      $refnos =~ s/ /,/g; # necessary if was broken across multiple lines
      my @refnos = split /,/, $refnos;
      @refnos = map { s/^ +//; s/ +$//; $_ } @refnos;
      @refnos = grep { $_ ne "" } @refnos;
      die "Empty reference list for PR line $prnum and id $id" if @refnos == 0;
      foreach my $refno (@refnos) {
        die "Cannot parse references $refno in PR line $prnum for id $id"
          unless $refno =~ m/^\d+$/;
      }
      # note that spaces were not removed from end above so should still be there, regardless of
      # if there is a trailing name or not
      # In rare cases the uniprot id contains lower-case characters but is correct
      unless ($rest =~ m/^(.*) ([A-Za-z][A-Za-z0-9_]+[0-9]) (SwissProt|UniProt) /) {
        # In rare cases, there is a UniProt statement but no UniProt id. I.e. the PR field is
        # PR      #114# Phaeobacter inhibens  UniProt <129>
        # or the identifier is not in UniProt format (mangled?), i.e.
	# "#41# Lysinibacillus sphaericus C5IFUO UniProt <73>"
        # Or, there is a comment about the UniProt taxonomy but there is no UniProt id.
        print STDERR "Warning: cannot parse remaining $rest for id $id PR $prnum\n";
        next;
      }

      my ($org, $protid, $db) = ($1,$2,$3);
      $protid = uc($protid);
      my @protid = ($protid);
      # And handle AND entries, which will show up as strings of " uniprotid AND" at the end of $org
      while($org =~ s/\s+(\S+)\s+AND\s*$//i) {
        my $newid = $1;
        die "Invalid protein id $newid for id $id PR $prnum\n" unless $newid =~ m/^[A-Z][A-Z0-9_]+$/;
        push @protid, $newid;
      }
      die "Cannot parse $value for id $id PR $prnum" if $org =~ m/\t/;
      $pr{$prnum} = { 'organism' => $org, 'protdb' => $db, 'protids' => \@protid, 'refnos' => \@refnos };
    } elsif ($field eq "RF") {
      $value =~ m/^<(\d+)> (.*)$/ || die "Cannot parse RF line: $value";
      my ($refno, $string) = ($1, $2);
      die "Duplicate entry for reference $refno in entry $id" if exists $rf{$refno};
      $rf{$refno} = $string;
    } elsif ($field eq "RN") {
      $enzname = $value;
    }
  } # end handling fields
  next if $id eq ""; # ignore transferred entries and such
  next if keys(%pr) == 0; # ignore if no characterized enzymes
  die "No enzyme name (RN field) for enzyme $id" unless $enzname;
  # print out all members of this enzyme class
  while (my ($prnum, $enz) = each %pr) {
    my $protids = $enz->{protids};
    die unless @$protids > 0;
    my $n = scalar(@$protids);
    my @pmids = ();
    foreach my $refno (@{ $enz->{'refnos'} }) {
      die "Invalid rf $refno for enzyme $id PR $prnum" unless exists $rf{$refno};
      my $rf = $rf{$refno};
      if ($rf =~ m/{Pubmed:(\d+)}/i) {
        push @pmids, $1;
      }
    }
    @pmids = sort { $a <=> $b } @pmids; # oldest first
    if (@pmids > 0) {
      foreach my $i (0..($n-1)) {
        my $i1 = $i+1;
        my $desc = $enzname;
        $desc .= " (subunit $i1/$n)" if $n > 1;
        push @{ $acc{ $protids->[$i] } }, [ $enz->{organism}, $desc, $id, join(",",@pmids) ];
        $nEntry++;
      }
    }
  }
}
print STDERR "Found $nEntry EC : UniProt links with literature support in BRENDA, for " . scalar(keys %acc) . " accessions\n";

local $/ = "\n";

my $nPrint = 0;
foreach my $faafile (@faa) {
  my $fh;
  if ($faafile =~ m/[.]gz$/) {
    open($fh, "zcat $faafile |") || die "Cannot run zcat on $faafile";
  } else {
    open($fh, "<", $faafile) || die "Cannot read $faafile";
  }
  my $state = {};
  while (my ($header,$seq) = ReadFastaEntry($fh,$state)) {
    # ignore non-full-length sequences because they are difficult to interpret for annotation
    next if $header =~ m/[(]Fragment[)]/;
    my ($ids) = split / /, $header;
    my @ids = split /[|]/, $ids;
    shift @ids; # ignore sp or tr at front
    die "Invalid header $header" if @ids == 0;
    foreach my $protid (@ids) {
      next unless exists $acc{$protid};
      my @descs = ();
      my %pmids = ();
      # Used to exclude duplicate EC #s
      my %ecSoFar = ();
      foreach my $row (@{ $acc{$protid} }) {
        my ($org, $desc, $ec, $pmids) = @$row;
        next if exists $ecSoFar{$ec};
        $ecSoFar{$ec} = 1;
        push @descs, "$desc (EC $ec)";
        foreach my $pmid (split /,/, $pmids) { $pmids{$pmid} = 1; }
      }
      my $org = $acc{$protid}[0][0];
      my @pmids = sort {$a <=> $b} keys %pmids;
      print join("\t", "BRENDA", $protid,
                   "", "", # secondary identifier and short name
                   join("; ", @descs), $org, $seq,
                   "", # comment
                   join(",", @pmids))."\n";
      $nPrint++;
    }
  }
  close($fh) || die "Error reading $faafile\n";
}
print STDERR "Wrote EC information from BRENDA for $nPrint accessions\n";

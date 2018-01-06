#!/usr/bin/perl -w
use strict;
use DBI;
use Getopt::Long;

my $usage = <<END
compare_dbs.pl -old old.db -new new.db -out comparison

    This script will compare two databases and will make a table of
    the paper-gene mappings that are present in only one database or
    the other.
END
;

my ($olddb, $newdb, $out);
die $usage unless GetOptions('old=s' => \$olddb,
                             'new=s' => \$newdb,
                             'out=s' => \$out)
  && @ARGV == 0
  && defined $olddb && defined $newdb && defined $out;

my $oldh = DBI->connect("dbi:SQLite:dbname=$olddb","","",{ RaiseError => 1 }) || die $DBI::errstr;
my $newh = DBI->connect("dbi:SQLite:dbname=$newdb","","",{ RaiseError => 1 }) || die $DBI::errstr;

# For simplicity, maps links by pubmed ids, and ignore papers without a pmId

my $nOldOnly = 0;
my $nNewOnly = 0;
open(OUT, ">", $out) || die "Cannot write to $out";

# Mined links:

# geneId => pmId => 1
my %old = ();
my %new = ();

my $oldhits = $oldh->selectall_arrayref(qq{ SELECT geneId, pmId FROM GenePaper WHERE pmId <> "" });
foreach my $row (@$oldhits) {
  my ($geneId, $pmId) = @$row;
  $old{$geneId}{$pmId} = 1;
}

my $newhits = $newh->selectall_arrayref(qq{ SELECT geneId, pmId FROM GenePaper WHERE pmId <> "" });
foreach my $row (@$newhits) {
  my ($geneId, $pmId) = @$row;
  $new{$geneId}{$pmId} = 1;
}

while (my ($geneId, $hashOld) = each %old) {
  foreach my $pmId (keys %$hashOld) {
    if (!exists $new{$geneId}{$pmId}) {
      print OUT join("\t", "old", $geneId, $pmId)."\n";
      $nOldOnly++;
    }
  }
}

while (my ($geneId, $hashNew) = each %new) {
  foreach my $pmId (keys %$hashNew) {
    if (!exists $old{$geneId}{$pmId}) {
      print OUT join("\t", "new", $geneId, $pmId)."\n";
      $nNewOnly++;
    }
  }
}

# And similarly for curated links; use db::protId as the gene id
%old = ();
%new = ();
$oldhits = $oldh->selectall_arrayref(qq{ SELECT db, protId, pmId FROM CuratedPaper });
foreach my $row (@$oldhits) {
  my ($db, $protId, $pmId) = @$row;
  $old{$db . "::" . $protId}{$pmId} = 1;
}

$newhits = $newh->selectall_arrayref(qq{ SELECT db, protId, pmId FROM CuratedPaper });
foreach my $row (@$newhits) {
  my ($db, $protId, $pmId) = @$row;
  $new{$db . "::" . $protId}{$pmId} = 1;
}

while (my ($geneId, $hashOld) = each %old) {
  foreach my $pmId (keys %$hashOld) {
    if (!exists $new{$geneId}{$pmId}) {
      print OUT join("\t", "old", $geneId, $pmId)."\n";
      $nOldOnly++;
    }
  }
}

while (my ($geneId, $hashNew) = each %new) {
  foreach my $pmId (keys %$hashNew) {
    if (!exists $old{$geneId}{$pmId}) {
      print OUT join("\t", "new", $geneId, $pmId)."\n";
      $nNewOnly++;
    }
  }
}

close(OUT) || die "Error writing to $out\n";

print STDERR "Wrote $nOldOnly old-only and $nNewOnly new-only entries to $out\n";


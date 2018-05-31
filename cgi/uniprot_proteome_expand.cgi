#!/usr/bin/perl -w
use strict;
use CGI;
use DBI;
my $limit = 50;
my $cgi = CGI->new();
print $cgi->header;
my $query = $cgi->param('query');
if (defined $query && $query ne "") {
  my $dbfile = "../static/uniprot_proteomes.db";
  die "No such file: $dbfile\n" unless -e $dbfile;
  my $dbh = DBI->connect("dbi:SQLite:dbname=$dbfile","","",{ RaiseError => 1 })
    || die $DBI::errstr;
  my $limit2 = $limit + 1;
  my $hits = $dbh->selectall_arrayref(qq{ SELECT name, isref, isrep FROM Proteome
                                          WHERE name LIKE ?
                                          LIMIT $limit2 }, {}, $query . "%");
  my @hits = @$hits;
  # If many hits, choose the interesting ones
  if (@hits > $limit) {
    my $hitsHi = $dbh->selectall_arrayref(qq{ SELECT name, isref, isrep FROM Proteome
                                               WHERE name LIKE ?
                                               AND (isref = "true" OR isrep = "true")
                                               LIMIT $limit }, {}, $query . "%");
    push @hits, @$hitsHi;
  }
  @hits = sort { $b->[1] cmp $a->[1]
                   || $b->[2] cmp $a->[2]
                     || $a->[0] cmp $b->[0] } @hits;
  my $nShown = 0;
  my %seen = ();
  foreach my $hit (@hits) {
    next if exists $seen{$hit->[0]};
    $seen{$hit->[0]} = 1;
    print join("\t", @$hit)."\n";
    $nShown++;
    last if $nShown >= $limit;
  }
}

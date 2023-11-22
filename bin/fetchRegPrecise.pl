#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use DBI;
use LWP::Simple;

my $mohost = "pub.microbesonline.org";
my $usage = <<END
Usage: fetchRegPrecise.pl -out regprecise > regprecise.curated_parsed

Saves the tables of binding sites to the regprecise/ directory.

Optional arguments:
-host $mohost -- the MicrobesOnline database host
END
;

my $outDir;
die $usage
  unless GetOptions('out=s' => \$outDir,
                    'host=s' => \$mohost)
  && @ARGV == 0
  && defined $outDir;

die "No such directory: $outDir\n" unless -d $outDir;

my ($db, $user, $passwd);
if ($mohost eq "localhost") {
  $db = "genomics_test";
  $user = "test";
  $passwd = "test";
} else {
  $db = "genomics";
  $user = "guest";
  $passwd = "guest";
}
my $dbh = DBI->connect("dbi:mysql:$db:$mohost", $user, $passwd)
  || die $DBI::errstr;

my $regulons = $dbh->selectall_arrayref(
   "SELECT DISTINCT regulonId,regulatorName,regulatorLocusId FROM Locus2RegPrecise;",
   { Slice => {} });

foreach my $regulon (@$regulons) {
  my $regulonId = $regulon->{regulonId};
  my $locusId = $regulon->{regulatorLocusId};
  # Fetch its sequence
  my ($seq) = $dbh->selectrow_array(qq{ SELECT sequence FROM Locus JOIN AASeq USING (locusId,version)
                                        WHERE priority = 1 AND locusId = ?; },
                                   {}, $locusId);
  # Get the locus tag (synonym type 1)
  my ($locusTag) = $dbh->selectrow_array(qq{ SELECT name FROM Synonym JOIN Locus USING (locusId,version)
                                             WHERE priority = 1 AND locusId = ? },
                                         {}, $locusId);
  $locusTag = "VIMSS".$locusId unless $locusTag;

  # Get the organism information
  my ($taxName) = $dbh->selectrow_array(qq{ SELECT shortName FROM Locus JOIN Scaffold USING (scaffoldId)
                                            JOIN Taxonomy USING (taxonomyId)
                                            WHERE priority = 1 AND locusId = ? AND isActive=1 },
                                        {}, $locusId);

  unless ($seq && $taxName) {
    print STDERR "Skipping $locusId $regulon->{regulatorName}\n";
    next;
  }

 # Svae the sites
  my $URL = "https://regprecise.lbl.gov/ExportServlet?type=site&regulonId=$regulonId";
  my $sites = get($URL);
  my $sitesFile = "$outDir/$regulonId.fna";
  if (!$sites) {
    print STDERR "Skipping regulon $regulonId, no sites\n";
    unlink($sitesFile);
    next;
  } else {
    open(my $fh, ">", $sitesFile) || die "Cannot write to $sitesFile\n";
    print $fh $sites;
    close($fh) || die "Error writing to $sitesFile\n";
  }
  # The regprecise page includes Regulator type, TF locus tag, Biological process, Effector (sometimes)
  # Do I want to scrape those to make a better description?
  $URL = "https://regprecise.lbl.gov/regulon.jsp?regulon_id=$regulonId";
  my $lines = get($URL);
  my $desc = $regulon->{regulatorName};
  unless ($lines) {
    print STDERR "Warning: cannot fetch $URL\n";
  } else {
    my $regulatorType = $1 if $lines =~ m!Regulator type:</td>\s*<td>([^<]+)</td>!;
    my $process = $1 if $lines =~ m!Biological process:</td>\s*<td>([^<]+)</td>!;
    my $mode = $1 if $lines =~ m!Regulation mode:</td>\s*<td>([^<]+)</td>!;
    $mode = "activator/repressor"
      if $mode && $mode =~ m/activator/ && $mode =~ m/repressor/;
    my $effector = $1 if $lines =~ m!Effector:</td>\s*<td>([^<]+)</td>!;
    $desc = "$regulon->{regulatorName}";
    $desc .= " regulator of $process" if $process;
    $desc .= ", effector $effector" if $effector;
    $desc .= " ($mode)" if $mode;
  }

  print join("\t", "regprecise", $regulonId, "VIMSS".$locusId,
             $regulon->{regulatorName},
             $desc, $taxName, $seq,
             # no comment or pubmedIds
             "", "")."\n";

}


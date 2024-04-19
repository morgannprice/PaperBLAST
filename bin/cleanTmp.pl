#!/usr/bin/perl -w
use strict;
use FindBin qw{$RealBin};
use Getopt::Long;
use lib "$RealBin/../lib";
use pbutils qw{ReadTable};

my $dir = "$RealBin/../tmp";
my $days = 120;

my $usage = <<END
Usage: cleanTmp.pl [ -tmp dir ] [ -days $days ]

Makes a list of commands to clean out subdirectories and files in
PaperBLAST's tmp directory that have not been accessed in more than
$days days. Files relevant to more recent GapMind analyses are not
slated for removal.

Default tmp directory: $dir
END
;

my $debug;
die $usage
  unless GetOptions('days=i' => \$days, 'tmp=s' => \$dir, 'debug' => \$debug)
  && @ARGV == 0;

my %oldDir = (); # directories to probably remove
my @checkOrg = (); # orgs file to check for subdirectories to be deleted

die "Invalid directory -- no downloaded subdirectory\n" unless -d "$dir/downloaded";

print "#\n# Commands to clean up data not accessed in $days days from\n# $dir\n";
opendir(my $dh, $dir) || die "Cannot open directory $dir\n";
# GapMind analyses to remove
while (my $subfile = readdir $dh) {
  if (-d "$dir/$subfile" && $subfile =~ m/__/) {
    # A regular directory with a GapMind analysis (gdb__gid), to maybe delete
    opendir(my $subdh, "$dir/$subfile") || die "Cannot open directory $dir/$subfile\n";
    my $old = 1;
    my $file;
    while ($old && ($file = readdir $subdh)) {
      next if $file eq "." || $file eq "..";
      my $path = "$dir/$subfile/$file";
      $old = 0 if (-A $path) <= $days;
    }
    closedir($subdh) || die "Error reading $dir/$subfile\n";
    if ($old) {
      $oldDir{$subfile} = 1;
    } else {
      print "# Leaving accessed directory: $dir/$subfile\n";
      push @checkOrg, "$dir/$subfile/orgs.org"
        if -e "$dir/$subfile/orgs.org";
    }
  } elsif (-d "$dir/$subfile" && -e "$dir/$subfile/orgs.org") {
    # A custom GapMind analysis
    push @checkOrg, "$dir/$subfile/orgs.org";
  }
}
closedir($dh) || die "Error reading $dir\n";

my %keep = (); # gdb => gid => 1

foreach my $file (@checkOrg) {
  my @rows = ReadTable($file, ["gdb","gid"]);
  foreach my $row (@rows) {
    $keep{ $row->{gdb} }{ $row->{gid} } = 1;
  }
}
if (defined $debug) {
  foreach my $gdb (sort keys %keep) {
    print "Keep $gdb " . scalar(keys %{ $keep{$gdb} }) . "\n";
  }
}

# Files to remove from downloaded/
opendir (my $ddh, "$dir/downloaded") || die "Cannot open directory $dir/downloaded";
while(my $dfile = readdir $ddh) {
  next if $dfile eq "." || $dfile eq "..";
  my $isDir = -d $dfile;
  my ($gdb, $gid);
  if ($dfile =~ m/^mogenome_(\d+)/) {
    $gdb = "MicrobesOnline";
    $gid = $1;
  } elsif ($dfile =~ m/^uniprot_([A-Za-z0-9-]+)[.]/) {
    $gdb = "UniProt";
    $gid = $1;
  } elsif ($dfile =~ m/^refseq_([A-Z]+_\d+[.]\d+)[.]/) {
    $gdb = "NCBI";
    $gid = $1;
  } elsif ($dfile =~ m/^[0-9a-z]+$/ && $isDir) {
    $gdb = "local";
    $gid = $dfile;
  } elsif ($isDir) {
    $gdb = "IMG";
    $gid = $dfile;
  }
  next unless defined $gdb && defined $gid && $gid ne "";
  unless (exists $keep{$gdb}{$gid} || -A "$dir/downloaded/$dfile" <= $days) {
    if ($isDir) {
      print "rm -Rf $dir/downloaded/$dfile\n";
    } else {
      print "rm $dir/downloaded/$dfile\n";
    }
  }
}

foreach my $subdir (sort keys %oldDir) {
  print "rm -Rf $dir/$subdir\n";
}

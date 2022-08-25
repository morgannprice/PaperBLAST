#!/usr/bin/perl -w
use strict;
use FindBin qw{$RealBin};
use Getopt::Long;
my $dir = "$RealBin/../tmp";
my $days = 120;

my $usage = <<END
Usage: cleanTmp.pl [ -tmp dir ] [ -days $days ]

Makes a list of commands to clean out subdirectories of PaperBLAST's
tmp directory that have not been cleaned out in more than $days $days.

Default tmp directory: $dir
END
;

die $usage
  unless GetOptions('days=i' => \$days, 'tmp=s' => \$dir)
  && @ARGV == 0;

my %oldDir = (); # directories to probably remove
my @checkOrg = (); # orgs file to check for subdirectories to be deleted

print "#\n# Commands to clean up data not accessed in $days days from\n# $dir\n";
opendir(my $dh, $dir) || die "Cannot open directory $dir\n";
while (my $subfile = readdir $dh) {
  if (-d "$dir/$subfile") {
    if ($subfile =~ m/^local__/
        || $subfile =~ m/^NCBI__/
        || $subfile =~ m/^IMG__/
        || $subfile =~ m/_FD$/
        || ($subfile =~ m/^[0-9a-f]+$/ && length($subfile) >= 30)) {
      # A regular directory, to maybe delete
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
      }
    } else {
      # A special directory, does it have an orgs.org file?
      my $orgFile = "$dir/$subfile/orgs.org";
      push @checkOrg, $orgFile
        if $subfile ne "." && $subfile ne ".." && -e $orgFile;
    }
  } else {
    # not a subdirectory
    if (-A "$dir/$subfile" > $days) {
      print "rm -f $dir/$subfile\n";
    } else {
      print "# Leaving accessed file: $dir/$subfile\n";
    }
  }
}
closedir($dh) || die "Error reading directory $dir\n";

foreach my $file (@checkOrg) {
  if (open(my $fh, "<", $file)) {
    while (my $line = <$fh>) {
      chomp $line;
      $line =~ s/\t.*//;
      my @subdir = ($line);
      if ($line =~ m/^local__/) {
        $line =~ s/^local__//;
        push @subdir, $line;
      }
      foreach my $subdir (@subdir) {
        if (exists $oldDir{$subdir}) {
          print "# Leaving directory due to reference: $dir/$line\n";
          delete $oldDir{$line};
        }
      }
    }
  }
}

foreach my $subdir (sort keys %oldDir) {
  print "rm -Rf $dir/$subdir\n";
}

    

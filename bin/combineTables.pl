#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my $usage = <<END
combineTables.pl -in table1 ... tableN > tableCombined
  Combines one or more tables with headers.
  Requires that all input table have identical header lines.
END
;

my @infiles;
die $usage
  unless GetOptions('in=s{1,}' => \@infiles)
  && @ARGV == 0;

my $headerFirst;
foreach my $file (@infiles) {
  open(my $fh, "<", $file) || die "Cannot read $file\n";
  my $headerThis = <$fh>;
  die "No header line in $file\n"
    unless defined $headerThis;
  if (defined $headerFirst) {
    die "Header line for $file does not match $infiles[0]\n"
      unless $headerThis eq $headerFirst;
  } else {
    print $headerThis;
    $headerFirst = $headerThis;
  }
  while(my $line = <$fh>) {
    print $line;
  }
  close($fh) || die "Error reading $file";
}

#!/usr/bin/perl -w
use strict;
use FindBin qw{$Bin};
use lib "$Bin/../lib";
use pbutils qw{SQLiteLine};

my $usage = <<END
Usage: Run as a filter on components.cif
Writes a tab-delimited file suitable for
import with sqlite3, with the fields
  ligandId, ligandDesc
END
;

die $usage if @ARGV > 0;

my $id = undef;
my %idSeen = ();
while(my $line = <STDIN>) {
  chomp $line;
  if ($line eq "#" || $line =~ m/^data/) {
    $id = undef;
  } elsif ($line =~ m/^_chem_comp[.]id +(\S+) *$/) {
    $id = $1;
    die "Duplicate id $id" if exists $idSeen{$id};
    $idSeen{$id} = 1;
  } elsif ($line =~ m/^_chem_comp[.] /) {
    die "Cannot handle split id line $line";
  } elsif ($line =~ m/^_chem_comp.name +(.*)$/) {
    my $name = $1;
    die "name line $line with no preceding id" unless defined $id;
    die "Duplicate name for id $id: $line" if $idSeen{$id} > 1;
    $idSeen{$id} = 2;
    $name =~ s/ *$//;
    if ($name eq "") {
      # empty name line implies following lines have the value
      # This could be a single line with the value
      # Or a line starting with "; and the value, followed by a line with ;
      # Or a line with just ";", then line(s) with the value, then a line with just ";"
      # or a single line has the value
      my $n = 0;
      my $inSemi = 0; # set if the first line started with ; (so we need to look for the next one)
      while ($line = <STDIN>) {
        $n++;
        chomp $line;
        my $hasSemi = 0;
        if ($line =~ m/^;/) {
          $hasSemi = 1;
          $line =~ s/^;//;
        }
        $inSemi = 1 if $n == 1 && $hasSemi;
        $name .= $line;
        last if !$inSemi || ($n > 1 && $hasSemi);
      }
      die "Still an empty name for $id after $n lines, inSemi = $inSemi" if $name eq "";
    }
    # occasional tab characters within names -- convert to white space
    $name =~ s/\s/ /g;
    # strip any leading or trailing whitespace
    $name =~ s/ *$//;
    $name =~ s/^ *//;
    if ($name =~ m/^"([^"]+)"$/) {
      $name = $1;
    } elsif ($name =~ m/^'([^']+)+'$/) {
      $name = $1;
    } elsif ($name =~ m/["]/) {
      print STDERR "Warning: quote characters in $name for $id\n";
    }
    print SQLiteLine($id,$name);
  } elsif ($line =~ m/^_chem_comp.name /) {
    die "Not sure how to handle name from $line";
  }
}


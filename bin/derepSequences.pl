#!/usr/bin/perl -w
use strict;
use FindBin qw{$Bin};
use Getopt::Long;
use IO::Handle; # for autoflush

my $usage = <<END
derepSequences.pl -dir dir

Given a directory that contains litsearch.faa and the sql database
(litsearch.db), compute the dereplicated sequences (uniq.faa) and
populate the SeqToDuplicate table.

Optionally -- with the -test argument, does not modify the database
END
    ;

my $test;
my $dir;
die $usage unless GetOptions('dir=s' => \$dir,
                             'test' => \$test)
    && @ARGV == 0
    && defined $dir;
die "Not a directory: $dir\n" unless -d $dir;
my $faaIn = "$dir/litsearch.faa";
die "No such file: $faaIn\n" unless -e $faaIn;
my $faaOut = "$dir/uniq.faa";
my $dbfile = "$dir/litsearch.db";
die "No such file: $dbfile\n" unless -e $dbfile;

my $usearch = "$Bin/usearch";
die "No such executable: $usearch\n" unless -x $usearch;

unlink("$dir/uc.derep");
system($usearch, "-quiet", "-fastx_uniques", $faaIn, "-uc", "$dir/uc.derep", "-fastaout", $faaOut) == 0
    || die "Error running usearch -fastx_uniques: $!";

open(UC, "<", "$dir/uc.derep") || die "Cannot read $dir/uc.derep";
open (OUT, ">", "$dir/SeqToDuplicate") || die "Cannot write to $dir/SeqToDuplicate";
my %seen = ();
while(my $line = <UC>) {
    chomp $line;
    my ($type, $cluster_id, undef, undef, undef, undef, undef, undef, $cluster_member, $cluster_name) = split /\t/, $line;
    die "Not enough fields in $line" unless defined $cluster_name;
    # Sometimes UniProt sequences are duplicated with the same accession -- ignore those duplicates
    # or duplicate matches to another sequence
    if ($type eq "H" && $cluster_name ne $cluster_member) {
        print OUT join("\t", $cluster_name, $cluster_member)."\n"
            unless exists $seen{$cluster_member};
        $seen{$cluster_member} = 1;
    }
}
close(UC) || die "Error reading from $dir/uc.derep";
close(OUT) || die "Error writing to $dir/SeqToDuplicate";

if (defined $test) {
    print STDERR "Test mode -- not writing to the database\n";
    exit(0);
}
# else
open(SQLITE, "|-", "sqlite3", $dbfile) || die "Cannot run sqlite3 on $dbfile";
autoflush SQLITE 1;
print SQLITE ".mode tabs\n";
my @tables = qw{SeqToDuplicate};
foreach my $table (@tables) {
    print SQLITE "DELETE FROM $table;\n";
    print SQLITE ".import $dir/$table $table\n";
}
    my $reportCmds = <<END
SELECT 'nDup', COUNT(*) FROM SeqToDuplicate;
END
    ;
print SQLITE $reportCmds;
close(SQLITE) || die "Error running sqlite3 commands\n";
print STDERR "Loaded SeqToDuplicate\n";


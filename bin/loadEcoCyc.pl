#!/usr/bin/perl -w
# This script is obsolete; see parse_ecocyc.pl instead
use strict;
use IO::Handle; # for autoflush
use Getopt::Long;

my $usage = <<END
loadEcoCyc.pl -dir dir -tab ecocyc.tab -seq protseq.fsa

Fills the EcoCyc and EcoCycToPubMed tables in dir/litsearch.db
Puts the sequences in protseq.fsa, in modified format, into
dir/ecocyc.faa

Optional argument: -test -- read the input and makes the files in dir
for importing into the database, but does not modify either database.

END
    ;

my $dir;
my $tabfile;
my $seqfile;
my $test;

die $usage unless GetOptions('dir=s' => \$dir,
                             'tab=s' => \$tabfile,
                             'seq=s' => \$seqfile,
                             'test' => \$test )
    && @ARGV == 0
    && defined $dir && defined $tabfile && defined $seqfile;
die "No such directory: $dir\n" unless -d $dir;
my $dbfile = "$dir/litsearch.db";
my $faaout = "$dir/ecocyc.faa";
die "No such file: $dbfile\n" unless -e $dbfile;
die "No such file: $tabfile\n" unless -e $tabfile;
die "No such file: $seqfile\n" unless -e $seqfile;

my %seq = (); # id like "gcl|ECOLI|PROTEIN-ID" to sequence
open(SEQ, "<", $seqfile) || die "Cannot read $seqfile";
my $seqname;
while(my $line = <SEQ>) {
    chomp $line;
    $line =~ s/\r$//; # an occasional issue
    if ($line =~ m/^>(\S+) /) {
        $seqname = $1;
        die "Duplicate sequence name $seqname" if exists $seq{$seqname};
        $seq{$seqname} = "";
        die "Invalid sequence name $seqname" unless $seqname =~ m/^gnl[|]ECOLI[|]./;
    } else {
        die "Invalid sequence line $line" unless $line =~ m/^[A-Z*]+$/;
        $line =~ s/[*]//g; # happens in a few pseudogenes
        die "No header in fasta file $seqfile" unless defined $seqname;
        $seq{$seqname} .= $line;
    }
}
close(SEQ) || die "Error reading $seqfile";
print STDERR "Read " . scalar(keys %seq) . " sequences from $seqfile\n";

my %prot = (); # protein_id to list of bnumber, short name, description
my %protPub = (); # protein to list of pubmed ids
my %protLen = (); # protein_id => length, or missing if it has no sequence
my %protSkip = (); # protein_id => 1 if missing
open(TAB, "<", $tabfile) || die "Error reading $tabfile";
while(my $line = <TAB>) {
    chomp $line;
    my ($protein_id, $bnumber, $name, $desc, $pmid_string) = split /\t/, $line, -1;
    die "Not enough fields in $line" unless defined $pmid_string;
    die "Duplicate entry for $protein_id in $tabfile" if exists $prot{$protein_id};
    $prot{$protein_id} = [ $bnumber, $name, $desc ];
    my @pmids = split /,/, $pmid_string;
    $protPub{$protein_id} = \@pmids;
    my $seq_id = "gnl|ECOLI|" . $protein_id;
    if (exists $seq{$seq_id}) {
        $protLen{$protein_id} = length($seq{$seq_id});
    } else {
        $protSkip{$protein_id} = 1;
    }
}
close(TAB) || die "Error reading $tabfile";
print STDERR "Read " . scalar(keys %prot) . " rows from $tabfile, skipping " . scalar(keys %protSkip) . " proteins with no sequence\n";

# Set up EcoCyc and EcoCycToPubMed
open(OUT, ">", "$dir/EcoCyc") || die "Cannot write to $dir/EcoCyc";
open(PUB, ">", "$dir/EcoCycToPubMed") || die "Cannot write to $dir/EcoCycToPubMed";
foreach my $protein_id (sort keys %prot) {
    next if exists $protSkip{$protein_id};
    print OUT join("\t", $protein_id, @{ $prot{$protein_id} }, $protLen{$protein_id})."\n";
    my $list = $protPub{$protein_id};
    foreach my $pubmedId (@$list) {
        print PUB join("\t", $protein_id, $pubmedId)."\n";
    }
}
close(OUT) || die "Error writing to $dir/EcoCyc";
close(PUB) || die "Error writing to $dir/EcoCycToPubMed";
print STDERR "Wrote $dir/EcoCyc\n";
print STDERR "Wrote $dir/EcoCycToPubMed\n";

# Write the ecocyc.faa file
open(FAA, ">", $faaout) || die "Cannot write to $faaout";
my $nSeqWritten = 0;
while (my ($seq_id, $seq) = each %seq) {
    # only save it if it has metadata
    $seq_id =~ m/^gnl\|ECOLI\|(.*)$/ || die "Bad sequence id $seq_id";
    my $protein_id = $1;
    next unless exists $prot{$protein_id};
    print FAA ">$seq_id\n$seq\n";
    $nSeqWritten++;
}
close(FAA) || die "Error writing to $faaout";
print STDERR "Wrote $nSeqWritten sequences to $faaout\n";
print STDERR "Warning! Wrong number of sequences!\n"
    unless $nSeqWritten == scalar(keys %prot) - scalar(keys %protSkip);

if (defined $test) {
    print STDERR "In test mode -- exiting with no changes to the databases\n";
    exit(0);
}

# Update database
open(SQLITE, "|-", "sqlite3", $dbfile) || die "Cannot run sqlite3 on $dbfile";
autoflush SQLITE 1;
print SQLITE ".mode tabs\n";
my @tables = qw{EcoCyc EcoCycToPubMed};
foreach my $table (@tables) {
    print SQLITE "DELETE FROM $table;\n";
    print SQLITE ".import $dir/$table $table\n";
}
    my $reportCmds = <<END
SELECT 'nEcoCyc', COUNT(*) FROM EcoCyc;
SELECT 'nEcoCycWithPubs', COUNT(DISTINCT protein_id) FROM EcoCycToPubmed;
END
    ;
print SQLITE $reportCmds;
close(SQLITE) || die "Error running sqlite3 commands\n";
print STDERR "Loaded EcoCyc tables\n";

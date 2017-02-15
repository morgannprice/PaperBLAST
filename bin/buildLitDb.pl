#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw{$Bin};
use File::Which;
use IO::Handle; # for autoflush

my $usage = <<END
buildLitDb.pl -dir dir [ -sprot sprot.char.tab ] prefix1 ... prefixN
        [ -snippets snippetsFile] [ -ecocyc datadir ]
	[ -blastdir blast ]

where each prefix is the output of parseEuropePMCHits, and the
(optional) sprot file is the output of sprotCharacterized.pl
snippets file is the output of buildSnippets.pl or combineSnippets.pl
	[the snippets.access file must exist as well]

The -blast argument specifies the path to the blastall and formatdb executables
(or else there must be a blast/ directory of the working directory or the executable directory).

The -ecocyc argument specifies the directory with the protseq.fsa and
proteins.dat files.

Creates dir/litsearch.db, dir/litsearch.faa, and formats the BLAST database
END
    ;

# sqlite3 expects CVS format, not exactly tab delimited format
# So, need to replace any " with "" and surround the field with quotes.
sub csv_quote($);

{
    my ($dir, $sprotchar, $snippetsFile);
    my $ecocycdir;
    my $blastdir;
    die $usage
        unless GetOptions('dir=s' => \$dir,
                          'sprot=s' => \$sprotchar,
                          'snippets=s' => \$snippetsFile,
                          'ecocyc=s' => \$ecocycdir )
        && defined $dir
        && defined $snippetsFile
        && @ARGV > 0;
    die "Not a directory: $dir\n" unless -d $dir;
    die "No such file: $sprotchar" if defined $sprotchar && ! -e $sprotchar;
    die "No such file: $snippetsFile\n" if ! -e $snippetsFile;
    die "No such file: $snippetsFile.access\n" if ! -e "$snippetsFile.access";

    if (defined $ecocycdir) {
        die "No such directory: $ecocycdir\n" unless -d $ecocycdir;
        die "No such file: $ecocycdir/proteins.dat\n" unless -e "$ecocycdir/proteins.dat";
        die "No such file: $ecocycdir/protseq.fsa\n" unless -e "$ecocycdir/protseq.fsa";
    }

    if (!defined $blastdir) {
        if (-d "blast") {
            $blastdir = "blast";
        } elsif (-d "$Bin/blast") {
            $blastdir = "$Bin/blast";
        } else {
            die "No -blastdir argument and not in . or $Bin\n";
        }
    }
    die "No such executable: $blastdir/blastall\n" unless -x "$blastdir/blastall";
    die "No such executable: $blastdir/formatdb\n" unless -x "$blastdir/formatdb";

    my $schema = "$Bin/litsearch.sql";
    die "Cannot find schema: not in $schema" unless -e $schema;

    # Check that all prefix files exist
    my @prefixes = @ARGV;
    foreach my $prefix (@prefixes) {
        die "No such file: $prefix.papers\n" unless -e "$prefix.papers";
        die "No such file: $prefix.queries\n" unless -e "$prefix.queries";
        die "No such file: $prefix.faa\n" unless -e "$prefix.faa";
        
    }

    my $sqldb = "$dir/litsearch.db";
    unlink($sqldb);

    my $faafile = "$dir/litsearch.faa";
    open(UNIPROT, ">", "$dir/UniProt") || die "Cannot write to $dir/Gene";

    system("sqlite3 $sqldb < $schema");

    if (defined $ecocycdir) {
        system("$Bin/parseEcoCycProteins.pl < $ecocycdir/proteins.dat > $dir/ecocyc.tab") == 0
            || die "Error running parseEcoCycProteins.pl: $!";
        system("$Bin/loadEcoCyc.pl", "-dir", $dir, "-tab", "$dir/ecocyc.tab", "-seq", "$ecocycdir/protseq.fsa") == 0
            || die "Error running loadEcoCyc.pl";
    }

    # These are to prevent duplicates
    my %locusLen = (); # geneId => length of protein sequence
    my %aaseq = (); # geneId => sequence if already saw a sequence
    my %isUniprot = (); # geneId => 1 if from uniprot
    if (defined $sprotchar) {
        open(SPROTIN, "<", $sprotchar) || die "Cannot read $sprotchar";
        while (my $line = <SPROTIN>) {
            chomp $line;
            my ($acc, $desc, $organism, $seq, $cc) = split /\t/, $line;
            die "Invalid sprot input:\n$line\n" unless defined $cc && $cc ne "";
            $aaseq{$acc} = $seq;
            $isUniprot{$acc} = 1;
            print UNIPROT join("\t", $acc, $desc, $organism, &csv_quote($cc), length($seq))."\n";
        }
        close(SPROTIN) || die "Error reading $sprotchar";
    }
    close(UNIPROT) || die "Error writing $dir/UniProt";

    open(GENETAG, ">", "$dir/Gene") || die "Cannot write to $dir/Gene";
    foreach my $prefix (@prefixes) {
        open (QUERYIN, "<", "$prefix.queries") || die "Cannot read $prefix.queries";
        while (my $line = <QUERYIN>) {
            chomp $line;
            my ($geneId, $organism, $protlen, $desc) = split /\t/, $line;
            die "Invalid input" unless defined $protlen && $protlen =~ m/^\d+$/;
            die "Mismatching length for $geneId in $prefix"
                if exists $locusLen{$geneId} && $locusLen{$geneId} ne $protlen;
            next if exists $locusLen{$geneId};
            $locusLen{$geneId} = $protlen;
            $desc = "" if !defined $desc;
            print GENETAG join("\t", $geneId, $organism, $protlen, $desc) . "\n";
        }
        close(QUERYIN) || die "Error reading $prefix.queries";

        open(FAAIN, "<", "$prefix.faa");
        while(my $header = <FAAIN>) {
            chomp $header;
            die "Invalid faa file $prefix.faa" unless $header =~ m/^>(.*)$/;
            my $geneId = $1;
            die "Unexpected gene $geneId in $prefix.faa" unless exists $locusLen{$geneId};
            my $seq = <FAAIN>;
            die "Invalid faa file $prefix.faa" unless defined $seq;
            chomp $seq;
            die "Incorrect length for $geneId in $prefix.faa" unless length($seq) == $locusLen{$geneId};
            if (exists $aaseq{$geneId}) {
                print STDERR "Warning: non-matching sequences for $geneId in $prefix.faa\n"
                    if $aaseq{$geneId} ne $seq;
            } else {
                $aaseq{$geneId} = $seq;
            }
        }
        close(FAAIN) || die "Error reading $prefix.faa";

    }
    close(GENETAG) || die "Error writing to $dir/Gene";

    open(OUT, ">", "$dir/Snippet") || die "Cannot write to $dir/Snippet";
    my %hasSnippet = (); # pmcId::pmId => queryId => 1
    my $nSuppressedByCase = 0;
    open(IN, "<", $snippetsFile) || die "Cannot read $snippetsFile";
    while(my $line = <IN>) {
        chomp $line;
        my ($pmcId, $pmId, $queryTerm, $queryId, $snippet) = split /\t/, $line;
        die "Not enough fields in $snippetsFile\n$line" unless defined $snippet;
        if ($queryTerm =~ m/^[a-zA-Z]\d+$/ && $snippet !~ m/$queryTerm/) {
            $nSuppressedByCase++;
        } else {
            print OUT join("\t", $queryId, $queryTerm, $pmcId, $pmId, &csv_quote($snippet))."\n";
            $hasSnippet{ join("::", $pmcId, $pmId) }{ $queryId } = 1;
        }
    }
    print STDERR "Wrote $dir/Snippet, suppressed $nSuppressedByCase snippets with risky locus tags and the wrong case\n";
    close(OUT) || die "Error writing to $dir/Snippet";

    my %full = (); # paperId (pmcId::pmId) => full text
    open(OUT, ">", "$dir/PaperAccess") || die "Cannot write to $dir/PaperAccess";
    open(IN, "<", "$snippetsFile.access") || die "Cannot read $snippetsFile.access";
    my %seen = ();
    while(my $line = <IN>) {
        chomp $line;
        my ($pmcId, $pmId, $access) = split /\t/, $line;
        die "Not enough fields in $snippetsFile.access\n$line" unless defined $access;
        my $key = $pmcId."::".$pmId;
        die "Duplicate row in $snippetsFile.access for $pmcId $pmId"
            if exists $seen{$key};
        $seen{$key} = 1;
        print OUT join("\t", $pmcId, $pmId, $access)."\n";
        $full{$key} = 1 if $access eq "full";
    }
    close(IN) || die "Error reading $snippetsFile.access";
    close(OUT) || die "Error writing to $dir/PaperAccess";

    # store which gene/paper combinations were left in
    my %geneHasPaper = ();
    my $nIgnoreNoSnippet = 0;
    open(GENEPAPER, ">", "$dir/GenePaper") || die "Cannot write to $dir/GenePaper";
    foreach my $prefix (@prefixes) {
        open(PAPERIN, "<", "$prefix.papers") || die "Cannot read $prefix.papers";
        while(my $line = <PAPERIN>) {
            chomp $line;
            my ($queryId, $queryTerm, $pmcId, $pmId, $doi, $title, $authors, $journal, $year, $isOpen) = split /\t/, $line;
            die "Invalid papers line in $prefix.papers:\n$line\n" unless defined $isOpen;
            die "Unexpected gene $queryId in $prefix.papers" unless exists $aaseq{$queryId};
            # Do not bother to remove duplicates
            my $key = join("::", $pmcId, $pmId);
            if (exists $full{$key} && !exists $hasSnippet{$key}{$queryId}) {
                $nIgnoreNoSnippet++;
            } else {
                $geneHasPaper{$queryId} = 1;
                print GENEPAPER join("\t", $queryId, $queryTerm, $pmcId, $pmId, $doi, &csv_quote($title), $authors, $journal, $year, $isOpen)."\n";
            }
        }
        close(PAPERIN) || die "Error reading $prefix.papers";
    }
    close(GENEPAPER) || die "Error writing to $dir/GenePaper";
    print STDERR "Read papers: ignored $nIgnoreNoSnippet cases with full text and no snippet. " . scalar(keys %geneHasPaper) . " cases left\n";

    my $nSkipGene = 0;
    open(FAA, ">", $faafile) || die "Cannot write to $faafile";
    while (my ($geneId, $seq) = each %aaseq) {
        if (exists $geneHasPaper{$geneId} || $isUniprot{$geneId}) {
            print FAA ">$geneId\n$seq\n";
        } else {
            $nSkipGene++;
        }
    }
    close(FAA) || die "Error writing to $faafile";
    print STDERR "Wrote $faafile, with $nSkipGene genes skipped because no remaining links to papers\n";

    open(SQLITE, "|-", "sqlite3", "$sqldb") || die "Cannot run sqlite3 on $sqldb";
    autoflush SQLITE 1;
    print SQLITE ".mode tabs\n";
    my @tables = qw{GenePaper Gene UniProt Snippet PaperAccess};
    foreach my $table (@tables) {
        print STDERR "Loading table $table\n";
        print SQLITE ".import $dir/$table $table\n";
    }

    my $reportCmds = <<END
SELECT 'nGenesWithPubs', COUNT(DISTINCT geneId) FROM Gene;
SELECT 'nPubMedIds', COUNT(DISTINCT pmId) FROM GenePaper;
SELECT 'nPaperLinks', COUNT(*) FROM GenePaper;
SELECT 'nSnippets', COUNT(*) FROM Snippet;
SELECT 'nUniProts', COUNT(*) FROM UniProt;
SELECT 'nFullPapers', COUNT(*) FROM PaperAccess WHERE access="full";
END
    ;
    print SQLITE $reportCmds;
    close(SQLITE) || die "Error running sqlite3 commands\n";

    if (defined $ecocycdir) {
        system("cat $dir/ecocyc.faa >> $dir/litsearch.faa") == 0
            || die "Appending to $dir/litsearch.faa failed: $!";
    }

    system("$Bin/derepSequences.pl", "-dir", $dir) == 0
        || die "derepSequences.pl failed: $!";

    system("$blastdir/formatdb", "-p", "T", "-i", "$dir/uniq.faa", "-o", "T") == 0
        || die "formatdb failed";

    print STDERR "Success\n";
}

sub csv_quote($) {
    my ($in) = @_;
    return $in unless $in =~ m/"/;
    $in =~ s/"/""/g;
    return '"' . $in . '"';
}

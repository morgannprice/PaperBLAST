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
(optional) snippets file is the output of buildSnippets.pl or combineSnippets.pl
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
        && @ARGV > 0;
    die "Not a directory: $dir\n" unless -d $dir;
    die "No such file: $sprotchar" if defined $sprotchar && ! -e $sprotchar;
    die "No such file: $snippetsFile\n" if defined $snippetsFile && ! -e $snippetsFile;
    die "No such file: $snippetsFile.access\n" if defined $snippetsFile && ! -e "$snippetsFile.access";

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

    my $faadb = "$dir/litsearch.faa";
    open(FAA, ">", $faadb) || die "Cannot write to $faadb";
    open(UNIPROT, ">", "$dir/UniProt") || die "Cannot write to $dir/Gene";

    system("sqlite3 $sqldb < $schema");

    if (defined $ecocycdir) {
        system("$Bin/parseEcoCycProteins.pl < $ecocycdir/proteins.dat > $dir/ecocyc.tab") == 0
            || die "Error running parseEcoCycProteins.pl: $!";
        system("$Bin/loadEcoCyc.pl", "-dir", $dir, "-tab", "$dir/ecocyc.tab", "-seq", "$ecocycdir/protseq.fsa") == 0
            || die "Error running loadEcoCyc.pl";
    }

    if (defined $sprotchar) {
        open(SPROTIN, "<", $sprotchar) || die "Cannot read $sprotchar";
        while (my $line = <SPROTIN>) {
            chomp $line;
            my ($acc, $desc, $organism, $seq, $cc) = split /\t/, $line;
            die "Invalid sprot input:\n$line\n" unless defined $cc && $cc ne "";
            print FAA ">$acc\n$seq\n";
            print UNIPROT join("\t", $acc, $desc, $organism, &csv_quote($cc), length($seq))."\n";
        }
        close(SPROTIN) || die "Error reading $sprotchar";
    }
    close(UNIPROT) || die "Error writing $dir/UniProt";

    # These are to prevent duplicates
    my %locusLen = (); # geneId => length of protein sequence
    my %aaseq = (); # geneId => 1 if already put out the sequence
    open(GENEPAPER, ">", "$dir/GenePaper") || die "Cannot write to $dir/GenePaper";
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
            print FAA ">$geneId\n$seq\n" unless exists $aaseq{$geneId};
            $aaseq{$geneId} = 1;
        }
        close(FAAIN) || die "Error reading $prefix.faa";

        open(PAPERIN, "<", "$prefix.papers") || die "Cannot read $prefix.papers";
        while(my $line = <PAPERIN>) {
            chomp $line;
            my ($queryId, $queryTerm, $pmcId, $pmId, $doi, $title, $authors, $journal, $year, $isOpen) = split /\t/, $line;
            die "Invalid papers line in $prefix.papers:\n$line\n" unless defined $isOpen;
            die "Unexpected gene $queryId in $prefix.papers" unless exists $aaseq{$queryId};
            # Do not bother to remove duplicates
            print GENEPAPER join("\t", $queryId, $queryTerm, $pmcId, $pmId, $doi, &csv_quote($title), $authors, $journal, $year, $isOpen)."\n";
        }
        close(PAPERIN) || die "Error reading $prefix.papers";
    }
    close(GENEPAPER) || die "Error writing to $dir/GenePaper";
    close(GENETAG) || die "Error writing to $dir/Gene";
    close(FAA) || die "Error writing to $faadb";

   
    open(OUT, ">", "$dir/Snippet") || die "Cannot write to $dir/Snippet";
    if (defined $snippetsFile) {
        open(IN, "<", $snippetsFile) || die "Cannot read $snippetsFile";
        while(my $line = <IN>) {
            chomp $line;
            my ($pmcId, $pmId, $queryTerm, $queryId, $snippet) = split /\t/, $line;
            die "Not enough fields in $snippetsFile\n$line" unless defined $snippet;
            print OUT join("\t", $queryId, $queryTerm, $pmcId, $pmId, &csv_quote($snippet))."\n";
        }
        
    }
    close(OUT) || die "Error writing to $dir/Snippet";

    open(OUT, ">", "$dir/PaperAccess") || die "Cannot write to $dir/PaperAccess";
    if (defined $snippetsFile) {
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
        }
        close(IN) || die "Error reading $snippetsFile.access";
    }
    close(OUT) || die "Error writing to $dir/PaperAccess";

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

    system("$blastdir/formatdb", "-p", "T", "-i", "$faadb", "-o", "T") == 0
        || die "formatdb failed";

    print STDERR "Success\n";
}

sub csv_quote($) {
    my ($in) = @_;
    return $in unless $in =~ m/"/;
    $in =~ s/"/""/g;
    return '"' . $in . '"';
}

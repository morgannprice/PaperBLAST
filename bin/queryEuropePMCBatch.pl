#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use POSIX; # for ceil
use FindBin qw{$Bin};
use IO::Handle; # for autoflush

my $baseURL = undef;
my $maxPerSec = 4.5;
my $parallel = 10;
my $dir = ".";
my $debug = undef;

my $usage = <<END
queryEuropePMCBatch.pl -in queryprot -out out.tab [ -dir $dir ]
   [ -maxPerSec $maxPerSec ] [ -parallel $parallel ]
   [ -URL alternate_url ]

Splits the input into parts, runs queryEuropePMC.pl on each, combines
the results, and makes a list of failures

The dir argument specifies where to put the temporary files.
END
    ;

sub ForkCommandsAndWait;

{
    my ($infile,$outfile);
    die $usage
        unless GetOptions('in=s' => \$infile,
                          'out=s' => \$outfile,
                          'maxPerSec=f' => \$maxPerSec,
                          'parallel=i' => \$parallel,
                          'URL=s' => \$baseURL,
                          'dir=s' => \$dir,
                          'debug' => \$debug)
        && defined $infile && defined $outfile
        && defined $dir
        && @ARGV == 0;

    die "No such directory: $dir\n" unless -d $dir;
    die "No such file: $infile\n" unless -e $infile;
    die "-parallel must be at least 1" unless $parallel > 0;
    die "No such executable: $Bin/queryEuropePMC.pl" unless -x "$Bin/queryEuropePMC.pl";

    autoflush STDERR 1;

    my ($nLinesIn,undef) = split / /, `wc -l $infile`;
    $parallel = 1 if $nLinesIn <= 1;
    my $minWait = sprintf("%.2f",$parallel / $maxPerSec);
    my $nPerPart = ceil($nLinesIn / $parallel);
    my $pre = "$dir/qepmcb.$$";
    print STDERR "Chopping $nLinesIn input lines into pieces with at least $nPerPart\n";
    print STDERR "Prefix $pre\n";
    system("rm -f $pre.*");
    system("split", "-l", $nPerPart, $infile, "$pre.") == 0
        || die "split failed: $!";
    my @parts  = glob("$pre.*");
    print join("\t", "MadeParts", @parts)."\n" if defined $debug;

    my @cmds = ();
    foreach my $part (@parts) {
        my $cmd = "$Bin/queryEuropePMC.pl -in $part -out $part.out -wait $minWait ";
        $cmd .= " -URL $baseURL" if defined $baseURL;
        $cmd .= " >& $part.log";
        push @cmds, $cmd;
    }

    ForkCommandsAndWait(@cmds);
    system("cat $pre*log");

    open(OUT, ">", $outfile) || die "Cannot write to $outfile\n";
    my %failed = (); # queryId => 1 if should retry
    foreach my $part (@parts) {
        open(PART, "<", "$part.out") || die "Error reading $part.out";
        foreach my $line (<PART>) {
            chomp $line;
            my ($queryId, $json) = split /\t/, $line;
            if ($json =~ m/hitCount/) {
                print OUT $line."\n";
            } else {
                $failed{$queryId} = 1;
            }
        }
        close(PART) || die "Error reading $part.out";
    }
    close(OUT) || die "Error writing to $outfile";
    my $nFail = scalar(keys %failed);
    print STDERR "Wrote $outfile except for $nFail failed queries\n";

    # failed queries
    my $failfile = "$outfile.fail";
    open(FAIL, ">", $failfile) || die "Cannot write to $failfile\n";
    foreach my $queryId (sort keys %failed) {
        print FAIL "$queryId\n";
    }
    close(FAIL) || die "Error writing to $failfile";
    print STDERR "Wrote $failfile\n";
    
    system("rm -f $pre.*");
}

sub ForkCommandsAndWait {
    my @cmds = @_;
    my $nChildren = 0;
    foreach my $cmd (@cmds) {
        my $pid = fork;
        if ($pid == 0) {
            print STDERR "Running $cmd\n" if defined $debug;
            system($cmd) == 0
                || die "Child process failed: $!";
            exit(0);
        }
        # else this is the parent
        $nChildren++;
    }        
    while ($nChildren > 0) {
        $nChildren--;
        my $pid = wait;
        die "No child process\n" if $pid < 0;
        my $status = $? >> 8;
        die "Child process failed: $?" if $status != 0;
    }
    print STDERR "All " . scalar(@cmds) . " children succeeded\n"
        if defined $debug;
}

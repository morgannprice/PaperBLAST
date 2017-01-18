#!/usr/bin/perl -w
# submitter.pl -- submit jobs and wait for them to finish
use strict;
use POSIX qw(sys_wait_h);

my $debug = 0;

my $wait = 2;
my $verbose = 0;
my $nCPUs = $ENV{MC_CORES} || (`egrep -c '^processor' /proc/cpuinfo`);
my $maxload = $ENV{SUBMITTER_MAXLOAD} || $ENV{MC_CORES} || $nCPUs+1;
my $n = $ENV{SUBMITTER_NJOBS} || $ENV{MC_CORES} || int(0.6 * ($nCPUs+1));
$n = 1 if $n < 1;
my $shell = $ENV{SUBMITTER_SHELL} || "ssh";
my $hosts = $ENV{SUBMITTER_HOSTS} || undef;

my $cd = "";
my $logPrefix = "";
my $envFile = "";

my $usage = <<END
Usage: submitter.pl [-shell rsh|ssh] [-host X,Y,Z] [-wait $wait]
    [-n $n ] [ -maxload $maxload ] [-cd dir] [-verbose 0]
    [-log name] [-env file] commandsFile

host is the list of machine(s) to run the jobs on (or this machine if none)

shell (ignored if no host) is ssh by default. You can use the
environment variable SUBMITTER_SHELL instead of -shell to set this.

n is the maximum number of jubs to run in parallel -- you can also use
the environment variable SUBMITTER_NJOBS instead of -n to set this.

maxload is the maximum load on the current machine or the ssh-into
machine. If uptime returns a higher value then it waits. You can set
this parameter with the environment variable SUBMITTER_MAXLOAD instead
of with the -maxload option.

cd is the directory to change to before running each job

-log is the prefix of the filename (relative to the cd dir) to use for each log
(it defaults to commandFile)

env should be a list of lines of the form VARNAME=VALUE (lines starting with ! or # ignored)

The commands in commandsFile will be executed by either
ssh host '(cd dir; var1=value1...; time LINE) >& dir/log-CMDNO' (if host(s) are specified)
or by
(cd dir; var1=value1...; time LINE) >& dir/log-CMDNO

Lines in commandFile that start with ! or # are also ignored.

END
    ;

use Getopt::Long;
die $usage unless GetOptions('debug=i' => \$debug,
			     'shell=s' => \$shell,
			     'host=s' => \$hosts,
			     'wait=i' => \$wait,
			     'n=i' => \$n,
			     'maxload=i' => \$maxload,
			     'cd=s' => \$cd,
			     'verbose=i' => \$verbose,
			     'logPrefix=s' => \$logPrefix,
			     'env=s' => \$envFile) && @ARGV==1;
my $commandsFile = shift @ARGV;
			     
my @commands = ();
my @envNames = ();
my @envValues = ();

my @hosts = (""); # a single host with no name if running jobs locally
@hosts = split /,/, $hosts if defined $hosts;

$shell = "rsh -n" if $shell eq "rsh";

$logPrefix = $commandsFile if ($logPrefix eq "");

open(COMMANDS, "<", $commandsFile) || die "Cannot open commands file $commandsFile";
while(<COMMANDS>) {
    $_ =~ s/[\r\n]+$//;
    if ($_ ne "" && !m/^[!#]/ && !m/^[ \t]+$/) {
       push(@commands, $_);
    }
}
close(COMMANDS) || die "Cannot close commands file $commandsFile";

my $envCommand = "";

if ($envFile ne "") {
    open(ENVFILE, "<", $envFile) || die "Cannot open environment file $envFile";
    while(<ENVFILE>) {
	$_ =~ s/[\r\n]+$//;
	if ($_ ne "" && !m/^[!#]/) {
	    die "Cannot parse environment line $_\n" unless m/^(\S+)=(.*)$/;
            push(@envNames, $1);
	    push(@envValues, $2);
        }
    }
    close(ENVFILE) || die "Cannot close environment file $envFile";
    $envCommand = join("; ", map("$envNames[$_]=$envValues[$_]",0..$#envNames)) . ";";
}


    my $cdPrefix = ($cd eq "") ? "" : "$cd/";
    my $cdCmd = ($cd eq "") ? "" : "cd $cd;";

    if ($verbose) {
	print "HOSTS=",join(",",@hosts),"\n";
	print "SHELL=$shell\n" if defined $hosts;
	print "Wait=$wait\nJobs=$n\nCD=$cd\nmax_load=$maxload\n";
	print "ENV: $envCommand\n";
	print "CMDS:\n",join("\n",@commands),"\n";
    }
    
    my $ncmd = 0;
    my %jobids = ();
    foreach my $host (@hosts) {
	my @list = map(-1, 0..($n-1));
        $jobids{$host} = \@list;
    }
    my @cmdJobs = map(-1, 0..$#commands);
    my %pidToNCmd = ();

    my $failed = 0;

    while ($ncmd < scalar @commands) {
	foreach my $host (@hosts) {
	    my @joblist = @{ $jobids{$host} };
	    my $shellString = $host eq "" ? "" : "$shell $host ";
	    foreach my $slot (0..($n-1)) {
		# See if job in this slot is done yet
		if ($joblist[$slot] != -1) {
		    my $child = waitpid($joblist[$slot], WNOHANG);
		    if ($child == $joblist[$slot]) {
			my $ncmdOld = $pidToNCmd{$child};
			if ($? != 0) {
			    print "Failed child $child job $ncmdOld status $? : $commands[$ncmdOld]\n";
			    $failed++;
			} else {
			    print "Completed child $child job $ncmdOld : $commands[$ncmdOld]\n" if $verbose > 0;
			}
			$joblist[$slot] = -1;
		    }
		}
		if ($joblist[$slot] == -1 && $ncmd < scalar @commands) {
		    my $upstring = `$shellString uptime`;
		    if (!$upstring) {
			die "Cannot run uptime on $host\n" if (scalar @hosts < 2);
			#else
			print "Cannot run uptime on $host\n";
		    }
		    my $up;
		    if ($upstring) {
			my @uparray = split(' ',$upstring);
			$up = $uparray[$#uparray-2];
			$up =~ s/,//;
			print "host $host slot $slot up $up ncmd $ncmd\n" if $verbose;
		    }
		    if ($upstring && $up < $maxload) {
			my $logstart = "$cdPrefix$logPrefix";
			$logstart = $logPrefix if $cdPrefix =~ m!/$! && $logPrefix =~ m!^/!;
			my $cmd = "($cdCmd $envCommand time $commands[$ncmd]) >& $logstart-$ncmd.log";
			my $execute = $host eq "" ? $cmd : "$shellString '$cmd'";
			if ($debug) {
			    print "Would execute $execute\n";
			} else {
			    my $child = fork();
			    if ($child == 0) {
				# we are the child;
				system($execute);
				print "job $ncmd FAILED $? : $commands[$ncmd]\n" if $? != 0;
				exit(0);
			    } else {
				$joblist[$slot] = $child;
				$cmdJobs[$ncmd] = $child;
				$pidToNCmd{$child} = $ncmd;
				print "Spawned $child for job $ncmd : $execute\n" if $verbose > 0;
			    }
			}
			$ncmd++;
		    } else {
			sleep($wait);
		    }
		}
	    }
	    $jobids{$host} = \@joblist;
	}
    }

    # And now wait for all the children to finish
    foreach my $host (@hosts) {
	my @joblist = @{ $jobids{$host} };
	foreach my $slot (0..($n-1)) {
	    if ($joblist[$slot] != -1) {
		my $ncmdOld = $pidToNCmd{$joblist[$slot]};
		my $child = waitpid($joblist[$slot], 0);
		if ($child ne $joblist[$slot]) {
		    print "Ack!!! waitpid for child $joblist[$slot] failed\n";
		} elsif ($? != 0) {
		    print "Failed child $child job $ncmdOld status $? : $commands[$ncmdOld]\n";
		    $failed++;
		} else {
		    print "Completed child $child job $ncmdOld : $commands[$ncmdOld]\n" if $verbose > 0;
		}
	    }
	}
    }
    if (!$debug) {
	if ($failed == 0) {
	    print "Done: all " . scalar(@commands) . " jobs completed successfully\n";
	    
	} else {
	    print "Done: $failed of " . scalar(@commands) . " jobs failed\n";
	}
	sleep(1);
    }
    exit($failed == 0 ? 0 : 1);

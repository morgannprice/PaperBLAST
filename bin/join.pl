#!/usr/bin/perl -w
# Join unsorted files, using lots of memory
#
# Example:
#
# ./join.pl -match Gene1ID=Gene1ID,Gene2ID=Gene2ID -keep 1.Group,1.Gene1ID,1.Gene1Name,1.Gene2ID,1.Gene2Name -rename 12.FisherLarge=flGitton&flOurs,12.LinearR=rGitton&rOurs /data/Gitton2002/clusters.Mm.comparisons /data/mouse-unigene-nov02/Gitton2002/clusters.Mm.comparisons

my $usage = "Usage: join.pl [-case false] [-header false] [-match [1.]Field1X=[2.]Field2X,...]\n"
    . "   [-only 10 | 120 | 1 | 2 | 12] [-uniq true ] [-rename [1|2|12.]X=Y,...] \n"
    . "   [ -ignore [1.]X,... | -keep [1|2.]X,...] | ] [-scan 1] file1 | - file2 | -\n"
    . "Joins two tab-separated tables to STDOUT. -rename implies -keep.\n";

# Field1X can be an integer (base 1) or a name, ditto for 1.X, 1.Y, etc.

my $keptOption = 0;
my $only = -1; # -1 means JOIN, 10 means LEFT JOIN, 1 means 1 not 2, 2 means 2 not 1, 12 is exclusive or, 120 is OUTER JOIN
my $header = 1;
my @matchNames1 = ();
my @matchNames2 = ();
my %rename = ();
my @ignore = ();
my @keep = ();
my $file1 = "";
my $file2 = "";
my @names1 = ();
my @names2 = ();
my @data1;
my @data2;
my @iMatch1 = ();
my @iMatch2 = ();

my @bOut1 = ();			# true if output field is from 1
my @iOut = ();			# index of each output field
my @bShow = ();

my @newnames = ();

my @empty = ();
my @emptyempty = ( \@empty );
my $case = 1; # default is case sensitive
my $uniq = 0; # default is no uniqueness requirement

my %index1 = ();
my %index2 = ();

my $scan = 0; # scan first file instead of loading into memory

sub myuc($) {
    return $_[0] if $case;
    return uc($_[0]);
}

sub addRename($$) {
    my ($from, $to) = @_;
    die "Repeated field to rename" if exists $rename{$from};
    $rename{$from} = $to;
    push(@keep, $from);
}

while(scalar @ARGV > 2) {
    my $option = shift @ARGV;
    my $value = shift @ARGV;
    if ($option eq "-header") {
	$header = ($value eq "true" || $value eq "TRUE" || $value eq "yes" || $value eq "YES" || $value eq "1");
    } elsif ($option eq "-match") {
	foreach (split(/,/,$value)) {
	    my @matches = split(/=/, $_);
	    die $usage if scalar @matches != 2;
	    push(@matchNames1, $matches[0]);
	    push(@matchNames2, $matches[1]);
	}
    } elsif ($option eq "-rename") {
	my @renames = split(/,/,$value);
	foreach (@renames) {
	    my @fields = split(/=/,$_,-1);	# no limit
	    die $usage if scalar @fields != 2 || $fields[0] eq "";
	    my $from = $fields[0];
	    my $to = $fields[1];
	    if ($from =~ m/^12\.(.*)$/) {
		$fromBoth = $1;
		my @tos = split(/\+/, $to);
		die "Cannot parse dual rename $to" if scalar @tos != 2;
		&addRename("1." . $fromBoth, $tos[0]);
		&addRename("2." . $fromBoth, $tos[1]);
	    } else {
		&addRename($from,$to);
	    }
	}
    } elsif ($option eq "-keep") {
	push(@keep, split(/,/,$value));
	$keptOption = 1;
    } elsif ($option eq "-ignore") {
	push(@ignore, split(/,/,$value) );
    } elsif ($option eq "-only") {
	die $usage if $value != 1 && $value != 2 && $value != 12 && $value != 10 && $value != 120;
	$only = $value;
    } elsif ($option eq "-case") {
	$case = ($value eq "true" || $value eq "TRUE" || $value eq "yes" || $value eq "YES" || $value eq "1");
    } elsif ($option eq "-uniq") {
	$uniq = ($value eq "true" || $value eq "TRUE" || $value eq "yes" || $value eq "YES" || $value eq "1");
    } elsif ($option eq "-scan") {
	$scan = ($value eq "true" || $value eq "TRUE" || $value eq "yes" || $value eq "YES" || $value eq "1");
    } else {
	die "Cannot recognize option $option:\n$usage";
    }
}

# We can add to keep via -rename, but that is ignored if there was no actual keep list
$#keep = -1 unless $keptOption;

die $usage if scalar @keep > 0 && scalar @ignore > 0;
die $usage if !$header && scalar (keys %rename) > 0;

die "Cannot use -scan 1 and -header 0" if $scan && !$header;
die "Cannot use -scan 1 and -only" if $scan && $only ne "-1";

die $usage if scalar @ARGV != 2;
$file1 = shift @ARGV;
$file2 = shift @ARGV;
die $usage if $file1 eq $file2 && $file1 eq "-";

sub load($$$$) {
    my ($file, $nameref, $dataref,$scan) = @_;
    if ($file eq "-") {
	open(FILE, "<&STDIN") || die "Cannot dup STDIN";
    } else {
	open(FILE, "<", $file) || die "Cannot read $file";
    }
    my $first = 1;
    while(<FILE>) {
	$_ =~ s/[\r\n]+$//;
	my @F = split(/\t/, $_, -1);
	if ($first) {
	    push( @{ $nameref }, @F ); # save first line as names even if no header for length checking
	}
	if ($header ? !$first : 1) {
	    push( @{ $dataref }, \@F );	# first line is data if no header
	}
	die "Too few columns in $file" if (!$first && scalar @F < scalar @{ $nameref });
	$first = 0;
	return if $scan;
    }
    close(FILE) || die "Cannot close $file";
}

sub getColumn($$$) {
    my ($name, $prefix, $namesref) = @_;
    if ($name =~ m/^(.*)\.(.*)$/) {
	return -1 if $1 ne $prefix;
	$name = $2;
    }
    return $name - 1 if ($name =~ m/^[0-9]+/ && $name > 0 && $name <= scalar @{$namesref} );
    return -1 if !$header;
    my @names = @{$namesref};
    map { return $_ if myuc($name) eq myuc($names[$_]) } (0..$#names);
    return -1;
}

sub joinFields($$) {
    my ($fieldsref, $indexfieldsref) = @_;
    return join("\t", map $fieldsref->[$_], @{ $indexfieldsref } );
}

sub addToIndex($$$) {
    my ($indexref,$fieldsref,$indexfieldsref) = @_;
    my $match = &joinFields($fieldsref, $indexfieldsref);
    my @hits = ();
    @hits = @{ $indexref->{myuc($match)} } if exists $indexref->{myuc($match)};
    push(@hits, $fieldsref);
    $indexref->{myuc($match)} = \@hits;
}

sub printMatches($$) {
    my ($list1, $list2) = @_;
    my @F1s = @{ $list1 };
    my @F2s = @{ $list2 };
    foreach (@F1s) {
	my @F1 = @{ $_ };
	my $nDups = ($uniq ? (scalar @F1s)*(scalar @F2s) : 1);
	foreach (@F2s) {
	    my @F2 = @{ $_ };
	    my $first = 1;
	    # Finally we get to print
	    foreach (0..$#bShow) {
		if ($bShow[$_]) {
		    my $b1 = $bOut1[$_];
		    my $i = $iOut[$_];
		    my $field;
		    if ($b1 && scalar @F1 == 0) {
			$field = "";
		    } elsif (!$b1 && scalar @F2 == 0) {
			$field = "";
		    } elsif ($b1) {
			$field = $F1[$i];
		    } else {
			$field = $F2[$i];
		    }
		    if ($nDups > 1) {
			print STDERR "\t" if (!$first);
			print STDERR $field;
		    } else {
			print "\t" if (!$first);
			print $field;
		    }
		    $first = 0;
		}
	    }
	    if ($nDups > 1) { print STDERR "\t$nDups\n"; }  else { print "\n"; }
	}
    }
}

&load($file2,\@names2,\@data2,0);
&load($file1,\@names1,\@data1,$scan); # FILE still open if $scan

@iMatch1 = map { my $i = &getColumn($_, 1, \@names1);
		 die "Cannot find match name for $_" if $i < 0;
		 $i;
	     } @matchNames1;
push(@iMatch1, 0) if scalar @iMatch1 == 0;

@iMatch2 = map { my $i = &getColumn($_, 2, \@names2);
		 die "Cannot find match name for $_" if $i < 0;
		 $i;
	     } @matchNames2;
push(@iMatch2, 0) if scalar @iMatch2 == 0;


if (scalar @keep > 0) {
    foreach my $out (@keep) {
	my $i = &getColumn($out, 1, \@names1);
	my $b1 = 1;
	if ($i == -1) {
	    $i = &getColumn($out, 2, \@names2);
	    $b1 = 0;
	}
	die "Cannot recognize keep $out" if $i < 0;
	push(@bOut1, $b1);
	push(@iOut, $i);
	push(@bShow,1);
    }
} else {
    map { push(@bOut1, 1); push(@iOut, $_); } (0..$#names1);
    map { push(@bOut1, 0); push(@iOut, $_) } (0..$#names2);
    map { push(@bShow, 1); } (0..$#iOut);

    foreach my $ign (@ignore) {
	my $i = &getColumn($ign, 1, \@names1);
	my $b1 = 1;
	if ($i == -1) {
	    $i = &getColumn($ign, 2, \@names2);
	    $b1 = 0;
	}
	die "Cannot recognize ignore $ign" if $i < 0;
	$i += scalar @names1 if !$b1;
	$bShow[$i] = 0;
    }
}

if ($header) {
    foreach (0..$#iOut) {
	my $i = $iOut[$_];
	my $b1 = $bOut1[$_];
	my $oldname = $b1 ? $names1[$i] : $names2[$i];
	my $newname = $oldname;
	my $prefix = $b1 ? "1" : "2";
	if (exists $rename{$oldname}) {
	    $newname = $rename{$oldname};
	} elsif (exists $rename{"$prefix.$oldname"}) {
	    $newname = $rename{"$prefix.$oldname"};
	}
	push(@newnames, $newname);
    }
}

# Build an index on the match field field for 1 and 2

if ($only == -1 || $only == 2 || $only == 12 || $only == 120 ) {
    map { addToIndex(\%index1, $_, \@iMatch1); } @data1;
}
if ($only == -1 || $only == 1 || $only == 12 || $only == 10 || $only == 120) {
    map { addToIndex(\%index2, $_, \@iMatch2); } @data2;
}

# Output

if ($header) {
    my $first = 1;
    foreach (0..$#newnames) {
	if ($bShow[$_]) {
	    print "\t" if (!$first);
	    print $newnames[$_];
	    if ($uniq) {
		print STDERR "\t" if (!$first);
		print STDERR $newnames[$_];
	    }
	    $first = 0;
	}
    }
    print "\n";
    print STDERR "\tnDuplicates\n" if $uniq;
}

if ($only == -1) {
    if ($scan) {
	while(<FILE>) {
	    $_ =~ s/[\r\n]+$//;
	    my @F = split(/\t/, $_, -1);
	    die "Too few columns in $file1" if scalar @F < scalar @names1;
	    my @F1s = ( \@F );
	    my $match = joinFields(\@F, \@iMatch1);
	    &printMatches(\@F1s, $index2{myuc($match)}) if exists $index2{myuc($match)};
	}
	close (FILE) || die "Cannot close $file1";
    } else {
	foreach (@data1) {
	    my @F1 = @{ $_ };
	    my @F1s = ( \@F1 );
	    my $match = joinFields(\@F1, \@iMatch1);
	    &printMatches(\@F1s, $index2{myuc($match)}) if exists $index2{myuc($match)};
	}
    }
}

if($only == 1 || $only == 12) {
    foreach (@data1) {
	my @F1 = @{ $_ };
	my @F1s = ( \@F1 );
	my $match = joinFields(\@F1, \@iMatch1);
	&printMatches(\@F1s, \@emptyempty) if !exists $index2{myuc($match)};
    }
}

if ($only == 10 || $only == 120) {
    foreach (@data1) {
	my @F1 = @{ $_ };
	my @F1s = ( \@F1 );
	my $match = joinFields(\@F1, \@iMatch1);
	&printMatches(\@F1s, (exists $index2{myuc($match)} ? $index2{myuc($match)} : \@emptyempty));
    }
}

if ($only == 2 || $only == 12 || $only == 120) {
    foreach (@data2) {
	my @F2 = @{ $_ };
	my @F2s = ( \@F2 );
	my $match = joinFields(\@F2, \@iMatch2);
	&printMatches(\@emptyempty, \@F2s) if !exists $index1{myuc($match)};
    }
}


#!/usr/bin/perl -w
# Modified from
# download_all_sets.pl
# from the BioLip web site
use strict;

chdir "BioLiP_updated_set" || die "No such subdirectory: BioLiP_updated_set";

my $head= "http://zhanglab.ccmb.med.umich.edu/BioLiP/weekly";
my $address="http://zhanglab.ccmb.med.umich.edu/BioLiP/weekly.html";
unlink("log");
unlink("weekly.html");
system("wget -o log -c $address") == 0 or die "System call failed: $!";


my @rst=`cat weekly.html`;
my @all=();
foreach my $r(@rst)
{
    if($r =~ /\<tr\>\<td\>(\S+)\<\/td\>/)
    {
	#print "$1\n";
	push(@all, $1);
    }
}
my $tot=@all;
print "\n====================================================\n";
print "In total, there are $tot weeks to update.\n\n";

my $annotation="BioLiP_UP.txt";
open(OUT, ">$annotation") || die;
close(OUT) || die;
my $annotation1="BioLiP_UP_nr.txt";
open(OUT, ">$annotation1") || die;
close(OUT) || die;

foreach my $r(@all)
{        
#    my $rec="receptor_$r.tar.bz2";
#    my $rec1="receptor1_$r.tar.bz2";
#    my $lig="ligand_$r.tar.bz2";
    my $ano="BioLiP_$r.txt";
    my $ano1="BioLiP_$r\_nr.txt";

    if(-s $ano && -s $ano1)
    {
	print "The week $r was updated before, skip this one\n";	
    }
    else
    {
	print "Dowload redundant set for the week $r...\n";    
#	system("wget -o log -c $head/$rec") == 0 or die "System call failed: $!";
#	system("tar -xvf $rec >log")== 0 or die "System call failed: $!";
#	system("wget -o log -c $head/$rec1") == 0 or die "System call failed: $!";
#	system("tar -xvf $rec1>log")== 0 or die "System call failed: $!";
#	system("wget -o log -c $head/$lig") == 0 or die "System call failed: $!";
#	system("tar -xvf $lig >log")== 0 or die "System call failed: $!";
	system("wget -o log -c $head/$ano") == 0 or die "System call failed: $!"; 
	
	
	print "Dowload non-redundant set for the week $r...\n";
#	$rec="receptor_$r\_nr.tar.bz2";
#	$rec1="receptor1_$r\_nr.tar.bz2";
#	$lig="ligand_$r\_nr.tar.bz2";
#	system("wget -o log -c $head/$rec") == 0 or die "System call failed: $!";
#	system("tar -xvf $rec >log")== 0 or die "System call failed: $!";
#	system("wget -o log -c $head/$rec1") == 0 or die "System call failed: $!";
#	system("tar -xvf $rec1>log")== 0 or die "System call failed: $!";
#	system("wget -o log -c $head/$lig") == 0 or die "System call failed: $!";
#	system("tar -xvf $lig >log")== 0 or die "System call failed: $!";
	system("wget -o log -c $head/$ano1") == 0 or die "System call failed: $!"; 
    }

    system("cat $ano >> $annotation") == 0 || die "cat failed: $!";
    system("cat $ano1 >> $annotation1") == 0 || die "cat failed: $!";

}

print "Cheers! All updates are done.\n";
print "====================================================\n\n";

print "Please download the old sets manually at 
http://zhanglab.ccmb.med.umich.edu/BioLiP/download.html

Please read http://zhanglab.ccmb.med.umich.edu/BioLiP/download/readme.txt 
about explanation of the annotation file.
";

print "Please feel free to contact me (Jianyi, yangji\@umich.edu) if you have any problems with BioLiP.\n
Thanks for using the BioLiP database!

---------------------------------------------
Please cite the following paper if you use BioLiP in your projects:

Jianyi Yang, Ambrish Roy, Yang Zhang, BioLiP: a semi-manually curated database for biologically 
relevant ligand-protein interactions, Nucleic Acids Research, 41:D1096-D1103, 2013. 

";


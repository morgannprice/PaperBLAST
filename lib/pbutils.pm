# Utilities for PaperBLAST
package pbutils;
require Exporter;
use strict;

our (@ISA,@EXPORT);
@ISA = qw(Exporter);
@EXPORT = qw(read_list wget ftp_html_to_files write_list mkdir_if_needed);

sub read_list($) {
  my ($file) = @_;
  open(LIST, "<", $file) || die "Cannot read $file";
    my @lines = <LIST>;
    close(LIST) || die "Error reading $file";
    @lines = map { chomp; $_ } @lines;
    return(@lines);
}

# Retry wget calls to increase odds of success
my $nFailures = 0;
my $nFailMax = 5;

sub wget($$) {
    my ($url, $file) = @_;
    while(system("wget", "-nv", "-O", $file, $url) != 0) {
      $nFailures++;
      die "Failed to load $url\n" if $nFailures >= $nFailMax;
      print STDERR "Waiting and retrying\n";
      sleep(1);
    }
}

sub ftp_html_to_files($) {
    my ($listfile) = @_;
    my @files = ();
    open(IN, "<", $listfile) || die;
    while(<IN>) {
        chomp;
        next unless m!>([A-Za-z0-9_.-]+)<!;
        my $file = $1;
        push @files, $file;
    }
    close(IN) || die "Error reading $listfile";
    return @files;
}

sub write_list($$) {
    my ($files, $outfile) = @_;
    open(LIST, ">", $outfile) || die "Cannot write to $outfile";
    foreach my $file (@$files) {
        print LIST "$file\n";
    }
    close(LIST) || die "Error writing to $outfile";
    print STDERR "Wrote $outfile\n";

}

sub mkdir_if_needed($) {
    my ($subdir) = @_;
    (-d $subdir) || mkdir($subdir) || die "Cannot mkdir $subdir";
}

1;

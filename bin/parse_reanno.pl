#!/usr/bin/perl -w
# Parse reannotations from the Fitness Browser
use strict;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use pbutils;

my $usage = <<END
parse_reanno.pl anno.tab aaseqs orginfo > reanno.curated_parsed

The annotation table must be tab-delimited and must include the fields
orgId, locusId, sysName, new_annotation, and comment.

The aaseqs file must be in fasta format with sequence identifiers of
the form orgId:locusId.

The organism information table must include the fields name or orgId,
genus, species, and strain.

END
;

die $usage unless @ARGV == 3;
my ($annofile, $seqsfile, $orgfile) = @ARGV;
foreach my $file ($annofile,$seqsfile,$orgfile) {
  die "No such file: $file\n" unless -e $file;
}

my @anno = ReadTable($annofile, qw(orgId locusId sysName new_annotation comment));
print STDERR "Read " . scalar(@anno) . " reannotations\n";
my %anno = map { $_->{orgId} . ":" . $_->{locusId} => $_ } @anno;

my @orgs = ReadTable($orgfile, qw(genus species strain));
die "No organisms\n" unless @orgs > 0;
if (exists $orgs[0]{name}) {
  foreach my $org (@orgs) {
    $org->{orgId} = $org->{name};
  }
}
die "organism file must include orgId or name field\n"
  unless exists $orgs[0]{orgId};
my %orgs = map { $_->{orgId} => $_ } @orgs;

foreach my $anno (@anno) {
  die "Unknown organism id $anno->{orgId}"
    unless exists $orgs{$anno->{orgId}};
}

open(my $fh, "<", $seqsfile) || die "Cannot read $seqsfile";
my $state = {};
while(my ($header, $seq) = ReadFastaEntry($fh, $state)) {
  $header =~ s/ .*//; # ignore any description
  $anno{$header}{seq} = $seq if exists $anno{$header};
}
close($fh) || die "Error reading $seqsfile";

foreach my $anno (@anno) {
  if (!exists $anno->{seq}) {
    print STDERR "Warning, no sequence for $anno->{orgId}:$anno->{locusId}\n";
  } else {
    my $org = $orgs{ $anno->{orgId} };
    die unless defined $org;
    my $comment = $anno->{comment};
    $comment =~ s/SEED_correct//;
    $comment =~ s/KEGG_correct//;
    $comment =~ s/[(][)]$//;
    $comment =~ s/ $//;
    $comment =~ s/;$//;
    print join("\t", "reanno", $anno->{orgId} . ":" . $anno->{locusId},
               $anno->{sysName} || $anno->{locusId},
               "", # short name
               $anno->{new_annotation},
               join(" ", $org->{genus}, $org->{species}, $org->{strain}),
               $anno->{seq},
               $comment,
               "" # pubmed ids
              )."\n";
  }
}

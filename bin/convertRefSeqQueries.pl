#!/usr/bin/perl -w
use strict;

my $usage = <<END
Run as a filter.

Input should be tab-delimited with a variable number of fields. The
first field has the query term (usually the locus_tag, the
old_locus_tag, or the protein id without the version number, but it
could also be the gene or gene_synonym). The other fields are of the
form key=value.

The output is suitable for input to queryEuropePMC.pl and includes the
organism, the query term, the locus_tag or protein_id, sequence, and description or "product".
END
    ;

die $usage if @ARGV > 0;

my %queries = (); # query => list of hash of field values

my $nRead = 0;
while(my $line = <STDIN>) {
    chomp $line;
    my @F = split /\t/, $line;
    my $query = shift @F;
    my %attr = ();
    foreach my $f (@F) {
        $f =~ m/^([^=]+)=(.*)$/ || die "Cannot parse key=value from $f";
        $attr{$1} = $2;
    }
    $nRead++;
    my @missing = ();
    foreach my $field (qw{organism protein_id translation}) {
        push @missing, $field if !exists $attr{$field};
    }
    if (scalar(@missing) > 0) {
        print STDERR "Skipping -- missing " . join(",",@missing) . " :\t" . $line . "\n"
            unless exists $attr{exception};
        next;
    }
    push @{ $queries{$query} }, \%attr;
}
print STDERR "Read $nRead queries\n";

my $nOut = 0;
while (my ($query, $list) = each %queries) {
    # If there are multiple entries, do they all belong to the same genus and have the same protein id?
    # (This can be legitimate if there are multiple isoforms, but it is not supported)
    my @prot = map $_->{protein_id}, @$list;
    my %prot = map { $_ => 1 } @prot;
    if (scalar(keys %prot) > 1) {
        print STDERR "Skipping non-unique term $query for multiple protein ids\n";
        next;
    }
    my @genus = map $_->{organism}, @$list;
    @genus = map s/ .*//, @genus;
    my %genus = map { $_ => 1 } @genus;
    if (scalar(keys %genus) > 1) {
        print STDERR "Skipping non-unique term $query in mutiple genera " . join(" ", sort keys %genus) . "\n";
        next;
    }
    # Arbitrarily use the first entry
    my $attr = $list->[0];
    print join("\t",
               $attr->{organism},
               $query,
               $attr->{protein_id},
               $attr->{translation},
               $attr->{product} || "" ) . "\n";
    $nOut++;
}
print STDERR "Wrote $nOut queries\n";

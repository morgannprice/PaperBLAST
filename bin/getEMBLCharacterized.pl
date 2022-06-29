#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;

die <<END
Usage: getEMBLCharacterized.pl < ENA_entry.embl > experiment.tab

For entries that have at least one pubmed identifier, saves any CDS
features that have the /experiment tag. The output is tab-delimited
with no header line and with fields entry id, organism, pubmed ids
(comma separated), protein id, gene name, description, location in
genbank format (such as "<1..267" or "complement(34109..34717)"), and
amino acid sequence.

Example output line:
L17342.1	Haemophilus parainfluenzae	2183189,7514149	AAA20481.1	hpaIIM	HpaII modification methyltransferase	1594..2670	MKDVLDDNLLEEPAAQYSLFEPESNPNLREKFTFIDLFAGIGGFRIAMQNLGGKCIFSSEWDEQAQKTYEANFGDLPYGDITLEETKAFIPEKFDILCAGFPCQAFSIAGKRGGFEDTRGTLFFDVAEIIRRHQPKAFFLENVKGLKNHDKGRTLKTILNVLREDLGYFVPEPAIVNAKNFGVPQNRERIYIVGFHKSTGVNSFSYPEPLDKIVTFADIREEKTVPTKYYLSTQYIDTLRKHKERHESKGNGFGYEIIPDDGIANAIVVGGMGRERNLVIDHRITDFTPTTNIKGEVNREGIRKMTPREWARLQGFPDSYVIPVSDASAYKQFGNSVAVPAIQATGKKILEKLGNLYD
END
  unless @ARGV == 0;


# get 1st value for tag, or empty string
sub getTag($$);

my $seqIO = Bio::SeqIO->new(-fh => \*STDIN, -format => 'EMBL');
my $nSeq = 0;
while (my $seq = $seqIO->next_seq) {
  $nSeq++;
  print STDERR "Processed $nSeq entries\n" if $nSeq > 0 && ($nSeq % 20000) == 0;
  my @pmIds;
  foreach my $ref ($seq->annotation->get_Annotations('reference')) {
    push @pmIds, $ref->medline if $ref->medline;
    push @pmIds, $ref->pubmed if $ref->pubmed;
  }
  next unless @pmIds > 0;
  my $id = $seq->display_id;
  my $version = $seq->seq_version;
  my $species = $seq->species->node_name;
  foreach my $cds ($seq->get_SeqFeatures('CDS')) {
    next unless $cds->has_tag('experiment') && $cds->has_tag('translation');
    my $gene = getTag($cds, 'gene');
    my $desc = getTag($cds, 'product');
    my $protId = getTag($cds, 'protein_id');
    my $aa = getTag($cds, 'translation');
    next unless $aa;
    print join("\t", $id.".".$version, $species, join(",",@pmIds),
               $protId, $gene, $desc, $cds->location->to_FTstring, $aa)."\n";
  }
}
print STDERR "Finished processing $nSeq entries\n";

sub getTag($$) {
  my ($ft, $tag) = @_;
  return "" unless $ft->has_tag($tag);
  my @values = $ft->get_tag_values($tag);
  return $values[0];
}

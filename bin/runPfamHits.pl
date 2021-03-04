#
# runPfamHits.pl -- build the pfam.hits.tab file for a set, in tmp/path.set

[[ $# -eq 0 ]] && echo "Usage: runPfamHits.pl set" && exit 0
set=$1
[[ ! -d tmp/path.$set ]] && echo Directory tmp/path.$set does not exist! && exit 1

for file in hmm/pfam.tab bin/hmmsearch bin/hmmfetch bin/submitter.pl hmm/Pfam-A.hmm; do
[[ ! -f $file ]] && echo file $file does not exist! && exit 1
done

echo Computing tmp/path.$set/pfam.hits.tab
mkdir /tmp/pfam.hits.$$
tail -n +2 hmm/pfam.tab | cut -f 1 > /tmp/pfam.hits.$$/list
for i in `cat /tmp/pfam.hits.$$/list`; do bin/hmmfetch hmm/Pfam-A.hmm $i > /tmp/pfam.hits.$$/$i.hmm; done
echo Fetched hmms
(for i in `cat /tmp/pfam.hits.$$/list`; do echo "bin/hmmsearch --cut_tc --domtblout /tmp/pfam.hits.$$/$i.domtbl --cpu 1 -o /dev/null /tmp/pfam.hits.$$/$i.hmm tmp/path.$set/curated.faa"; done) > /tmp/pfam.hits.$$/cmds
bin/submitter.pl /tmp/pfam.hits.$$/cmds >& /tmp/pfam.hits.$$/cmds.log
(for i in `cat /tmp/pfam.hits.$$/list`; do perl -ane 'next if m/^[#-]/; ($ids, $dash1, $qlen, $hmmName, $hmmAcc, $hmmLen, $seqEval, $seqBits, $seqBias, $domI, $domN, $evalC, $evalE, $bits, $bias, $hmmFrom, $hmmTo, $seqFrom, $seqTo) = split / +/; die $_ unless $dash1 eq "-" && $ids =~ m/:/; print join("\t", $ids, $hmmName, $hmmAcc, $evalC, $bits, $seqFrom, $seqTo, $qlen, $hmmFrom, $hmmTo, $hmmLen)."\n";' < /tmp/pfam.hits.$$/$i.domtbl; done) | sort > tmp/path.$set/pfam.hits.tab
echo Wrote tmp/path.$set/pfam.hits.tab
rm -Rf /tmp/pfam.hits.$$

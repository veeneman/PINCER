#!/bin/bash
#code to generate genome-wide provean-based conservation scores
#not yet tested for reproducibility
BIN=`dirname "$0"`;
REF="/hpc/grid/oncology/veeneb/ref";

SPECIES="hg38" #"mm10"

cd /home/veeneb/hpc/ref/$SPECIES/anno; #mm10 too
mkdir provean;
cd provean;

cut -f 3 ../appris.txt | tail -n +2 | sort > NP.txt
$BIN/lib/split_fasta.pl $REF/$SPECIES/anno/refseq.aa.gz NP.txt;
head -n 5000 NP.txt > NP.1.txt
tail -n +5001 NP.txt | head -n 5000 > NP.2.txt
tail -n +10001 NP.txt | head -n 5000 > NP.3.txt
tail -n +15001 NP.txt > NP.4.txt
rm NP.txt

bsub -q general -app large -n 64 -W 16:00 -o "1.out" -e "1.err" "set -e; cat NP.1.txt | xargs -n 1 -P 8 $BIN/lib/run_provean.pl";
bsub -q general -app large -n 64 -W 16:00 -o "2.out" -e "2.err" "set -e; cat NP.2.txt | xargs -n 1 -P 8 $BIN/lib/run_provean.pl";
bsub -q general -app large -n 64 -W 16:00 -o "3.out" -e "3.err" "set -e; cat NP.3.txt | xargs -n 1 -P 8 $BIN/lib/run_provean.pl";
bsub -q general -app large -n 64 -W 16:00 -o "4.out" -e "4.err" "set -e; cat NP.4.txt | xargs -n 1 -P 8 $BIN/lib/run_provean.pl";

cat 1.out 2.out 3.out 4.out > bsub.out
rm 1.out 2.out 3.out 4.out
cat 1.err 2.err 3.err 4.err > bsub.err
rm 1.err 2.err 3.err 4.err

for i in `cat NP.1.txt NP.2.txt NP.3.txt NP.4.txt`; do cat $i.log; done > provean.log
cat NP.1.txt NP.2.txt NP.3.txt NP.4.txt | \
 perl -ne 'chomp; print "$_.log2\n";' | xargs grep . | \
 grep -v "Selenocysteine" > provean.log2

for i in `cat NP.1.txt NP.2.txt NP.3.txt NP.4.txt`; do cat $i.bed; done > provean.aa.bed
for i in `cat NP.1.txt NP.2.txt NP.3.txt NP.4.txt`; do rm $i.bed $i.log $i.log2; done
rm NP.{1..4}.txt

#look at files, then
rm provean.log provean.log2 bsub.out bsub.err

#convert aa.bed to nt.bedgraph
$BIN/lib/process.conservation.R provean.aa.bed provean.bedgraph \
 $REF/$SPECIES/anno/cdd.aa.gff2.gz $REF/$SPECIES/anno/refseq.to.genome.gff3.gz

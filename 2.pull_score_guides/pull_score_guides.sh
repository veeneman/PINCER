#!/bin/bash
#Note - this code has Not been tested for reproducibility beyond its original run
# For bugs or questions, contact brendan veeneman.
#REQUIREMENTS:
# Currently implemented LSF queue (bsub)
# In path: R; samtools; perl; python 2.7;
#          batmis, batman, batdecode, etc. from Batmis
set -e;

#reference filepath, from previous step (1.download_annotations)
REF="/hpc/grid/oncology/veeneb/ref";
G="hg38"; #"mm10"; #

BIN=`dirname "$0"`;
BSUB_OPTS="-q general -app medium -n 1";
BSUB_4g="-R \"rusage[mem=4194304]\" -M 4194304";
BSUB_6g="-R \"rusage[mem=6291456]\" -M 6291456";
BSUB_16g="-R \"rusage[mem=16777216]\" -M 16777216";

FA="$REF/$G/batmis/$G.fa";
CDS="$REF/$G/anno/gencode.CDS.Rdata";
OT="5000";
C="2000"; #chunks
INPUT=$FA;

mkdir log;
bsub $BSUB_OPTS -W 4:00 $BSUB_4g -J $G.1 -o log/$G.01.o -e log/$G.01.e \
 "$BIN/lib/pull_guides.pl $INPUT $G.base.sam $G.30nt $G.fa";

mkdir tmp;
cd tmp;
bsub $BSUB_OPTS -W 4:00 $BSUB_4g -J $G.2 -o ../log/$G.02.o -e ../log/$G.02.e -w "done($G.1)" \
 "set -e; $BIN/lib/split_file.pl ../$G.30nt .30nt $C 1; rm ../$G.30nt;";
bsub $BSUB_OPTS -W 4:00 $BSUB_4g -J $G.3 -o ../log/$G.03.o -e ../log/$G.03.e -w "done($G.1)" \
 "set -e; $BIN/lib/split_file.pl ../$G.fa .fa $C 2; rm ../$G.fa;";

bsub $BSUB_OPTS -W 4:00 $BSUB_16g -J "$G.4[1-$C]" -o "%I.04.o" -e "%I.04.e" -w "done($G.2)" \
 "set -e; i=\$LSB_JOBINDEX;\
  if [ -e \$i.azm ]; then exit 0; fi;\
  $BIN/lib/score_azimuth.py \$i.30nt \$i.azm;\
  rm \$i.30nt;";

bsub $BSUB_OPTS -W 48:00 $BSUB_6g -J "$G.5[1-$C]" -o "%I.05.o" -e "%I.05.e" -w "done($G.3)" \
 "set -e; i=\$LSB_JOBINDEX;\
  if [ -e \$i.fa ]; then \
   batman -n 5 --nismismatch -m $OT -g $FA -l /dev/null --mishits=\$i.um -q \$i.fa -o \$i.bat;\
   if [ -s \$i.um ]; then >&2 echo 'UNMAPPED GUIDE'; exit 1; fi;\
   rm \$i.fa \$i.um; fi;\
  if [ -e \$i.bat ]; then \
   batdecode --misplus -m $OT -g $FA -L /dev/null -i \$i.bat -o /dev/stdout | \
    $BIN/lib/clean_sam.pl \$i.lim $OT | \
    samtools view -b - -o \$i.bam;\
   rm \$i.bat; fi;\
  if [ -e \$i.bam ]; then \
   $BIN/lib/score_hsu13.R $CDS $OT \$i.bam > \$i.hsu;\
   samtools view -C -T $FA \$i.bam -o \$i.cram;\
   rm \$i.bam; fi;\
  if [ ! -e \$i.hsu ]; then >&2 echo 'NO RESULT'; exit 1; fi;";

bsub $BSUB_OPTS -W 48:00 $BSUB_4g -J $G.6 -o ../log/$G.06.o -e ../log/$G.06.e \
 -w "done($G.4[1-$C])&&done($G.5[1-$C])" \
 "set -e;\
  for i in 04 05; do for j in o e; do \
   grep . \`eval echo \"{1..$C}.\$i.\$j\"\` | perl -ne 'print' > ../log/$G.\$i.\$j;\
   rm \`eval echo \"{1..$C}.\$i.\$j\"\`; done; done;\
  for i in azm lim hsu; do \
   cat \`eval echo \"{1..$C}.\$i\"\` > ../$G.\$i;\
   rm \`eval echo \"{1..$C}.\$i\"\`; done;\
  samtools view -H 1.cram > ../$G.samh;\
  samtools cat \`eval echo \"{1..$C}.cram\"\` > ../$G.ot.cram;\
  rm \`eval echo \"{1..$C}.cram\"\`;";

cd ..;
bsub $BSUB_OPTS -W 24:00 $BSUB_4g -J $G.7 -o log/$G.07.o -e log/$G.07.e -w "done($G.6)" \
 "set -e; rmdir tmp;\
  $BIN/lib/final_join.pl $G.samh $G.base.sam $G.azm $G.lim $G.hsu > $G.sam;\
  samtools view $G.sam -bo $G.bam;\
  rm $G.base.sam $G.samh $G.azm $G.lim $G.hsu $G.sam;";

#reposition guides onto cut sites, and coordinate sort
#XXX - THIS DIDN'T WORK, because repositioned guides were still 20nt long
#bsub -q priority -app large -W 48:00 -n 8 -R "rusage[mem=67108864]" -M 67108864 "./lib/bam2cssortbam.sh hg38.bam hg38.CS.sort.bam"
#bsub -q priority -app large -W 48:00 -n 8 -R "rusage[mem=67108864]" -M 67108864 "./lib/bam2cssortbam.sh mm10.bam mm10.CS.sort.bam"
#samtools view -L ~/hpc/ref/hg38/anno/gencode.CDS.bed hg38.CS.sort.bam -bo hg38.CSCDS.bam 
#samtools view -L ~/hpc/ref/mm10/anno/gencode.CDS.bed mm10.CS.sort.bam -bo mm10.CSCDS.bam 
#rm hg38.CS.sort.bam mm10.CS.sort.bam

#samtools view -h -L ~/hpc/ref/hg38/anno/refseq.CDS.bed hg38.bam | samtools sort -o hg38.CDS.bam
#samtools view -h -L ~/hpc/ref/mm10/anno/refseq.CDS.bed mm10.bam | samtools sort -o mm10.CDS.bam

#filter for guides touching CDS, and expand out arrayed tags (which R has issues with)
for i in hg38 mm10; do
 samtools view -h -L $REF/$i/anno/refseq.clean.uniq_cds.bed $i.bam | perl -ne \
  's/MG:B:S,(\d+),(\d+),(\d+),(\d+),(\d+)/G0:i:$1\tG1:i:$2\tG2:i:$3\tG3:i:$4\tG4:i:$5/; \
   s/MP:B:S,(\d+),(\d+),(\d+),(\d+),(\d+)/P0:i:$1\tP1:i:$2\tP2:i:$3\tP3:i:$4\tP4:i:$5/; \
   print;' | samtools sort -o $i.CDS.bam;
done


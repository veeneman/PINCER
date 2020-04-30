#!/bin/bash
#This code is reproducible as of September 4th, 2018.
#CDD, and parts of Refseq (the genome map, and AA fasta) are not
# versioned by NCBI, and are therefore difficult to manage. We
# cleaned up the version we froze on, but further versions would need
# further curation.
#REQUIREMENTS:
# In path: samtools, perl, R, bigBedToBed
# R libs: rtracklayer, ?
# readseq.jar
# batman/batmis aligner

#URLs - Firewall blocked FTP
UCSC_URL="http://hgdownload.cse.ucsc.edu";
GENCODE_URL="http://ftp.sanger.ac.uk/pub/gencode";
MAPPING_URL="https://raw.githubusercontent.com/dpryan79/ChromosomeMappings/master";
NCBI_URL="http://ftp.ncbi.nlm.nih.gov";
APPRIS_URL="http://apprisws.bioinfo.cnio.es/pub/current_release/datafiles";

#Software
READSEQJAR="/home/veeneb/afs/software/bin/readseq.jar";
BATMIS="/home/veeneb/afs/software/src/BatMis-3.00";

#Cleaning scripts
SCRIPT_BASE=`dirname "$0"`"/cleaning_scripts";
#"/home/veeneb/afs/software/guideDB/download_annotations/cleaning_scripts"
GENCODE2UCSC="$SCRIPT_BASE/gencode2ucsc.pl";
GFF2CDS="$SCRIPT_BASE/gff2cds.R";
CLEANGENCODE="$SCRIPT_BASE/clean.gencode.R";
REFSEQ2UCSC="$SCRIPT_BASE/refseq2ucsc.pl";
CLEANREFSEQ="$SCRIPT_BASE/clean.refseq.R";
CLEANGENBANK="$SCRIPT_BASE/clean.genbank.pl";
CLEANCDD="$SCRIPT_BASE/clean.CDD.R";
CLEANAPPRIS="$SCRIPT_BASE/clean.appris.R";
PROCESS3P="$SCRIPT_BASE/process.3p.R";
CLEANDBSNP="$SCRIPT_BASE/clean.dbSNP.R";
PULLESE="$SCRIPT_BASE/pull.ESEs.R";
SCORENMD="$SCRIPT_BASE/score.NMD.R";
CLEANUNIPROT="$SCRIPT_BASE/clean.uniprot.R";

#1. Download hg38 and mm10 fasta, and build batmis index
for g in hg38 mm10; do
 mkdir -p $g && pushd "$_";
 wget --quiet -l 1 -e robots=off --wait 1 -r $UCSC_URL/goldenPath/$g/chromosomes #FTP is currently blocked
 pushd hgdownload.cse.ucsc.edu/goldenPath/$g/chromosomes/
 rm -f md5sum.txt README.txt *_alt.fa.gz index.html*
 for i in {1..22} X Y M; do if [ -e chr$i.fa.gz ]; then zcat chr$i.fa.gz; fi; done > ../../../../$g.fa
 for i in `ls | grep "random\|Un"`; do zcat $i; done >> ../../../../$g.fa
 popd
 rm -rf hgdownload.cse.ucsc.edu
 samtools dict $g.fa -o $g.samh
 samtools faidx $g.fa;
 
 mkdir batmis && pushd "$_";
 ln -s ../$g.fa $g.fa;
 $BATMIS/scripts/build_index $g.fa;
 samtools faidx $g.fa;
 popd
 
 popd
done

#2. Download and clean annotation
for g in hg38 mm10; do
 mkdir -p $g/anno && pushd "$_";
 case $g in
  "hg38") o="human"; s="Homo_sapiens"; s2="H_sapiens"; s3="homo_sapiens"; i="GRCh38"; p="7"; r="27";
          a="/ARCHIVE/ANNOTATION_RELEASE.108"; gcf="GCF_000001405.33"; ap="rs108v26";;
  "mm10") o="mouse"; s="Mus_musculus"; s2="M_musculus"; s3="mus_musculus"; i="GRCm38"; p="4"; r="M15";
          a=""; gcf="GCF_000001635.24"; ap="rs106v26"; esac;
 
 #GENCODE genes
 wget --quiet $GENCODE_URL/Gencode_$o/release_$r/gencode.v$r.primary_assembly.annotation.gff3.gz;
 wget --quiet $MAPPING_URL/$i\_gencode2UCSC.txt;
 zcat gencode.v$r.primary_assembly.annotation.gff3.gz | $GENCODE2UCSC $i\_gencode2UCSC.txt | gzip > gencode.gff3.gz;
 rm gencode.v$r.primary_assembly.annotation.gff3.gz $i\_gencode2UCSC.txt;
 $GFF2CDS gencode.gff3.gz gencode.CDS.Rdata gencode.CDS.bed;
 $CLEANGENCODE gencode.gff3.gz gencode.clean.Rdata;
 
 #Refseq genes
 #THE DEAL WITH REFSEQ:
 #1) CDD is not included in refseq releases. It is updated every week.
 #   It includes refseq proteins (aa fasta), which I have to use so the domains are correct
 #2) Refseq gene / mRNA / protein annotations (gtf), are released in releases.  Current : 109.
 #3) Appris is released in releases.  Current : 108.
 wget --quiet $NCBI_URL/genomes/$s$a/GFF/ref_$i.p$p\_top_level.gff3.gz; 
 wget --quiet $NCBI_URL/genomes/$s$a/Assembled_chromosomes/chr_accessions_$i.p$p;
 $REFSEQ2UCSC chr_accessions_$i.p$p ref_$i.p$p\_top_level.gff3.gz | gzip > refseq.gff3.gz
 $CLEANREFSEQ refseq.gff3.gz refseq.clean.gff3 refseq.clean.Rdata refseq.clean.uniq_cds.bed \
              refseq.isoform_collapse.bed refseq.isoform_collapse.Rdata
 grep -v ^# refseq.clean.gff3 | sort -t$'\t' -k1,1 -k4,4n -k5,5n | gzip > refseq.clean.gff3.gz
 rm chr_accessions_$i.p$p ref_$i.p$p\_top_level.gff3.gz refseq.clean.gff3
 
 #APPRIS principal isoforms
 #human
 wget $APPRIS_URL/$s3/$ap/appris_data.appris.txt
 $CLEANAPPRIS appris_data.appris.txt refseq.clean.Rdata $g appris.txt
 rm appris_data.appris.txt;
 
 #Pre-calculate 3' position in protein
 $PROCESS3P refseq.clean.Rdata appris.txt ThreePrime.Rdata
 
 #CDD - Updated weekly, and not archived
 #hg38: Run Feb 16, 2018
 #mm10: Run Feb 23, 2018
 mkdir -p cdd && pushd "$_";
 #THESE HAVE MITOCHONDRIAL GENES
 wget --quiet -l 1 -e robots=off --wait 1 -r $NCBI_URL/refseq/$s2/mRNA_Prot/;
 mv ftp.ncbi.nlm.nih.gov/refseq/$s2/mRNA_Prot/$o.* .;
 rm -r ftp.ncbi.nlm.nih.gov $o.files.installed;
 chunks=1;
 while [ -f $o.$[$chunks+1].protein.gpff.gz ]; do chunks=$[$chunks+1]; done;
 zcat `eval echo $o.{1..$chunks}.protein.gpff.gz` | $CLEANGENBANK | gzip > cdd.gpff.gz
 zcat `eval echo $o.{1..$chunks}.protein.faa.gz`  | gzip > ../refseq.aa.gz
 rm `eval echo $o.{1..$chunks}.{protein.gpff,rna.gbff,protein.faa,rna.fna}.gz`;
 java -jar $READSEQJAR -f=GFF -o=cdd.aa.gff2.o cdd.gpff.gz
 perl -ane 'print if($_ =~ m/^\#/ || $F[2] ne "variation");' cdd.aa.gff2.o | gzip > ../cdd.aa.gff2.gz
 rm cdd.gpff.gz cdd.aa.gff2.o
 popd
 rmdir cdd
 #THESE ARE MISSING MITOCHONDRIAL GENES
 wget --quiet $NCBI_URL/genomes/$s$a/Assembled_chromosomes/chr_accessions_$i.p$p
 wget --quiet $NCBI_URL/refseq/$s2/alignments/$gcf\_knownrefseq_alignments.gff3 -O map.gff3
 gzip map.gff3
 $REFSEQ2UCSC chr_accessions_$i.p$p map.gff3.gz | gzip > refseq.to.genome.gff3.gz
 rm chr_accessions_$i.p$p map.gff3.gz
 $CLEANCDD cdd.aa.gff2.gz refseq.to.genome.gff3.gz appris.txt \
  cdd.all.igv.gff3 cdd.isoform_collapse.igv.gff3 \
  cdd.principal.domains.igv.gff3 cdd.principal.sites.igv.gff3 \
  cdd.all.Rdata cdd.isoform_collapse.Rdata \
  cdd.principal.domains.Rdata cdd.principal.sites.Rdata;
 
 #DBSNP (hg38 only)
 if [ $g == hg38 ]; then
  wget --quiet $UCSC_URL/goldenPath/hg38/database/snp147Common.txt.gz;
  zcat snp147Common.txt.gz | cut -f 2-5,12,23,25 | gzip > dbSNP.common.147.txt.gz;
  $CLEANDBSNP dbSNP.common.147.txt.gz dbSNP.common.147.Rdata;
  rm dbSNP.common.147.txt.gz snp147Common.txt.gz;
 fi
 
 #Exonic splicing enhancers
 $PULLESE ../$g.fa refseq.clean.Rdata INT3.ESE.Rdata
 
 #Predict NMD
 $SCORENMD refseq.clean.Rdata appris.txt ../$g.fa NMD.bedgraph
 
 #Uniprot annotations (secondary structures)
 mkdir uniprot && pushd "$_";
 wget -r --no-parent -nd $UCSC_URL/gbdb/$g/uniprot/
 v=`cat version.txt | cut -d' ' -f 4`; #2017_09
 rm index.html* robots.txt version.txt
 for i in *.bb; do bigBedToBed $i /dev/stdout | \
  perl -ne 's/ /./g; s/["\047]//g; s/\t(?=[\t\n])/\tNA/g; print' > $i.bed; done
 rm *.bb
 cat *.bed | sort -k1,1 -k2,2n -k3,3n > ../uniprot.IGV.bed;
 $CLEANUNIPROT;
 mv uniprot.Rdata ..;
 rm *.bed;
 popd;
 rmdir uniprot;
 
 popd;
done

#missing - 
#PFAM is from the UCSC table browser
#Blast NR 2011 database
#aligner index builds

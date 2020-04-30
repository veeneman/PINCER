#!/usr/bin/env Rscript
#./clean.appris.R appris_data.appris.txt refseq.clean.Rdata output.txt
#This code uses Appris and Refseq to identify a single isoform, per gene, in this order
# APPRIS > Longest isoform > Lowest NP_#

library("stringr", "stringi", "rtracklayer", "GenomicAlignments", "GenomicFeatures");
argv = commandArgs(trailingOnly = TRUE);
#argv = c("appris_data.appris.txt", "refseq.clean.Rdata", "mm10", "appris.txt")
appris.fp = argv[1];
load(argv[2]);
species = argv[3];

#Load and set up appris
a = read.table(appris.fp, fill = TRUE);
a = a[a[,5] == "TRANSLATION",]; #drop XR, NR
a = a[,c(2:4,8,20)]; #drop columns
colnames(a) = c("gene","NM","NP","CCDS","APPRIS");
a = a[!(grepl("^XM",a$NM) | grepl("^XP",a$NP) | grepl("MINOR",a$APPRIS) | grepl("ALTERNATIVE",a$APPRIS)),];
a = unique(a);

#For each gene APPRIS didn't include, add all isoforms
g = unique(refseq$gene)
g = setdiff(g, a$gene)
nm = mcols(refseq[refseq$gene %in% g & refseq$type == "mRNA"])[,c("gene","Name","ID")];
np = setNames(refseq[refseq$type == "CDS"]$Name, refseq[refseq$type == "CDS"]$Parent)
tmp = setNames(data.frame(nm$gene, nm$Name, np[nm$ID], NA, NA), colnames(a))
a = rbind(a,tmp)

#Split APPRIS into unique and non-unique "principal" isoforms
a$NP.N = as.numeric(str_replace(a$NP, "^[^_]+_", ""));
a = a[order(a$gene, a$NP.N),];
dup = (a$gene %in% a$gene[duplicated(a$gene)]);
tmp = a[dup,];
a = a[!dup,];

#Drop non-unique "principal" isoforms which have exact same CDS as a
# lower-numbered principal isoform for the same gene
p = refseq[refseq$Name %in% tmp$NM];
p = p[!(p$Name %in% p$Name[duplicated(p$Name)]) | seqnames(p) == "chrX"]; #drop PAR on chrY
cds = refseq[refseq$Parent %in% p$ID & refseq$type == "CDS"];
cds = GRanges(paste(seqnames(cds), cds$gene), ranges(cds), strand(cds), mcols(cds));
cds$NP.N = as.numeric(str_replace(cds$Name, "^[^_]+_", ""));
cds.split = split(cds[,NULL], cds$NP.N);
o = findOverlaps(cds.split, type = "within", maxgap = 0, minoverlap = 1);
drop = as.numeric(names(cds.split)[unique(from(o[isRedundantHit(o)]))]);
tmp = tmp[!(tmp$NP.N %in% drop),]
dup = (tmp$gene %in% tmp$gene[duplicated(tmp$gene)]);
a = rbind(a, tmp[!dup,]);
tmp = tmp[dup,]

#Finally, for genes with multiple "principal" isoforms,
# which don't have the exact same CDS,use 
# either a) the longest, or b) the lowest numbered NP_#, in that order
tmp$length = sum(width(cds.split))[as.character(tmp$NP.N)]
tmp = tmp[order(tmp$gene, -tmp$length, tmp$NP.N),]
a = rbind(a, tmp[!duplicated(tmp$gene),][,colnames(a)])
rm(tmp)

#Remove readthrough transcripts, which are Not genes
if(species == "hg38") {
  rt = grepl("-", a$gene) & #find all hyphenated gene names
       !grepl("-\\d+$|^HLA|-AS\\d+", a$gene) & #exclude numbered genes, hla genes, antisense genes
       !a$gene %in% c("WI2-2373I1.2","GS1-259H13.2","CH17-360D5.1","MIR1-1HG") #specific non-readthroughs
  a = a[!rt,];
} else if(species == "mm10") {
  rt = unique(refseq[grepl("readthrough",refseq$description)]$gene);
  a = a[!a$gene %in% rt,];
}

#Write out output
a = a[order(a$gene, a$NP),c("gene","NM","NP","APPRIS")];
write.table(a, file = argv[4], row.names = F, quote = F, sep = "\t");

#CHECKING WHETHER APPRIS NMS AND NPS MATCH:
##1: Trying with NMs
#p = refseq[refseq$Name %in% b$NM]
#p = p[!(p$Name %in% p$Name[duplicated(p$Name)]) | seqnames(p) == "chrX"] #drop PAR on chrY
#v1 = refseq[refseq$Parent %in% p$ID & refseq$type == "CDS"]
##2: Trying with NPs
#v2 = refseq[refseq$Name %in% b$NP]
#i = do.call(rbind,str_split(unique(paste(seqnames(v2), v2$Name))," "))
#d = apply(i[i[,2] %in% i[,2][duplicated(i[,2])] & i[,1] == "chrY",], 1, paste, collapse = " ")
#v2 = v2[!(paste(seqnames(v2), v2$Name) %in% d)]
#identical(v1, v2)
##[1] TRUE #YES

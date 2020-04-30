#!/usr/bin/env Rscript
library("Biostrings","GenomicRanges", "rtracklayer");

#Usage: ./pull.ESEs.R genome.fa refseq.clean.Rdata output.rdata
argv = commandArgs(trailingOnly = TRUE);
ese.fp = "~/afs/software/lib/ESE.INT3.txt";
genome.fp = argv[1];
refseq.fp = argv[2];
output.fp = argv[3];

g = readDNAStringSet(genome.fp);
g = g[names(g) %in% paste0("chr",c(1:22,"X","Y","M"))];

x = PDict(read.table(ese.fp)[,1]);
y = lapply(g, function(i) { unlist(IRangesList(matchPDict(x, i))); });
z = unlist(IRangesList(y));

result = GRanges(names(z), z, "*");
result = sortSeqlevels(result);
result = result[order(as.character(seqnames(result)), start(result), end(result))];

load(refseq.fp);
refseq = refseq[refseq$type == "exon"];

o = findOverlaps(result, refseq, type = "within", maxgap = 0, minoverlap = 1);
result = result[unique(from(o))];
result = reduce(result);
ese = result;

save(ese, file = output.fp);

#NB - this is file specific
library("stringr", "stringi", "rtracklayer", "GenomicAlignments", "GenomicFeatures");
load("cdd.principal.domains.Rdata")

#try2
first    = !duplicated(cdd.p$domain.id);
last     = !duplicated(cdd.p$domain.id, fromLast = T);
starts   = setNames(start(cdd.p)[first] -1, cdd.p$domain.id[first])
ends     = setNames(end(cdd.p)[last], cdd.p$domain.id[last])
seqnames = setNames(as.character(seqnames(cdd.p)[first]), cdd.p$domain.id[first])
strands  = setNames(as.character(strand(cdd.p)[first]), cdd.p$domain.id[first])
names    = setNames(cdd.p$domain.name[first], cdd.p$domain.id[first])
cdd.p$blockSize = width(cdd.p)
cdd.p$blockStart = (start(cdd.p) -1) - starts[as.character(cdd.p$domain.id)]
blockCounts = sapply(split(cdd.p$blockSize, cdd.p$domain.id), length)
blockSizes  = sapply(split(cdd.p$blockSize, cdd.p$domain.id), paste, collapse = ",")
blockStarts = sapply(split(cdd.p$blockStart, cdd.p$domain.id), paste, collapse = ",")
u = unique(as.character(cdd.p$domain.id))

x = data.frame("chr" = seqnames[u], "start" = starts[u], "end" = ends[u],
    "domain.name" = names[u], "score" = rep(1000, length(u)),
    "strand" = strands[u], "thickstart" = starts[u], "thickend" = ends[u],
    "color" = rep("0,0,0", length(u)),
    "blockCounts" = blockCounts[u], "blockSizes" = blockSizes[u], "blockStarts" = blockStarts[u])
write.table(x, "cdd.principal.domains.igv.bed", row.names = F, col.names = F, sep = "\t", quote = F)
q()

library("stringr", "stringi", "rtracklayer", "GenomicAlignments", "GenomicFeatures");
load("refseq.clean.Rdata")

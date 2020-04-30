#!/usr/bin/env Rscript
#./gff2cds.R input.gff3.gz output.cds.Rdata
suppressMessages(library(rtracklayer));
argv = commandArgs(trailingOnly = TRUE);
gff = as(import.gff3(argv[1]), "GRanges");
cds = gff[gff$type == "CDS"];
cds = reduce(cds); #collapse all CDSs together
cds = cds +2; #splice sites, known bug - start/end codon also
save(cds, file = argv[2]);
df = data.frame(seqnames(cds), start(cds)-1, end(cds), ".", ".", strand(cds));
write.table(df, file = argv[3], quote = F, sep = "\t", row.names = F, col.names = F);
q();

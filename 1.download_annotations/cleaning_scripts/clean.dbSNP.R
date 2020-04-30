#!/usr/bin/env Rscript
library("stringr", "stringi", "rtracklayer", "GenomicAlignments", "GenomicFeatures");
argv = commandArgs(trailingOnly = TRUE);

x = read.table(argv[1], sep = "\t", header = F)
dbSNP = GRanges(x[,1], IRanges(x[,2] +1, x[,3]), "*",
                DataFrame("id" = x[,4], "type" = x[,5], "alleles" = x[,6], "freq" = x[,7]));
dbSNP$freq = NumericList(strsplit(dbSNP$freq, ","));
dbSNP$alleles = CharacterList(strsplit(dbSNP$alleles, ","));

save(dbSNP, file = argv[2])
q()

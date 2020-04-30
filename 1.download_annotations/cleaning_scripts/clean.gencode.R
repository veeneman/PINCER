#!/usr/bin/env Rscript
#./clean.gencode.R input.gff3.gz output.Rdata
calc.stats = function(x) {
  genes = length(x[x$type == "gene",]);
  transcripts = length(x[x$type == "transcript",]);
  nt.u = sum(width(reduce(x[x$type == "exon", ])));
  return(c(genes, transcripts, nt.u))
}

suppressMessages(library(rtracklayer))
argv = commandArgs(trailingOnly = TRUE);

#Load Gencode genes
gff = import.gff3(argv[1]); 
result = data.frame(matrix(nrow=0,ncol=3,dimnames=list(NULL,c("genes","transcripts","nt.u"))))
result["full",] = calc.stats(gff)

#filter for chr[1-22,X,Y,M]
chrs = paste0("chr",c(1:22,"X","Y","M"))
gff = gff[seqnames(gff) %in% chrs,]
result["primary.chrs",] = calc.stats(gff)

#filter for gencode-basic transcripts, and their parents (genes)
basic.transcripts = gff[gff$type == "transcript" & any(gff$tag %in% "basic"),]
g = unique(basic.transcripts$gene_id)
t = unique(basic.transcripts$transcript_id)
gff = gff[(gff$type == "gene" & gff$gene_id %in% g) | gff$transcript_id %in% t,]
result["basic",] = calc.stats(gff)
rm(g, t)

#filter for protein-coding transcripts, and their parents (genes)
pc.transcripts = gff[gff$type == "transcript" & gff$transcript_type %in% "protein_coding",]
g = unique(pc.transcripts$gene_id)
t = unique(pc.transcripts$transcript_id)
gff = gff[(gff$type == "gene" & gff$gene_id %in% g) | gff$transcript_id %in% t,]
result["protein.coding",] = calc.stats(gff)
rm(g, t)

#filter for TSL1 transcripts, and their parents (genes)
tsl1.transcripts = gff[gff$type == "transcript" & !(gff$transcript_support_level %in% c(2:5)),]
g = unique(tsl1.transcripts$gene_id)
t = unique(tsl1.transcripts$transcript_id)
gff = gff[(gff$type == "gene" & gff$gene_id %in% g) | gff$transcript_id %in% t,]
result["tsl1",] = calc.stats(gff)
rm(g, t)

#calculate unique CDS nucleotides for summary, and print it to stdout
result["cds",] = result[nrow(result),]
result["cds","nt.u"] = sum(width(reduce(gff[gff$type == "CDS", ])));
#print(result);

gencode = gff
save(gencode, file = argv[2]);

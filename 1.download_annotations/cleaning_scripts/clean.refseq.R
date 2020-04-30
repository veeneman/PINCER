#!/usr/bin/env Rscript
#./clean.refseq.R in.gff3.gz out.clean.gff3 out.clean.Rdata out.clean.uniq_cds.bed \
#                 out.isoform_collapse.bed out.isoform_collapse.Rdata
#in order to work correctly, this
# relies on NM_* mapping correctly to the mRNA's name, and
# relies on parent<>child ID relationships being correct
suppressMessages(library(rtracklayer));
argv = commandArgs(trailingOnly = TRUE);
argv = c("refseq.gff3.gz",
         "refseq.clean.gff3", "refseq.clean.Rdata", "refseq.clean.uniq_cds.bed",
         "refseq.isoform_collapse.bed", "refseq.isoform_collapse.Rdata")
gff = import.gff3(argv[1]);

#filter for chr[1-22,X,Y,M]
chrs = paste0("chr",c(1:22,"X","Y","M"))
gff = gff[seqnames(gff) %in% chrs,]
seqlevels(gff) = levels(factor(seqnames(gff)))

#filter for NM_ and YP_ transcripts, which means both curated and protein-coding
g = unique(gff[gff$type == "gene",]$ID)
nm.t = gff[any(gff$Parent %in% g) & grepl("^NM_|^YP_",gff$Name),]$ID
nm.g = unique(unlist(gff[any(gff$Parent %in% g) & grepl("^NM_|^YP_",gff$Name),]$Parent))
gff = gff[gff$ID %in% c(nm.g, nm.t) | any(gff$Parent %in% nm.t),]

#drop useless metadata
columns.keep = c("source", "type", "phase", "ID",
                 "Dbxref", "Name", "gbkey", "description",
                 "gene", "gene_biotype", "Parent", "product", 
                 "transcript_id", "gene_synonym", "protein_id", "Note")
mcols(gff) = mcols(gff)[,columns.keep]

#fix Parent field
if(any(lengths(gff$Parent) > 1)) { write("Error: Parent field", stderr()); q(); }
gff$Parent[lengths(gff$Parent) == 0] = "NA";
gff$Parent = unlist(gff$Parent);

#write PRI/PC/Curated Refseq (aka Clean)
refseq = gff;
export(refseq, argv[2], format = "GFF3");
save(refseq, file = argv[3]);
rm(refseq);

#write full CDS file, usable in off-target calculation
# (though I used gencode before)
cds = gff[gff$type == "CDS"];
cds = reduce(cds); #collapse all CDSs together
cds = cds +2; #splice sites, known bug - start/end codon also
df = data.frame(seqnames(cds), start(cds)-1, end(cds), ".", ".", strand(cds));
write.table(df, file = argv[4], quote = F, sep = "\t", row.names = F, col.names = F);

#collapse genes by isoform-constitutive CDS,
# and save bed and Rdata
# nb - we are Not targeting splice sites
gene = unlist(split(gff[gff$type == "gene"]$Name, gff[gff$type == "gene"]$ID))
mrna = split(gff[gff$type == "mRNA"]$ID, gff[gff$type == "mRNA"]$Parent)
cds = split(gff[gff$type == "CDS",NULL], gff[gff$type == "CDS"]$Parent)
refseq.shared = unlist(GRangesList(sapply(mrna, function(m) { Reduce(intersect, cds[m]); })));
mcols(refseq.shared)$ID = names(refseq.shared)
mcols(refseq.shared)$gene = gene[names(refseq.shared)]
names(refseq.shared) = NULL

#rescue 13 mitochondrial genes (for which CDS's parent is the gene directly):
chrM = gff[seqnames(gff) %in% "chrM" & gff$type == "CDS"][,c(11,9)]
colnames(mcols(chrM)) = c("ID","gene")
refseq.shared = c(refseq.shared,chrM)

df = as.data.frame(refseq.shared)[,c(1:3,5,7)];
write.table(df, file = argv[5], quote = F, sep = "\t", row.names = F, col.names = F)
save(refseq.shared,file = argv[6])
q();

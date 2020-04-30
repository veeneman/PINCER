#!/usr/bin/env Rscript
#./score.NMD.R refseq.clean.Rdata appris.txt genome.fa
library("stringr", "stringi", "rtracklayer", "GenomicAlignments", "GenomicFeatures", "Rsamtools");

warn = function(msg) {
  write(paste("Warning:",msg), stderr());
}
drop.par = function(g) {
  gene.lookup = unique(as.data.frame(g)[,c("seqnames","gene.id","gene.name")])
  gene.ct = table(gene.lookup$gene.name)
  par = names(gene.ct[gene.ct > 1])
  return(g[!(seqnames(g) == "chrX" & g$gene.name %in% par)]);
}
load.refseq.p.cds = function(REFSEQ.ALL.FP, APPRIS.FP) {
  load(REFSEQ.ALL.FP);
  appris = read.table(APPRIS.FP, header = T);
  appris$gene.id = as.character(1:nrow(appris));
  
  refseq = refseq[refseq$type == "CDS" & refseq$Name %in% appris$NP]
  id.lookup = setNames(appris$gene.id, appris$NP)
  refseq$gene.name = refseq$gene;
  refseq$gene.id = id.lookup[refseq$protein_id]
  refseq = refseq[,c("gene.id","gene.name","phase")]
  refseq = drop.par(refseq);
  refseq$gene.id = as.numeric(refseq$gene.id);
  return(refseq);
}
load.refseq.p.exons = function(REFSEQ.ALL.FP, APPRIS.FP) {
  load(REFSEQ.ALL.FP);
  appris = read.table(APPRIS.FP, header = T);
  appris$gene.id = as.character(1:nrow(appris));
  
  refseq = refseq[refseq$type == "exon" & refseq$transcript_id %in% appris$NM]
  id.lookup = setNames(appris$gene.id, appris$NM)
  refseq$gene.name = refseq$gene;
  refseq$gene.id = id.lookup[refseq$transcript_id]
  refseq = refseq[,c("gene.id","gene.name")]
  refseq = drop.par(refseq);
  refseq$gene.id = as.numeric(refseq$gene.id);
  return(refseq);
}
shift.seqnames = function(g, id) {
  g$chr.o = as.character(seqnames(g));
  g = GRanges(paste(mcols(g)[,id], g$chr.o), ranges(g), strand(g), mcols(g));
  return(g);
}
unshift.seqnames = function(g, id) {
  g = GRanges(mcols(g)[,id], ranges(g), strand(g), mcols(g)[,!colnames(mcols(g)) == "chr.o"]);
  return(g);
}
cDNA.last.sj = function(refseq.exons) { #1-indexed
  #exon widths
  widths = split(width(refseq.exons), refseq.exons$gene.id);
  #sum of exon widths, ignoring the last exon
  last.sjs = sapply(widths, function(i) { sum(i[-length(i)]); });
  last.sjs = last.sjs +1;
  return(last.sjs);
}
cDNA.offsets = function(refseq.exons) { #0-indexed
  #set cDNA offsets
  refseq.exons$id.o = 1:length(refseq.exons);
  refseq.exons = refseq.exons[order(refseq.exons$gene.id, refseq.exons$id.o)];
  tmp = split(width(refseq.exons), refseq.exons$gene.id);
  tmp2 = lapply(tmp, function(i) { c(0, cumsum(i[-length(i)])); });
  refseq.exons$cDNA.offset = unname(do.call(c, tmp2));
  mcols(refseq.exons) = mcols(refseq.exons)[,!colnames(mcols(refseq.exons)) == "id.o"]
  return(refseq.exons);
}
cDNA.cds.s = function(refseq.exons, refseq.cds) { #1-indexed
  #compute cDNA coordinates of CDS starts and ends
  
  #find genomic coordinates of CDS starts
  cds.sg = refseq.cds[!duplicated(refseq.cds$gene.id)];
  end(  cds.sg[strand(cds.sg) == "+"]) = start(cds.sg[strand(cds.sg) == "+"])
  start(cds.sg[strand(cds.sg) == "-"]) =   end(cds.sg[strand(cds.sg) == "-"])
  
  #shift CDS starts from genome to cDNA
  o = findOverlaps(cds.sg, refseq.exons, maxgap = 0, minoverlap = 1)
  cds.s = (start(cds.sg[from(o)]) - start(refseq.exons[to(o)])) + refseq.exons[to(o)]$cDNA.offset; #+
  tmp =   (end(refseq.exons[to(o)]) - end(cds.sg[from(o)]))     + refseq.exons[to(o)]$cDNA.offset; #-
  i = as.logical(strand(cds.sg[from(o)]) == "-");
  cds.s[i] = tmp[i];
  names(cds.s) = cds.sg[from(o)]$gene.id;
  cds.s = cds.s +1; #0-index to 1-index
  return(cds.s);
}
cDNA.cds.e = function(refseq.exons, refseq.cds) { #1-indexed
  #find genomic coordinates of CDS ends
  cds.eg = refseq.cds[!duplicated(refseq.cds$gene.id, fromLast = TRUE)];
  start(cds.eg[strand(cds.eg) == "+"]) =   end(cds.eg[strand(cds.eg) == "+"]);
  end(  cds.eg[strand(cds.eg) == "-"]) = start(cds.eg[strand(cds.eg) == "-"]);
  
  #shift CDS ends from genome to cDNA
  o = findOverlaps(cds.eg, refseq.exons, maxgap = 0, minoverlap = 1);
  cds.e = (start(cds.eg[from(o)]) - start(refseq.exons[to(o)])) + refseq.exons[to(o)]$cDNA.offset; #+
  tmp =   (end(refseq.exons[to(o)]) - end(cds.eg[from(o)]))     + refseq.exons[to(o)]$cDNA.offset; #-
  i = as.logical(strand(cds.eg[from(o)]) == "-");
  cds.e[i] = tmp[i];
  names(cds.e) = cds.eg[from(o)]$gene.id;
  cds.e = cds.e -3; #Drop the stop codon
  cds.e = cds.e +1; #0-index to 1-index
  return(cds.e);
}

argv = commandArgs(trailingOnly = TRUE);
REFSEQ.ALL.FP = argv[1];
APPRIS.FP     = argv[2];
GENOME.FP     = argv[3];
OUTPUT.FP     = argv[4];

#load refseq cds, and mRNA
refseq.cds   = load.refseq.p.cds(REFSEQ.ALL.FP, APPRIS.FP);
refseq.exons = load.refseq.p.exons(REFSEQ.ALL.FP, APPRIS.FP);
refseq.exons = c(refseq.exons, refseq.cds[seqnames(refseq.cds) == "chrM",1:2]);
refseq.cds   = shift.seqnames(refseq.cds, "gene.id");
refseq.exons = shift.seqnames(refseq.exons, "gene.id");
refseq.exons = cDNA.offsets(refseq.exons);

#pull cDNA coordinates
last.sjs = cDNA.last.sj(refseq.exons);
cds.s = cDNA.cds.s(refseq.exons, refseq.cds);
cds.e = cDNA.cds.e(refseq.exons, refseq.cds);
refseq.exons = unshift.seqnames(refseq.exons, "chr.o");

#create reverse cDNA map
#per-range offset, adjusted for stranded cDNA positions
offset = start(refseq.exons) - refseq.exons$cDNA.offset;
i = as.character(strand(refseq.exons)) == "-";
offset[i] = end(refseq.exons)[i] + refseq.exons$cDNA.offset[i];
map = GRanges(refseq.exons$gene.id,
              IRanges(refseq.exons$cDNA.offset +1,
                  refseq.exons$cDNA.offset + width(refseq.exons)),
              "+",
              DataFrame(chr = seqnames(refseq.exons),
                        offset = offset,
                        str = strand(refseq.exons)));

#pull cDNA sequences, minus last exons, which are large and unnecessary
genome = FaFile(GENOME.FP);
i = duplicated(refseq.exons$gene.id, fromLast = TRUE);
cDNA.split = GRangesList(split(refseq.exons[i], refseq.exons[i]$gene.id))
cDNA.all = extractTranscriptSeqs(genome, cDNA.split); #gets correct strand

#handle single-exon genes, which we'll concatenate on later
single.exon = setdiff(unique(refseq.exons$gene.id), names(cDNA.all)); #dropped when we drop last exon
single.exon = refseq.cds[refseq.cds$gene.id %in% single.exon];
single.exon = unshift.seqnames(single.exon, "chr.o");
single.exon$score = 0;
single.exon = single.exon[,"score"];
end(  single.exon[strand(single.exon) == "+"]) = end(  single.exon[strand(single.exon) == "+"]) -3;
start(single.exon[strand(single.exon) == "-"]) = start(single.exon[strand(single.exon) == "-"]) +3;
rm(refseq.cds);

#run actual NMD prediction
ids = names(cDNA.all);
hits = gregexpr(cDNA.all, pattern = "TAA|TGA|TAG");
names(hits) = names(cDNA.all);

x = DataFrame("hits" = IntegerList(hits[ids]),
              "last.sjs" = last.sjs[ids],
              "cds.s" = cds.s[ids],
              "cds.e" = cds.e[ids],
              row.names = ids);
x$pos = x$hits[x$hits < (x$last.sjs - 50) & x$hits > x$cds.s & x$hits < x$cds.e]
x$phase = (x$pos - x$cds.s) %% 3

#warn about genes with in-frame PTCs
ptc = row.names(x[any(x$phase == 0),]);
ptc = unique(refseq.exons[refseq.exons$gene.id %in% ptc]$gene.name);
warn(paste(c("These genes have in-frame PTCs:",ptc), collapse = " "));

#find last frame1 and frame2 PTCs, and sort them by coord
x$ptc1 = max(x$pos[x$phase == 1]);
x$ptc2 = max(x$pos[x$phase == 2]);
x[x$ptc1 < 0,"ptc1"] = x[x$ptc1 < 0,"cds.s"];
x[x$ptc2 < 0,"ptc2"] = x[x$ptc2 < 0,"cds.s"];
x[x$ptc1 > x$ptc2, c("ptc1","ptc2")] = x[x$ptc1 > x$ptc2, c("ptc2","ptc1")]; #reorder

#debug: tag non-cds
#cov = GRanges(rep(row.names(x), 4),
#              IRanges(c(rep(1, nrow(x)), x$cds.s, x$ptc1, x$ptc2),
#                      c(x$cds.s -1, x$ptc1 -1, x$ptc2 -1, x$cds.e)),
#              "+",
#              DataFrame("score" = c(rep(4,nrow(x)), rep(3,nrow(x)), rep(2,nrow(x)), rep(1,nrow(x)))));
cov = GRanges(rep(row.names(x), 3),
              IRanges(c(x$cds.s, x$ptc1, x$ptc2),
                      c(x$ptc1 -1, x$ptc2 -1, x$cds.e)),
              "+",
              DataFrame("score" = c(rep(1,nrow(x)), rep(0.5,nrow(x)), rep(0,nrow(x)))));
o = findOverlaps(cov, map, maxgap = 0, minoverlap = 1)
val = pintersect(cov[from(o)], map[to(o)])
mcols(val) = cbind(mcols(cov[from(o)]), mcols(map[to(o)]))
result = GRanges(val$chr,
                 IRanges(start(val) + val$offset -1, end(val) + val$offset -1),
                 val$str,
                 DataFrame("score" = val$score));
neg.g = val$str == "-";
ranges(result[neg.g]) = IRanges(val[neg.g]$offset - end(val[neg.g]) +1,
                                val[neg.g]$offset - start(val[neg.g]) +1)
result = c(result, single.exon);
result = sortSeqlevels(result)
result = result[order(as.character(seqnames(result)), start(result), end(result))]
export(result, OUTPUT.FP, format = "bedgraph")

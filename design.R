#!/usr/bin/env Rscript
#TODO - switch to primer3 tm
#TODO - gene overlap handling for the three protein data tracks: nmd, 3', and conservation
date();
library("stringr", "stringi", "grid", "gridBase", "gridExtra",
        "rtracklayer", "GenomicAlignments", "GenomicFeatures", "HELP");

#Parameters & Filepaths
argv = commandArgs(trailingOnly = TRUE);
SPECIES            = argv[1]; #"mm10"; #"hg38"; #	
GUIDES.PER.GENE    = 6;
MAX.GUIDE.PERMUT   = 15; #heuristic, doesn't affect result unless it's too low
MIN.GUIDE.DIST     = 3;
POLYN              = setNames(c(5, 5, 5, 3), c("A", "C", "G", "T")); #minimum dis-allowed homopolymer length
RESTRICTION.MOTIF  = "CGTCTC";       #restriction endonuclease cut site, for filtering out guides
GUIDE.PAD.LEFT     = "CGTCTCTCACCG"; #sequence immediately left  of the guide in the synthesized oligo
GUIDE.PAD.RIGHT    = "GTTTGGAGACG";  #sequence immediately right of the guide in the synthesized oligo
GUIDES.FP          = paste0("~/hpc/guideDB/",SPECIES,".CDS.bam");
REFSEQ.ALL.FP      = paste0("~/hpc/ref/",SPECIES,"/anno/refseq.clean.Rdata");
REFSEQ.ISO.FP      = paste0("~/hpc/ref/",SPECIES,"/anno/refseq.isoform_collapse.Rdata");
APPRIS.FP          = paste0("~/hpc/ref/",SPECIES,"/anno/appris.txt");
THREEPRIME.FP      = paste0("~/hpc/ref/",SPECIES,"/anno/ThreePrime.Rdata");
CDD.C.FP           = paste0("~/hpc/ref/",SPECIES,"/anno/cdd.isoform_collapse.Rdata");
CDD.P.FP           = paste0("~/hpc/ref/",SPECIES,"/anno/cdd.principal.domains.Rdata");
CDD.S.FP           = paste0("~/hpc/ref/",SPECIES,"/anno/cdd.principal.sites.Rdata");
PFAM.FP            = paste0("~/hpc/ref/",SPECIES,"/anno/pfam.ucsc.bed.gz");
UNIPROT.FP         = paste0("~/hpc/ref/",SPECIES,"/anno/uniprot.Rdata");
PROVEAN.FP         = paste0("~/hpc/ref/",SPECIES,"/anno/provean.bedgraph");
INT3.ESE.FP        = paste0("~/hpc/ref/",SPECIES,"/anno/INT3.ESE.Rdata");
NMD.FP             = paste0("~/hpc/ref/",SPECIES,"/anno/NMD.bedgraph");
SNP.FP             = paste0("~/hpc/ref/",SPECIES,"/anno/dbSNP.common.147.Rdata");
GENOME.FP          = paste0("~/hpc/ref/",SPECIES,"/",SPECIES,".fa");
VERSION            = "1.3c";

load.guides = function(fp) {
  tags = c("CS","AM","AZ","OL","HG","HP", paste0("G",0:4), paste0("P",0:4));
  guides = GRanges(import(fp, param = ScanBamParam(what = "seq", tag = tags)))
  neg = strand(guides) == "-";
  mcols(guides[neg])$seq = reverseComplement(mcols(guides[neg])$seq)
  seqlevels(guides) = levels(factor(seqnames(guides)))
  return(guides);
}
drop.par = function(g) {
  gene.lookup = unique(as.data.frame(g)[,c("seqnames","gene.id","gene.name")])
  gene.ct = table(gene.lookup$gene.name)
  par = names(gene.ct[gene.ct > 1])
  return(g[!(seqnames(g) == "chrX" & g$gene.name %in% par)]);
}
load.refseq.p = function(REFSEQ.ALL.FP, APPRIS.FP) {
  load(REFSEQ.ALL.FP);
  appris = read.table(APPRIS.FP, header = T);
  appris$gene.id = as.character(1:nrow(appris));
  refseq = refseq[refseq$type == "CDS" & refseq$Name %in% appris$NP]
  id.lookup = setNames(appris$gene.id, appris$NP)
  refseq$gene.name = refseq$gene;
  refseq$gene.id = id.lookup[refseq$protein_id];
  refseq$entrez.id = sapply(refseq$Dbxref, function(v) { v[grepl("^GeneID:", v)]; });
  refseq$entrez.id = str_replace(refseq$entrez.id, "^GeneID:", "");
  refseq = refseq[,c("gene.id","gene.name","entrez.id")]
  refseq = drop.par(refseq);
  refseq$gene.id = as.numeric(refseq$gene.id);
  refseq$entrez.id = as.numeric(refseq$entrez.id);
  colnames(mcols(refseq)) = paste0("principal.", colnames(mcols(refseq)));
  refseq$first = !duplicated(refseq$principal.gene.id);
  refseq$last  = !duplicated(refseq$principal.gene.id, fromLast = TRUE);
  return(refseq);
}
load.refseq.c = function(REFSEQ.ISO.FP, APPRIS.FP) {
  load(REFSEQ.ISO.FP);
  names(mcols(refseq.shared)) = c("gene.id", "gene.name");
  appris = read.table(APPRIS.FP, header = T); #redo IDs to match principal isoforms
  appris$gene.id = as.character(1:nrow(appris));
  id.lookup = setNames(appris$gene.id, appris$gene)
  refseq.shared$gene.id = id.lookup[refseq.shared$gene.name];
  refseq.shared = drop.par(refseq.shared);
  refseq.shared$gene.id = as.numeric(refseq.shared$gene.id);
  colnames(mcols(refseq.shared)) = paste0("isoform_shared.", colnames(mcols(refseq.shared)));
  return(refseq.shared);
}
load.threeprime = function(fp) {
  load(fp);
  return(ThreePrime[,-1]);
}
load.cdd.p = function(fp) {
  load(fp);
  mcols(cdd.p) = mcols(cdd.p)[,c("domain.name", "note", "domain.id")]
  colnames(mcols(cdd.p)) = paste0("cdd.p.", c("name", "description", "id"));
  return(cdd.p);
}
load.cdd.c = function(fp) {
  load(fp);
  mcols(cdd.c) = mcols(cdd.c)[,c("domain.name", "note", "domain.id")]
  colnames(mcols(cdd.c)) = paste0("cdd.c.", c("name", "description", "id"));
  return(cdd.c);
}
load.cdd.s = function(fp) {
  load(fp);
  mcols(sites.p) = DataFrame("cdd.p.site" = paste(sites.p$name, sites.p$note, sep = "; "))
  return(sites.p);
}
load.pfam = function(fp) {
  pfam.o = import.bed(fp);
  pfam.o = pfam.o[seqnames(pfam.o) %in% paste0("chr",c(1:22,"X","Y","M")),];
  seqlevels(pfam.o) = seqlevels(pfam.o)[1:25];
  b = blocks(pfam.o);
  names(b) = 1:length(b);
  pfam = unlist(b);
  pfam$domain.id = as.numeric(names(pfam));
  pfam$domain.name = pfam.o$name[pfam$domain.id]
  names(pfam) = NULL;
  colnames(mcols(pfam)) = c("pfam.id", "pfam.name");
  mcols(pfam) = mcols(pfam)[c("pfam.name", "pfam.id")];
  return(pfam);
}
load.uniprot.2nd = function(fp) {
  load(fp);
  uniprot = uniprot[grepl("Secondary.structure", uniprot$name)];
  uniprot$beta  = grepl("Secondary.structure:beta.strand", uniprot$name);
  uniprot$ccr   = grepl("Secondary.structure:coiled-coil.region", uniprot$name);
  uniprot$helix = grepl("Secondary.structure:helix", uniprot$name);
  uniprot$turn  = grepl("Secondary.structure:turn", uniprot$name);
  uniprot = uniprot[,c("beta","ccr","helix","turn")]
  colnames(mcols(uniprot)) = paste0("uniprot.",colnames(mcols(uniprot)));
  return(uniprot);
}
load.provean = function(fp) {
  provean = import(fp);
  colnames(mcols(provean)) = "provean"
  provean$provean = provean$provean * -1
  return(provean);
}
load.ese = function(fp) {
  load(INT3.ESE.FP);
  ese$ese = TRUE
  return(ese);
}
load.nmd = function(fp) {
  nmd = import(fp);
  colnames(mcols(nmd)) = paste0("nmd.",colnames(mcols(nmd)));
  return(nmd);
}
tag.restriction.sites = function(guides) {
  #Somehow this was faster than a single regex search
  # It shouldn't be - R's regexes are probably flawed somehow
  # I also couldn't get it working natively with DNAStrings (i.e., in 2bit)
  pad.L      = DNAString(GUIDE.PAD.LEFT);
  pad.R      = DNAString(GUIDE.PAD.RIGHT);
  restrict.F = RESTRICTION.MOTIF;
  restrict.R = as.character(reverseComplement(DNAString(RESTRICTION.MOTIF)));
  seq = as.character(xscat(pad.L, guides$seq, pad.R));
  ct.F = stri_count_fixed(seq, restrict.F, overlap = TRUE);
  ct.R = stri_count_fixed(seq, restrict.R, overlap = TRUE);
  guides$RS = (ct.F > 1 | ct.R > 1);
  return(guides);
}
tag.homopolymers = function(guides) {
  pad.L = DNAString(str_match(GUIDE.PAD.LEFT,  "(.)\\1*$")[,1]);
  pad.R = DNAString(str_match(GUIDE.PAD.RIGHT, "^(.)\\1*")[,1]);
  sequence = as.character(xscat(pad.L, guides$seq, pad.R));
  guides$PA = longestConsecutive(sequence, "A");
  guides$PC = longestConsecutive(sequence, "C");
  guides$PG = longestConsecutive(sequence, "G");
  guides$PT = longestConsecutive(sequence, "T");
  return(guides);
}
oligoTm.primer3 = function(guides) {
  x = tempfile();
  y = tempfile();
  write(as.character(guides$seq), x);
  system(paste0("for i in `cat ", x, "`; do oligotm -tp 1 -sc 1 $i; done >", y));
  Tm = read.table(y, header = F)[,1];
  if(length(Tm) != length(guides)) { write("ERROR in OligoTm", stderr()); q(); }
  unlink(c(x, y));
  guides$Tm.p3 = Tm;
  return(guides);
}
tag.snps = function(guides, refseq.p, GENOME.FP, SNP.FP) {
  guides$dbSNP.rs  = CharacterList(NULL);
  guides$dbSNP.vaf = NumericList(NULL);
  if(SPECIES != "hg38") { return(guides); }
  load(SNP.FP);
  
  #filter for SNPs within 30nt of CDS, for speed
  o = findOverlaps(dbSNP, refseq.p +30, maxgap = 0, minoverlap = 0, ignore.strand = TRUE);
  dbSNP = dbSNP[unique(from(o))];
  
  #figure out minor allele frequencies
  genome = FaFile(GENOME.FP);
  dbSNP$hg38 = as.character(getSeq(genome, dbSNP))
  dbSNP[width(dbSNP) == 0]$hg38 = "-";
  #this is complicated:
  # for SNPs, it sums the bases dbSNP says are at the site, that don't match hg38
  # for Deletions (in the individual relative to the reference), it sums EITHER
  #  alleles that are deletions, OR alleles that don't match hg38
  # for Insertions (in the individual relative to the reference), it sums
  #  any allele that isn't "-".  It's impossible to know if the
  #  reference instead matches the insertion
  dbSNP$vaf = sum(dbSNP$freq[dbSNP$alleles != dbSNP$hg38]);
  
  #PAMs also can't overlap guides
  pams = flank(flank(guides[,NULL], 3, start=FALSE), -2, start=F)
  
  #check whether 20nt guides, or their 20nt PAMs overlap common SNPs,
  # and merge the hits.  PAMs and guides have the same sort order so it works.
  # the degenerate N nucleotide isn't checked, on purpose.
  o1 = findOverlaps(guides, dbSNP, maxgap = 0, minoverlap = 0, ignore.strand = TRUE);
  o2 = findOverlaps(  pams, dbSNP, maxgap = 0, minoverlap = 0, ignore.strand = TRUE);
  o = unique(rbind(as.data.frame(o1), as.data.frame(o2)));
  
  rs  = CharacterList(split(dbSNP[o$subjectHits]$id,  o$queryHits));
  vaf =   NumericList(split(dbSNP[o$subjectHits]$vaf, o$queryHits));
  guides[as.numeric(names(rs))]$dbSNP.rs   = unname(rs);
  guides[as.numeric(names(vaf))]$dbSNP.vaf = unname(vaf);
  
  return(guides);
}
splice.sites = function(refseq.p) {
  f = refseq.p[strand(refseq.p) == "+"];
  r = refseq.p[strand(refseq.p) == "-"];
  
  i = !f$last;
  sd0.f = GRanges(seqnames(f[i]), IRanges(  end(f[i]) +1,   end(f[i]) +0), strand(f[i]));
  sd1.f = shift(sd0.f, +1);
  sd2.f = shift(sd0.f, +2);
  
  i = !f$first;
  sa0.f = GRanges(seqnames(f[i]), IRanges(start(f[i]) -0, start(f[i]) -1), strand(f[i]));
  sa1.f = shift(sa0.f, -1);
  sa2.f = shift(sa0.f, -2);
  
  i = !r$last;
  sd0.r = GRanges(seqnames(r[i]), IRanges(start(r[i]) -0, start(r[i]) -1), strand(r[i]));
  sd1.r = shift(sd0.r, -1);
  sd2.r = shift(sd0.r, -2);
  
  i = !r$first;
  sa0.r = GRanges(seqnames(r[i]), IRanges(  end(r[i]) +1,   end(r[i]) +0), strand(r[i]));
  sa1.r = shift(sa0.r, +1);
  sa2.r = shift(sa0.r, +2);
  
  mcols(sd0.f) = DataFrame("ss" = "sd0");
  mcols(sd1.f) = DataFrame("ss" = "sd1");
  mcols(sd2.f) = DataFrame("ss" = "sd2");
  
  mcols(sd0.r) = DataFrame("ss" = "sd0");
  mcols(sd1.r) = DataFrame("ss" = "sd1")
  mcols(sd2.r) = DataFrame("ss" = "sd2");
  
  mcols(sa0.f) = DataFrame("ss" = "sa0");
  mcols(sa1.f) = DataFrame("ss" = "sa1");
  mcols(sa2.f) = DataFrame("ss" = "sa2");
  
  mcols(sa0.r) = DataFrame("ss" = "sa0");
  mcols(sa1.r) = DataFrame("ss" = "sa1");
  mcols(sa2.r) = DataFrame("ss" = "sa2");
  
  result = sort(c(sd0.f, sd1.f, sd2.f, sd0.r, sd1.r, sd2.r, sa0.f, sa1.f, sa2.f, sa0.r, sa1.r, sa2.r))
  result$ss = factor(result$ss);
  return(result);
}
shift.ss = function(refseq.p) {
  refseq.p = refseq.p +2
  i = ((refseq.p$first & strand(refseq.p) == "+") | 
        (refseq.p$last  & strand(refseq.p) == "-"));
  start(refseq.p[i]) = start(refseq.p[i]) +2;
  i = ((refseq.p$last  & strand(refseq.p) == "+") |
        (refseq.p$first & strand(refseq.p) == "-"));
  end(refseq.p[i]) = end(refseq.p[i]) -2;
  return(refseq.p);
}
legible.columns = function(guides) {
  if(!is.data.frame(guides)) {
    guides = as.data.frame(guides);
  }
  if(is.null(guides$tier)) {
    guides$tier = as.numeric(NA);
    guides$pick = as.character(NA);
  }
  guides$CS = format(guides$CS, big.mark=",", scientific=FALSE, trim = TRUE, drop0trailing = TRUE);
  guides$CS = paste(guides$seqnames, guides$CS, sep = ":");
  guides$isohit = !is.na(guides$isoform_shared.gene.id);
  guides[is.na(guides$ese),"ese"] = FALSE;
  guides$OT.G = paste(guides$G0, guides$G1, guides$G2, guides$G3, guides$G4, sep = ";")
  guides$OT.P = paste(guides$P0, guides$P1, guides$P2, guides$P3, guides$P4, sep = ";")
  #guides$sgname = paste(guides$principal.gene.name, guides$cdd.p.name,
  #                      paste0("T",guides$tier), guides$pick, sep = ".")
  guides$sgname = paste(guides$principal.gene.name,
                        paste0("T",guides$tier), guides$pick, sep = ".")
  guides$homopolymer = paste(guides$PA, guides$PC, guides$PG, guides$PT, sep = ";")
  x = table(guides$seq)
  dup = names(x[x > 1])
  guides$Unique = TRUE
  guides[guides$seq %in% dup,"Unique"] = FALSE
  
  guides = guides[,c("sgname", "seq", "principal.gene.name", "principal.entrez.id", "CS", "tier",
                     "AZ", "provean", "ThreePrime.pct", "isohit",
                     "cdd.p.name", "cdd.p.description", "cdd.p.site", "pfam.name",
                     "homopolymer", "RS", "ese", "ss",
                     "Tm", "HG", "OT.G", "OT.P", "Unique")];
  colnames(guides) = c("Name", "Sequence", "Gene", "Entrez.ID", "Cut.site", "Tier",
                       "Cleavage.Efficacy", "Conservation", "Percent.peptide", "Hits.all.isoforms",
                       "CDD.domain", "CDD.description", "CDD.residue", "Pfam.domain",
                       "Poly:A;C;G;T", "Esp3I", "Splicing.enhancer", "Splice.site",
                       "Tm", "Specificity.score", "Genome.offtargets:0-4mm", 
                       "CDS.offtargets:0-4mm", "Unique");
  return(guides);
}
summary.statistics = function(guides, guides.picked, pdf, SPECIES, refseq.p) {
  pdf(pdf, 8.5, 11);
  par(mfcol = c(3,2));
  #Higher stats
  df = rbind(
   c("Species", SPECIES),
   c("Genes targeted", length(unique(refseq.p$principal.gene.id))),
   c("Genes with any guides, at all", length(unique(guides$principal.gene.id))),
   c("Genes with any guides picked", length(unique(guides.picked$principal.gene.id))),
   c("Genes with six guides picked", length(which(table(guides.picked$principal.gene.id) == 6))),
   c("Genes with any tier 1 guides picked", length(unique(guides.picked[guides.picked$tier == 1]$principal.gene.id))),
   c("Genes with six tier 1 guides picked", length(which(table(guides.picked[guides.picked$tier == 1]$principal.gene.id) == 6))),
   c("Average guides per gene", round(mean(table(guides.picked$principal.gene.name)),2)));
  plot.new();
  grid.table(data.frame(df), rows = NULL, cols = NULL,
             vp=baseViewports()$figure, theme = ttheme_default(base_size = 8));

  #Tier breakdown
  df1 = table(guides.picked[guides.picked$tier %in% 1:6]$tier);
  df2 = table(guides.picked[guides.picked$tier %in% c(7.1,7.2,7.3,7.4,7.5,7.6)]$tier);
  #df2.2 = as.data.frame(table(guides.picked[guides.picked$tier > 7.6 ]$tier));
  #df2.2["5",] = c(NA, NA);
  #df2.2["6",] = c(NA, NA);
  df3 = data.frame(df1, df2);
  colnames(df3) = c("Specific\nTier", "\nGuides", "Non-specific\nTier", "\nGuides");#
                   #, "Ultra-Non-specific\nTier", "\nGuides");
  plot.new();
  grid.table(df3, rows = c("All Features", "Not Domain", "Low Cleavage",
                           "Low Conservation", "C-terminus", "Homopolymers or Esp3I"),
             vp=baseViewports()$figure, theme = ttheme_default(base_size = 8));

  rmats = lapply(list(guides, guides.picked), function(g) {
   return(data.frame("Hits Domain" = !is.na(g$cdd.p.id),
                     "Cuts Well" = g$AZ > 0.4,
                     "Hits Conserved AA" = !is.na(g$provean) & g$provean >= 7,
                     "Avoids C-terminus" = g$ThreePrime.pct < 0.95,
                     "Avoids Homopolymer" = g$PA < POLYN["A"] & g$PC < POLYN["C"] &
                                            g$PG < POLYN["G"] & g$PT < POLYN["T"],
                     "Avoids Restriction Site" = !g$RS,
                     "Hits Target Specifically" = g$HG > 0.50 & g$G0 == 0 & g$G1 == 0,
                     "Avoids Common SNPs" = !any(g$dbSNP.vaf > 0.10)));
  });
  names(rmats) = c("All Guides", "Picked Guides");
  df4 = data.frame(cbind(
          "All Guides" = apply(rmats[["All Guides"]],    2, function(col) {
            n = length(which(col));
            d = length(guides);
            paste(n, "/", d, paste0("(",round(100*n/d, 1),"%)")); }),
          "Picked Guides" = apply(rmats[["Picked Guides"]], 2, function(col) {
            n = length(which(col));
            d = length(guides.picked);
            paste(n, "/", d, paste0("(",round(100*n/d, 1),"%)")); })));
  plot.new();
  grid.table(df4, cols = c("All Guides", "Picked Guides"),
             rows = c("Hits Domain", "Cuts Well", "Hits Conserved AA",
                      "Avoids C-Terminus", "Avoids Homopolymers", "Avoids Restriction Sites", 
                      "Hits Target Specifically", "Avoids Common SNPs"),
             vp=baseViewports()$figure, theme = ttheme_default(base_size = 8));

  my.list = data.frame("col" = c("AZ", "HG", "provean"),
                       "label" = c("Cleavage Efficacy", "Specificity", "Conservation"));
  par(pty = "s");
  
  invisible(apply(my.list, 1, function(i) {
    plot(density(mcols(guides.picked)[,i["col"]], na.rm = T), col = "red",
         main = i["label"], xlab = i["label"]);
    lines(density(mcols(guides)[,i["col"]], na.rm = T), col = "blue");
    legend("topleft", legend = c("All guides", "Picked guides"), col = c("blue", "red"), pch = "-");
    abline(h=0);
    return(NULL);
  }));
  graphics.off();  
}

date();
guides = load.guides(GUIDES.FP);                        #Load Guides
refseq.p = load.refseq.p(REFSEQ.ALL.FP, APPRIS.FP);     #Load Principal-Isoform Refseq
ss = splice.sites(refseq.p);                            #Find splice sites using Refseq
refseq.p = shift.ss(refseq.p);                          #Shift Refseq onto splice sites
refseq.c = load.refseq.c(REFSEQ.ISO.FP, APPRIS.FP);     #Load Isoform-collapsed Refseq
threeprime = load.threeprime(THREEPRIME.FP);            #Load 3'% (derived from APPRIS)
cdd.c = load.cdd.c(CDD.C.FP);                           #Load Isoform-collapsed CDD domains
cdd.p = load.cdd.p(CDD.P.FP);                           #Load Principal isoform CDD domains
cdd.s = load.cdd.s(CDD.S.FP);                           #Load Principal isoform CDD sites
pfam = load.pfam(PFAM.FP);                              #Load Pfam (from UCSC)
uniprot = load.uniprot.2nd(UNIPROT.FP);                 #Load Uniprot (from UCSC)
provean = load.provean(PROVEAN.FP);                     #Load Provean-based Conservation
ese = load.ese(INT3.ESE.FP);                            #Load exonic splicing enhancers
nmd = load.nmd(NMD.FP);                                 #Load NMD predictions
guides = tag.restriction.sites(guides);                 #Tag guides containing RE sites
guides = tag.homopolymers(guides);                      #Tag guides containing homopolymers
date();
guides = tag.snps(guides, refseq.p, GENOME.FP, SNP.FP); #Tag common SNPs for hg38. mm10 = C57BL/6J = "Black 6"
guides$Tm = calcTm(as.character(guides$seq));           #Tag Tm using HELP library
#guides = oligoTm.primer3(guides);                      #Tag Tm using Primer3. Too slow
date();

refseq.p$first = NULL;
refseq.p$last  = NULL;

#merge on data tracks
guides$ranges.o = ranges(guides);
ranges(guides) = IRanges(guides$CS, guides$CS -1);      #Shift guides onto cutsite
for(i in list(refseq.p, refseq.c, cdd.c, cdd.p, cdd.s, pfam, uniprot, nmd, threeprime, ese, ss)) { 
  o = findOverlaps(guides, i, maxgap = 0, minoverlap = 0, ignore.strand = TRUE);
  tmp = mcols(i)[1,,drop=F];
  tmp[1,] = NA;
  mcols(guides)[,colnames(mcols(i))] = tmp;
  mcols(guides[from(o)])[,colnames(mcols(i))] = mcols(i[to(o)])
}

#special merge for provean, to average multiple AA hits (e.g., cuts b/w AAs)
guides$provean = as.numeric(NA);
o = findOverlaps(guides, provean, maxgap = 0, minoverlap = 0);
mcols(guides[from(o)])[,colnames(mcols(provean))] = mcols(provean[to(o)])
dup = unique(from(o)[duplicated(from(o))]);
o = o[from(o) %in% dup];
avg = mean(NumericList(split(provean[to(o)]$provean, from(o))));
guides[as.numeric(names(avg))]$provean = unname(avg);

#Save completed guide database
save(guides, file = paste0(SPECIES, ".guideDB_v", VERSION, ".full.Rdata"))
date();

#Generate Tier list
guides = guides[!is.na(guides$principal.gene.id) & guides$OL == 0];
guides = guides[is.na(guides$ss) | guides$ss == "sd0" | guides$ss == "sa0"]; #Exon edges, but not +1 or +2

#Extremely non-specific guides
guides$tier = 8.5;
guides[guides$HG > 0.05]$tier = 8.4;
guides[guides$HG > 0.10]$tier = 8.3;
guides[guides$HG > 0.15]$tier = 8.2;
guides[guides$HG > 0.20]$tier = 8.1;
guides[guides$HG > 0.25]$tier = 7.6;
#Non-specific guides
guides[guides$tier == 7.6 & guides$HG <= 0.50 & !guides$RS &
       guides$PA < POLYN["A"] & guides$PC < POLYN["C"] &
       guides$PG < POLYN["G"] & guides$PT < POLYN["T"]]$tier = 7.5;
guides[guides$tier == 7.5 & guides$ThreePrime.pct < 0.95]$tier = 7.4;
guides[guides$tier == 7.4 & !is.na(guides$provean) & guides$provean >= 7]$tier = 7.3;
guides[guides$tier == 7.3 & guides$AZ > 0.4]$tier = 7.2;
guides[guides$tier == 7.2 & !is.na(guides$cdd.p.id)]$tier = 7.1;
#Specific guides
guides[guides$HG > 0.50 &
       guides$G0 == 0 &
       guides$G1 == 0 &
       !any(guides$dbSNP.vaf > 0.10)]$tier = 6;
guides[guides$tier == 6 & !guides$RS &
       guides$PA < POLYN["A"] & guides$PC < POLYN["C"] &
       guides$PG < POLYN["G"] & guides$PT < POLYN["T"]]$tier = 5;
guides[guides$tier == 5 & guides$ThreePrime.pct < 0.95]$tier = 4;
guides[guides$tier == 4 & !is.na(guides$provean) & guides$provean >= 7]$tier = 3;
guides[guides$tier == 3 & guides$AZ > 0.4]$tier = 2;
guides[guides$tier == 2 & !is.na(guides$cdd.p.id)]$tier = 1;
guides = guides[order(guides$tier, -guides$AZ)];
guide.gene = split(guides, guides$principal.gene.id);

#Pick Guides
#precomputed pick orders, based on validity of minimum spacing criterion
p = lapply(1:MAX.GUIDE.PERMUT, function(i) {
      if(i <= GUIDES.PER.GENE) { return(NULL); }
      p = t(combn(i,GUIDES.PER.GENE));
      return(p[order(rowMeans(p)),]);
    });
date();
r = lapply(guide.gene, function(gd) {
      if(length(gd) <= GUIDES.PER.GENE) { #fewer than 6 sgs/gene
        gd$pick = paste(1:length(gd), length(gd), sep = "/");
        return(gd);
      }
      force = gd$tier < gd$tier[GUIDES.PER.GENE +1]; #force include sgs above #7 in a higher tier
      d = as.matrix(dist(gd$CS)) < MIN.GUIDE.DIST; #overlap matrix
      diag(d) = FALSE; #remove self overlaps
      d[force, force] = FALSE; #ignore overlaps of forced sgs w/ each other
      p.id = min(length(gd), MAX.GUIDE.PERMUT); #id of pre-computed permutations
      v = apply(p[[p.id]], 1, function(set) { !any(d[set,set]); }) #find permitted sets
      r = NA;
      if(any(v)) {
        if(any(force)) { #best permitted set including all forced sgs
          n.force = max(which(force));
          r = gd[p[[p.id]][which(v & p[[p.id]][,n.force] == n.force)[1],]];
        } else {
          r = gd[p[[p.id]][which(v)[1],]];
        }
      } else { #no non-overlapping set - return top 6
        r = gd[1:GUIDES.PER.GENE];
      }
      r$pick = paste(1:length(r), length(r), sep = "/");
      return(r);
    });
date();
guides.picked = do.call(c, unname(r));
my.guides = legible.columns(as.data.frame(guides.picked));
summary.statistics(guides, guides.picked, paste0(SPECIES, ".guideDB_v", VERSION, ".pdf"), SPECIES, refseq.p)
save(guides.picked, file = paste0(SPECIES, ".guideDB_v", VERSION, ".Rdata"));
write.xlsx(my.guides, file = paste0(SPECIES, ".guideDB_v", VERSION, ".xlsx"))
#write.csv(my.guides, file = paste0(SPECIES, ".guideDB_v", VERSION, ".csv"), row.names = FALSE);
q();

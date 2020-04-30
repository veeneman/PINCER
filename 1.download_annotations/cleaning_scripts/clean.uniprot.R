#!/usr/bin/env Rscript
#./clean.uniprot.R
library("stringr", "stringi", "rtracklayer", "GenomicAlignments", "GenomicFeatures");

filter.pri = function(g) {
  g = g[seqnames(g) %in% paste0("chr",c(1:22,"X","Y","M")),];
  seqlevels(g) = seqlevels(g)[seqlevels(g) %in% paste0("chr",c(1:22,"X","Y","M"))];
  return(g);
}

data = list();
x = c("unipLocSignal", "unipChain",
      "unipLocCytopl", "unipLocTransMemb", "unipLocExtra",
      "unipDisulfBond",
      "unipModif",
      "unipMut",
      "unipOther",
      "unipDomain",
      "unipStruct",
      "unipRepeat");
for(i in x) { data[[i]] = import(paste0(i,".bb.bed")); }

x = data[["unipLocSignal"]];
x$desc = x$NA.3
x = x[,c("name", "desc", "blocks")];
y = data[["unipChain"]];
y$desc = "NA";
y$name = "Mature.protein";
y = y[,c("name", "desc", "blocks")];
uniprot.processing = filter.pri(c(x, y));
uniprot.processing$name = paste0("Processing:", uniprot.processing$name);

x = data[["unipLocCytopl"]];
x$desc = "NA";
y = data[["unipLocTransMemb"]];
y$desc = "NA";
z = data[["unipLocExtra"]];
z$desc = "NA";
uniprot.localization = filter.pri(c(x, y, z)[,c("name","desc","blocks")]);
uniprot.localization$name = paste0("Localization:", uniprot.localization$name);

x = data[["unipDisulfBond"]];
x$name = "Disulfide.bond";
x$desc = x$NA.3;
uniprot.disulfide = filter.pri(x[,c("name","desc","blocks")]);

x = data[["unipModif"]];
x$desc = x$NA.3;
key = c("acetyl", "glyco", "hydroxy", "Lipidation", "lipo", "methyl", "modif", "nitro", "phos", "sulfo")
value = c("Acetylation", "Glycosylation", "Hydroxylation", "Lipidation",
          "Lipoylation", "Methylation", "Other", "Nitrosylation", "Phosphorylation", "Sulfinylation")
lookup = setNames(value, key);
x$name = lookup[x$name];
uniprot.ptm = filter.pri(x[,c("name", "desc", "blocks")]);
uniprot.ptm$name = paste0("PTM:", uniprot.ptm$name);

x = data[["unipMut"]];
x$desc = paste(x$name, x$NA.4, x$NA.6, sep = ",");
x$name = x$NA.1;
x[x$name == "Naturally.occurring.sequence.variant"]$name = "Variant:Polymorphism";
x[x$name == "Experimental.mutation.of.amino.acids"]$name = "Variant:Mutation";
uniprot.mut = filter.pri(x[,c("name", "desc", "blocks")]);

x = data[["unipOther"]];
x$desc = paste(x$NA.1, x$name, x$NA.3, sep = ",");
x$name = "Other";
uniprot.other = filter.pri(x[,c("name", "desc", "blocks")]);

x = data[["unipDomain"]];
r = (x$NA.1 == "zinc.finger.region");
x[r]$NA.3 = paste0("zinc.finger,", x[r]$NA.3);
r = (x$NA.1 == "topological.domain");
x[r]$NA.3 = paste0("topological,", x[r]$NA.3);
x$name = "Domain";
x$desc = x$NA.3;
uniprot.domain = filter.pri(x[,c("name","desc","blocks")]);

x = data[["unipStruct"]];
x$name = paste0("Secondary.structure:", x$NA.1);
x$desc = "NA";
uniprot.structure = filter.pri(x[,c("name","desc","blocks")]);

x = data[["unipRepeat"]];
x$name = "Repeat";
x$desc = x$NA.3;
uniprot.repeat = filter.pri(x[,c("name","desc","blocks")]);

#concatenate uniprot tracks to a single track
uniprot.o1 = c(uniprot.processing,
              uniprot.localization,
              uniprot.disulfide,
              uniprot.ptm,
              uniprot.mut,
              uniprot.other,
              uniprot.domain,
              uniprot.structure,
              uniprot.repeat);
uniprot.o1[uniprot.o1$desc == "NA"]$desc = NA;

#expand out gapped ranges
tmp = as.data.frame(blocks(uniprot.o1), use.outer.mcols = TRUE);
uniprot.o2 = GRanges(tmp$seqnames, IRanges(tmp$start, tmp$end), tmp$strand,
                     DataFrame(tmp[,8:ncol(tmp)]));
#find unique blocks
uniprot = disjoin(uniprot.o2);
#identify original uniprot ranges hitting unique blocks (one->many)
o = findOverlaps(uniprot, uniprot.o2, maxgap = 0, minoverlap = 1);

tmp1 = split(mcols(uniprot.o2[to(o)])$name, from(o));
tmp2 = split(mcols(uniprot.o2[to(o)])$desc, from(o));
name.new = sapply(tmp1, function(j) { paste(sort(unique(j)), collapse = "; "); });
desc.new = sapply(tmp2, function(j) { paste(sort(unique(j)), collapse = "; "); });
mcols(uniprot)[as.numeric(names(name.new)), "name"] = unname(name.new);
mcols(uniprot)[as.numeric(names(desc.new)), "desc"] = unname(desc.new);
save(uniprot, file = "uniprot.Rdata");

#debug
#sapply(colnames(mcols(a)), function(i) { if(i == "thick" | i == "blocks") { return(NA); } else { return(head(sort(table(mcols(a)[,i]), decreasing = T))); } })
#sapply(colnames(mcols(x)), function(i) { if(i == "thick" | i == "blocks") { return(NA); } else { return(head(sort(table(mcols(x)[,i]), decreasing = T))); } })

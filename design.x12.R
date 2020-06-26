#!/usr/bin/env Rscript
argv = commandArgs(trailingOnly = TRUE);
SPECIES            = argv[1]; #"mm10"; # "hg38"; #
GUIDES.PER.GENE    = 12; #6; #
MAX.GUIDE.PERMUT   = 21; #15; #heuristic, doesn't affect result unless it's too low
print(MAX.GUIDE.PERMUT);
MIN.GUIDE.DIST     = 3;
VERSION            = "1.3c";
POLYN              = setNames(c(5, 5, 5, 3), c("A", "C", "G", "T")); #minimum dis-allowed homopolymer length

library("stringr")
library("stringi")
library("grid")
library("gridBase")
library("gridExtra")
library("rtracklayer")
library("GenomicAlignments")
library("GenomicFeatures")
library("HELP")
library("openxlsx")
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

load(paste0("../../submission.04302020/to.github/database.", SPECIES, ".Rdata"));

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
save(guides.picked, file = paste0(SPECIES, ".guideDB_v", VERSION, ".Rdata"));
write.xlsx(my.guides, file = paste0(SPECIES, ".guideDB_v", VERSION, ".xlsx"));
#write.csv(my.guides, file = paste0(SPECIES, ".guideDB_v", VERSION, ".csv"), row.names = FALSE);
q();



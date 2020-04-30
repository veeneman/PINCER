#!/usr/bin/env Rscript
#./process.conservation.R cdd.aa.gff2.gz refseq.to.hg38.gff3.gz output.gff3 output.Rdata
library("stringr", "stringi", "rtracklayer", "GenomicAlignments", "GenomicFeatures");

die = function(msg) {
  write(paste("Error:",msg), stderr());
  q(status = 1);
}
aa.to.cDNA = function(aa) {
  #aa.nt to cDNA.nt (multiply by three and deal with 1-indexing)
  ranges(aa) = IRanges((3 * (start(aa) -1)) +1,
                       (3 * (  end(aa) -0)) +0);
  
  map = aa[aa$type %in% "CDS"][,"coded_by"] #use CDS rows as the map
  rc = as.character(seqnames(map[grepl("complement",map$coded_by)])); #this applies to One Gene.  Thx Refseq
  map$coded_by = str_replace_all(map$coded_by,"^complement\\(|^join\\(|\\)$","")
  map = map[!grepl("^XM_", map$coded_by)] #drop predicted proteins, matches my refseq processing
  map.1 = map[!grepl(",", map$coded_by)] #aa->cDNA w/o gaps (the easy case)
  nms = do.call(rbind, strsplit(map.1$coded_by, ":|\\.\\.")) #NM_001:1..100 #human had no whitespace here
  map.1$NM = nms[,1]
  map.1$offset = as.numeric(nms[,2]) -1
  map.m = map[grepl(",", map$coded_by)] #aa->cDNA with gaps (the hard case)
  map.m = unlist(GRangesList(lapply(map.m, function(m) { 
                #split as before, but by commas and whitespace first
                tmp = do.call(rbind, strsplit(unlist(strsplit(m$coded_by, ",\\s*")), ":|\\.\\."));
                cDNA.s = as.numeric(tmp[,2]) -1;            #cDNA starts
                cDNA.e = as.numeric(tmp[,3]);               #cDNA ends
                widths = cDNA.e - cDNA.s;                   #ungapped chunk lengths
                aa.s = c(0, cumsum(widths[-nrow(tmp)])) +1; #aa starts
                aa.e = cumsum(widths);                      #aa ends
                #resolving the offsets is easy here because it's all + strand to + strand
                offsets = cumsum(c(cDNA.s[1], cDNA.s[-1] - cDNA.e[-nrow(tmp)]));
                mcols = DataFrame("NM" = tmp[,1], "offset" = offsets);
                return(unname(GRanges(seqnames(m), IRanges(aa.s, aa.e), "+", mcols)));
              })));
  map = c(map.1[,-1],map.m); #combine ungapped and gap-resolved maps
  #this is the actual mapping step - it intersects cDNA coords on cDNA sequences
  #  i.e. - using NM_* as the "chromosomes".
  overlap = findOverlaps(aa, map, ignore.strand=T); #find any overlaps (cDNA coords)
  cDNA = pintersect(aa[from(overlap)], map[to(overlap)], ignore.strand = T) #intersect overlaps (cDNA coords)
  mcols(cDNA) = cbind(mcols(aa[from(overlap)]), mcols(map[to(overlap)])) #tack on original metadata
  cDNA = shift(cDNA, cDNA$offset) #right-shift everything by array of offsets. only possible b/c of + strand
  keep.cols = c("NP", "NM", "type", "db_xref", "note") #"name", 
  cDNA = GRanges(cDNA$NM, ranges(cDNA), strand(cDNA), mcols(cDNA)[,keep.cols]) #drop junk
  strand(cDNA[cDNA$NP %in% rc]) = "-"; #fix strand for single negative strand mitochondrial gene
  return(cDNA);
}
load.cDNA.genome.map = function(fp) {
  map.o = import(fp);
  
  #filter for primary assembly.
  #note: Refseq mapped NM to the genome 1:2 for the Pseudoautosomal region (X & Y).
  #      I used both mappings here, but just use chrX for gene targeting.
  #      The "ID" column here could support logic re:PAR.
  map.o = map.o[seqnames(map.o) %in% paste0("chr",c(1:22,"X","Y","M")),];
  seqlevels(map.o) = seqlevels(map.o)[seqlevels(map.o) %in% paste0("chr",c(1:22,"X","Y","M"))]
  
  #deal with singleton alignment mistake in mouse
  id = which(map.o$Target == "NM_011093.1 993 1301 +" & map.o$Gap == "I1 M308")
  if(length(id) > 0) {
    map.o[id]$Target = "NM_011093.1 994 1301 +";
    map.o[id]$Gap = NA;
  }
  
  #split apart target string, eg: "NM_001 1 200", and restore as a GRange. Gap field encodes gaps
  tmp = do.call(rbind, strsplit(map.o$Target," "))
  map.o = GRanges(seqnames(map.o), ranges(map.o), strand(map.o),
      DataFrame("Gap" = map.o$Gap,
          "Target" = tmp[,1],
          "cDNA.s" = as.numeric(tmp[,2]),
          "cDNA.e" = as.numeric(tmp[,3])));
  
  #ungapped ranges
  map.1 = map.o[is.na(map.o$Gap),c("Target","cDNA.s","cDNA.e")]
  
  #split apart gapped ranges
  map.2 = unlist(GRangesList(lapply(map.o[!is.na(map.o$Gap)], function(x) {
                #split apart gapped CIGAR
                CIGAR = do.call(rbind, strsplit(unlist(strsplit(x$Gap,"\\s")),"(?<=\\D)", perl = T))
                m = as.numeric(CIGAR[CIGAR[,1] == "M",2,drop=F]) #chunk lengths
                offsets.g = m[-length(m)]; #genome offsets
                offsets.m = m[-length(m)]; #cDNA offsets (mRNA)
                ops = CIGAR[CIGAR[,1] != "M",,drop=F] #gaps
                offsets.g[ops[,1] == "D"] = offsets.g[ops[,1] == "D"] + as.numeric(ops[ops[,1] == "D",2]) #offsets
                offsets.m[ops[,1] == "I"] = offsets.m[ops[,1] == "I"] + as.numeric(ops[ops[,1] == "I",2]) #offsets
                starts.m = cumsum(c(x$cDNA.s -1, offsets.m)) #offsets
                if(as.character(strand(x)) == "-") {
                  starts.g = c(end(x), end(x) - cumsum(offsets.g)) +1 - m; #subtract cDNA coords from exon ends
                  ends.g   = c(end(x), end(x) - cumsum(offsets.g));        #subtract cDNA coords from exon ends
                } else {
                  starts.g = cumsum(c(start(x) -1, offsets.g)) +1; #add cDNA coords to exon starts, so re-adding works
                  ends.g = cumsum(c(start(x) -1, offsets.g)) +m;   #add cDNA coords to exon starts, so re-adding works
                }
                result = GRanges(seqnames(x), IRanges(starts.g, ends.g), strand(x),
                    DataFrame("Target" = x$Target,
                        "cDNA.s" = starts.m +1,
                        "cDNA.e" = starts.m +m))
                return(result);
              })));
  map.3 = c(map.1,map.2); #recombine gapped and ungapped ranges into single GRanges
  map.3 = map.3[order(map.3$Target, map.3$cDNA.s)]
  #set up per-range offset, which is adjusted for stranded cDNA positions
  offset = start(map.3) - map.3$cDNA.s;
  neg.g = as.character(strand(map.3)) == "-";
  offset[neg.g] = end(map.3)[neg.g] + map.3$cDNA.s[neg.g];
  map = GRanges(map.3$Target, IRanges(map.3$cDNA.s, map.3$cDNA.e), "+",
      DataFrame("chr" = as.character(seqnames(map.3)),
          "offset" = offset,
          "chr.s" = as.character(strand(map.3))))
  #manually add chrM: NC_012920.1 is human, NC_005089.1 is mouse
  map = suppressWarnings(c(map, GRanges(c("NC_012920.1","NC_005089.1"),
              IRanges(c(1,1), c(100000,100000)),
              c("+","+"),
              DataFrame("chr" = c("chrM","chrM"),
                  "offset" = c(1,1),
                  "chr.s" = c("+","+")))))
  #manually remove annotation errors related to human gene "NAIP".  mouse may also have some.eg:
  #cov = coverage(map.cDNA_genome)
  #result = sapply(cov, function(i) { return(any(i > 1)); })
  #multimapped = names(result[which(result)])
  #multimapped[!grepl("^NR",multimapped)]
  #x = data.frame(map.cDNA_genome[seqnames(map.cDNA_genome) %in% multimapped])
  #unique(mcols(cdd.cDNA[cdd.cDNA$NM %in% c("NM_001346870.1","NM_004536.2", "NM_022892.1")][,c("NM","NP")]))
  blacklist.hs = c("NM_001346870.1","NM_004536.2", "NM_022892.1");
  map = map[!(seqnames(map) %in% blacklist.hs & map$chr.s == "+")];
  blacklist.mm = c("NM_001033450.4", "NM_001301745.1");
  map = map[!(seqnames(map) %in% blacklist.mm & map$offset < 173800000)];
  return(map);
}
cDNA.to.genome = function(cdd.cDNA, map) {
  #drop NMs not in the cDNA (avoids warnings)
  map = map[as.character(seqnames(map)) %in% seqlevels(cdd.cDNA)]
  seqlevels(map) = unique(as.character(seqnames(map)))
  
  rc = unique(cdd.cDNA[as.character(strand(cdd.cDNA)) %in% "-"]$NP); #single mt gene
  strand(cdd.cDNA) = "+"; #single mt gene
  hits = findOverlaps(cdd.cDNA, map); #actually map cDNA to genome, using NP seqnames
  result.o = pintersect(cdd.cDNA[from(hits)], map[to(hits)]);
  mcols(result.o) = cbind(mcols(cdd.cDNA[from(hits)]), mcols(map[to(hits)]))
  result = GRanges(result.o$chr,
      IRanges(result.o$offset + start(result.o),
          result.o$offset + end(result.o)),
      result.o$chr.s,
      DataFrame(mcols(result.o)[,c("NP","NM","type","db_xref","note")],
          "cDNA.s" = start(result.o),
          "cDNA.e" = end(result.o)))
  neg.g = result.o$chr.s == "-";
  ranges(result[neg.g]) = IRanges(result.o[neg.g]$offset - end(result.o[neg.g]),
      result.o[neg.g]$offset - start(result.o[neg.g]))
  result = sortSeqlevels(result)
  result = result[order(as.character(seqnames(result)), start(result), end(result))]
  strand(result[result$NP %in% rc]) = "-";
  return(result);
}

#argv = c("input.aa.bed", "output.genome.bedgraph", "cdd.aa.gff2.gz", "refseq.to.genome.gff3.gz")
argv = commandArgs(trailingOnly = TRUE);
input.fp      = argv[1];
output.fp     = argv[2];
CDD.AA.fp     = argv[3];
Refseq.map.fp = argv[4];
if(!file.exists(input.fp)) { die(paste0("missing:", input.fp)); }
if(!file.exists(CDD.AA.fp)) { die(paste0("missing:", CDD.AA.fp)); }
if(!file.exists(Refseq.map.fp)) { die(paste0("missing:", Refseq.map.fp)); }

cdd.aa = import(CDD.AA.fp);
keep.cols= c("coded_by", "type", "db_xref", "region_name", "group", "site_type");
g = import(input.fp)
g = GRanges(str_replace(seqnames(g), "\\.\\d+$", ""), ranges(g), "*", DataFrame("note" = g$name));
mcols(g)[,keep.cols] = "NA";
g = c(g, cdd.aa[cdd.aa$type == "CDS"][,c("note",keep.cols)]);
g$NP = as.character(seqnames(g))

g.cDNA = aa.to.cDNA(g);

map.cDNA_genome = load.cDNA.genome.map(Refseq.map.fp);
g.genome = cDNA.to.genome(g.cDNA, map.cDNA_genome);
g.genome = g.genome[g.genome$type != "CDS"]

#FORK - MANUALLY SAVE #13 and #17 HERE
mcols(g.genome) = DataFrame("score" = as.numeric(g.genome$note))
strand(g.genome) = "*"
export(g.genome, output.fp, format = "bedgraph")

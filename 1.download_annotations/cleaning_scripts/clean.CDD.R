#!/usr/bin/env Rscript
#./clean.CDD.R cdd.aa.gff2.gz refseq.to.hg38.gff3.gz output.gff3 output.Rdata
library("stringr", "stringi", "rtracklayer", "GenomicAlignments", "GenomicFeatures");

load.aa = function(fp) {
  #load and remove
  # - "source" rows, which contain no useful information
  # - "Bond" rows, which are too sparse (28)
  # - "misc_feature" rows, which are too sparse (23)
  # - "proprotein" rows, which are too sparse (302), and redundant with "Protein"s
  # - "mat_peptide" rows, which are abundant (6204), but reflect post-translational processing
  #   these could be useful later, but I can't think how to use them now
  aa = import(fp);
  aa$type = as.character(aa$type);
  aa$name = as.character(NA);
  aa = aa[!aa$type %in% c("source", "Bond", "misc_feature", "proprotein", "mat_peptide")];
  
  #pull full gene names from "Protein" rows, assign them to "CDS" rows by NP_#
  rp = aa$type == "Protein";
  rc = aa$type == "CDS";
  tmp = aa[rp]$note
  names(tmp) = as.character(seqnames(aa[rp]))
  aa[rc]$note = tmp[as.character(seqnames(aa[rc]))]
  aa = aa[!rp];
  rm(rp, rc, tmp);
  
  #Reassign metadata for CDS rows to name and note fields
  r = aa$type == "CDS";
  aa[r]$name = aa[r]$gene;
  aa[r]$gene_synonym = str_replace(aa[r]$gene_synonym, ";", ",")
  aa[r]$note = paste0("gene_synonym=", aa[r]$gene_synonym, "; ", aa[r]$note)
  rm(r);
  
  #convert transit and signal peptide features into Regions
  r = aa$type %in% c("transit_peptide", "sig_peptide");
  aa[r]$region_name = aa[r]$type;
  aa[r]$type = "Region";
  rm(r);
  
  #Rename "Region"s to be unique within an NP_#, and assign to name field
  r = aa$type == "Region";
  tmp = make.unique(paste(as.character(seqnames(aa[r])), aa[r]$region_name));
  aa[r]$name = str_match(tmp, "^\\S+ (.+)$")[,2];
  rm(r, tmp)
  
  #Rename "Site"s to be unique within an NP_#, and assign to name field
  r = aa$type == "Site" & !is.na(aa$site_type);
  tmp = make.unique(paste(as.character(seqnames(aa[r])), aa[r]$site_type));
  tmp = str_match(tmp, "^\\S+ (.+)$")[,2];
  aa[r]$name = tmp;
  rm(r, tmp);
  
  #Fill out blank Site rows
  cols = c("name", "db_xref", "note");
  id = paste(as.character(seqnames(aa)), aa$group);
  r = aa$type == "Site" & !is.na(aa$site_type) & !is.na(aa$group);
  tmp = mcols(aa[r])[,cols];
  rownames(tmp) = id[r];
  r = aa$type == "Site" & !is.na(aa$group);
  mcols(aa[r])[,cols] = tmp[id[r],]
  rm(id, r, cols, tmp);
  
  #restructure columns
  keep.cols = c("coded_by", "type", "name", "db_xref", "note")
  mcols(aa) = mcols(aa)[,keep.cols]
  aa$NP = as.character(seqnames(aa))  
  return(aa);
}
aa.to.cDNA = function(aa) {
  #aa.nt to cDNA.nt (multiply by three and deal with 1-indexing)
  ranges(aa) = IRanges((3 * (start(aa) -1)) +1,
                       (3 * (  end(aa) -0)) +0)
  
  map = aa[aa$type %in% "CDS"][,1] #use CDS rows as the map
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
  overlap = findOverlaps(aa, map, ignore.strand=T, maxgap = 0, minoverlap = 1); #find any overlaps (cDNA coords)
  cDNA = pintersect(aa[from(overlap)], map[to(overlap)], ignore.strand = T) #intersect overlaps (cDNA coords)
  mcols(cDNA) = cbind(mcols(aa[from(overlap)]), mcols(map[to(overlap)])) #tack on original metadata
  cDNA = shift(cDNA, cDNA$offset) #right-shift everything by array of offsets. only possible b/c of + strand
  keep.cols = c("NP", "NM", "type", "name", "db_xref", "note")
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
  hits = findOverlaps(cdd.cDNA, map, maxgap = 0, minoverlap = 1); #actually map cDNA to genome, using NP seqnames
  result.o = pintersect(cdd.cDNA[from(hits)], map[to(hits)]);
  mcols(result.o) = cbind(mcols(cdd.cDNA[from(hits)]), mcols(map[to(hits)]))
  result = GRanges(result.o$chr,
                   IRanges(result.o$offset + start(result.o),
                           result.o$offset + end(result.o)),
                   result.o$chr.s,
                   DataFrame(mcols(result.o)[,c("NP","NM","type","name","db_xref","note")]))
  neg.g = result.o$chr.s == "-";
  ranges(result[neg.g]) = IRanges(result.o[neg.g]$offset - end(result.o[neg.g]),
                                  result.o[neg.g]$offset - start(result.o[neg.g]))
  result = sortSeqlevels(result)
  result = result[order(as.character(seqnames(result)), start(result), end(result))]
  strand(result[result$NP %in% rc]) = "-";
  
  #propagate gene names to domains
  gene.lookup = setNames(result[result$type == "CDS"]$name, result[result$type == "CDS"]$NP);
  result$gene = unname(gene.lookup[result$NP]);
  
  return(result);
}
format.for.browser = function(cdd, cols) {
  #give ranges UIDs by feature type (cds, region, site), by the specific site, chromosome, and NP
  cdd$ID = do.call(paste, args = c(data.frame(cdd)[,cols], sep = "."));
  genes = cdd[!duplicated(cdd$ID)]; #pull single rows for each feature (by UID)
  cdd$type = "CDS";    #original features become child "CDS" for visualization
  cdd$Parent = cdd$ID; #parents become the "genes"
  genes$type = "gene"; #parent features
  genes$Parent = NA;
  mcols(genes)[,colnames(mcols(genes)) %in% c("cDNA.s", "cDNA.e")] = NA; #blank these
  start(genes) = sapply(split(start(cdd), cdd$Parent), min)[genes$ID]; #set gene ranges to max range
  end(genes)   = sapply(split(  end(cdd), cdd$Parent), max)[genes$ID]; #set gene ranges to max range
  cdd = c(genes, cdd); #concatenate them back together
  cdd$ID = make.unique(cdd$ID); #this increments all the child UIDs, b/c they're after the parents in the list
  cdd = sortSeqlevels(cdd) #coord sort
  cdd = cdd[order(as.character(seqnames(cdd)), start(cdd), end(cdd))] #coord sort
  return(cdd);
}
split.domains = function(g, cols) {
  #Splits a GRanges (g) by pasted columns (cols),
  # into a GRangesList, and preserves metadata
  ids = do.call(paste, mcols(g)[,cols]);
  gl = split(g[,NULL], ids);
  m = unique(mcols(g));
  rownames(m) = do.call(paste, m[,cols]);
  mcols(gl) = m[names(gl),];
  return(gl);
}
rename.chrs = function(gl, paste.cols, split.cols) {
  df = as.data.frame(gl, use.outer.mcols = TRUE); #unlist the GRangesList, with mcols
  g = GRanges(do.call(paste, df[,paste.cols]),
              IRanges(df$start, df$end),
              df$strand,
              DataFrame(df[,8:ncol(df)])); #create GRanges, with CDD id in the seqnames now
  gl = split.domains(g, split.cols); #split GRanges
}
collapse.metadata = function(gl , o) {
  meta = split(as.matrix(mcols(gl)), 1:length(gl)); #list of vectors of metadata, split by row number
  x = unique(to(o)); #row numbers we will assign to
  y = split(c(x, from(o)), c(x, to(o))); #split "from" by "to", and prepend one copy of "to"
  #For each row number vector,
  # Bind together columns into vectors, and for each of those,
  # Take the unique values and paste them together.
  # Finally, rbind the list of results into a matrix, and wrap it in a DataFrame.
  z = DataFrame(do.call(rbind, lapply(y, function(i) {
        apply(do.call(cbind, meta[i]), 1, function(r) {
          paste(unique(unlist(strsplit(r, ";"))), collapse = ";"); }); })));
  colnames(z) = colnames(mcols(gl)); #reassign column names
  mcols(gl[as.numeric(rownames(z))]) = z; #"to" was used as the split list's names, and rbind preserves them
  return(gl);
}
drop.identical.or.embedded = function(gl, retain.metadata) {
  gl = rev(gl); #Combining reverse order with redundant hit filter means keeping the first instance
  o = findOverlaps(gl, type = "within", maxgap = 0, minoverlap = 1); #implies either "within" OR "equal"
  o = o[!isSelfHit(o) & !isRedundantHit(o)];
  if(length(o) == 0) { return(rev(gl)); }
  if(retain.metadata) { gl = collapse.metadata(gl, o); }
  return(rev(gl[-unique(from(o))]));
}
collapse.domains = function(cdd) {
  #filter for Regions (domains), with CDD ids
  # this excludes domains from UniprotKB (b/c of overlap issues),
  # and unannotated sources (eg: signal peptides, also b/c of overlap issues),
  # and Site features, which seem useful but harder to use
  #cdd = cdd[cdd$type == "Region" & grepl("^CDD\\:", cdd$db_xref)];
  
  #Split domain GRanges into individual domains, and reposition onto same gene for filtering
  cdd$chr.o = as.character(seqnames(cdd));
  cdd.split = split.domains(cdd, c("NP", "name", "chr.o"));
  cdd.split = rename.chrs(cdd.split, c("gene", "chr.o"), c("NP", "name", "chr.o")); #166382
  cdd.split = cdd.split[order(-sum(width(cdd.split)),mcols(cdd.split)$NP)]; #sort by size, then NP
  cdd.split = drop.identical.or.embedded(cdd.split, FALSE); #47419
  
  #Merge overlapping domains with same CDD ID
  #Reposition onto CDD ID for merging
  cdd.split = rename.chrs(cdd.split, c("gene", "chr.o", "db_xref"), c("NP", "name", "chr.o"));
  tmp = cdd.split[order(-sum(width(cdd.split)),mcols(cdd.split)$NP)]; #sort by size, then NP
  cdd.split = GRangesList();
  while(TRUE) {
    o = findOverlaps(tmp, maxgap = 0, minoverlap = 1); #domains touching each other
    o = o[from(o) < to(o)]; #remove self hits and reciprocal hits
    if(length(o) == 0) { cdd.split = c(cdd.split, tmp); break; }
    cdd.split = c(cdd.split, tmp[-unique(c(from(o), to(o)))]); #domains unmerged in next iteration
    tmp = union(tmp[from(o)], tmp[to(o)]); #pairwise union
    #You need to run this step because:
    # 3 mutually overlapping but not embedded domains will union to 3 infinitely, and
    # 4 mutually overlapping but not embedded domains will union to an infinite number of overlaps
    #I considered setting keep metadata to TRUE, but then names might stop being unique,
    # which will break the re-split later. I'm leaving it as FALSE unless there's feedback otherwise
    tmp = drop.identical.or.embedded(tmp, FALSE); #drop redundant domains (WITHIN CDD ID)
  } #44873
  
  #Merge domains overlapping by 80% w/o same CDD ID
  #Reposition back off of CDD IDs
  cdd.split = rename.chrs(cdd.split, c("gene", "chr.o"), c("NP", "name", "chr.o")); #44873
  tmp = cdd.split[order(-sum(width(cdd.split)), mcols(cdd.split)$NP)]; #sort by size, then NP
  mcols(tmp)$name = str_replace(mcols(tmp)$name, "\\.\\d+$", ""); #De-Unique Domain names
  tmp2 = drop.identical.or.embedded(tmp, TRUE); #special merge for embedded or identical #44756
  cdd.split = GRangesList();
  while(TRUE) {
    o = findOverlaps(tmp, maxgap = 0, minoverlap = 1); #domains touching each other
    o = o[from(o) < to(o)]; #remove self hits and reciprocal hits
    w = sum(width(intersect(tmp[from(o)], tmp[to(o)]))); #width of hits' intersections
    o = o[w / sum(width(tmp[from(o)])) > 0.80 |
          w / sum(width(tmp[to(o)])) > 0.80]; #exclude hits w/o 80% overlap
    if(length(o) == 0) { cdd.split = c(cdd.split, tmp); break; }
    cdd.split = c(cdd.split, tmp[-unique(c(from(o), to(o)))]); #domains unmerged in next iteration
    
    tmp = collapse.metadata(tmp, o);
    tmp = union(tmp[to(o)], tmp[from(o)]); #pairwise union, keeps "to" metadata, same as collapse.metadata
    tmp = drop.identical.or.embedded(tmp, TRUE); #special merge for embedded or identical
  } #41587	
  
  #Re-Unique domain names
  ids = paste(mcols(cdd.split)$gene, mcols(cdd.split)$name);
  mcols(cdd.split)$domain.name = str_match(make.unique(ids), "^\\S+\\s(.+)$")[,2];
  mcols(cdd.split)$domain.id = 1:length(cdd.split);
  cdd = as.data.frame(cdd.split, use.outer.mcols = TRUE);
  keep.cols = c("NP", "NM", "type","db_xref", "note", "gene", "domain.name", "domain.id");
  cdd = GRanges(cdd$chr.o, IRanges(cdd$start, cdd$end), cdd$strand, DataFrame(cdd[,keep.cols]));
  cdd = sortSeqlevels(cdd); #coord sort
  cdd = cdd[order(as.character(seqnames(cdd)), start(cdd), end(cdd))]; #coord sort
  return(cdd);
}
collapse.sites = function(cdd) {
  #adjust onto NP ids - don't want to collapse sites from different proteins
  cdd$chr.o = as.character(seqnames(cdd));
  cdd = GRanges(paste(as.character(seqnames(cdd)), cdd$NP),
                ranges(cdd), strand(cdd), mcols(cdd));
  
  #actual collapse step
  result = disjoin(cdd);
  o = findOverlaps(result, cdd, maxgap = 0, minoverlap = 1); #FROM == COLLAPSED ; TO == ORIGINAL SITES
  
  #blunt assign an arbitrary matching site to collapsed bins
  # this is correct for 1:1 sites, and all of the columns Except: name, db_xref, and note
  o.uniq = o[!duplicated(from(o))]; #single result hit
  if(!identical(from(o.uniq), 1:length(result))) { die("error in site collapse"); } #check correct order
  mcols(result) = mcols(cdd[to(o.uniq)]); #assign from original cdd sites
  
  #fix collapsed bins which match multiple input sites
  o.dup = o[from(o) %in% unique(from(o)[duplicated(from(o))])];
  
  tmp1 = split(mcols(cdd[to(o.dup)])$name,    from(o.dup));
  tmp2 = split(mcols(cdd[to(o.dup)])$db_xref, from(o.dup));
  tmp3 = split(mcols(cdd[to(o.dup)])$note,    from(o.dup));
  name.new    = sapply(tmp1, function(j) { paste(sort(unique(j)), collapse = "; "); });
  db_xref.new = sapply(tmp2, function(j) { paste(sort(unique(j)), collapse = "; "); });
  note.new    = sapply(tmp3, function(j) { paste(sort(unique(j)), collapse = "; "); });
  mcols(result)[as.numeric(names(name.new)),       "name"] = unname(name.new);
  mcols(result)[as.numeric(names(db_xref.new)), "db_xref"] = unname(db_xref.new);
  mcols(result)[as.numeric(names(note.new)),       "note"] = unname(note.new);
  
  #readjust back onto regular chromosomes
  result = GRanges(result$chr.o, ranges(result), strand(result),
                   mcols(result)[,!(colnames(mcols(result)) == "chr.o")]);
  return(result);
}

#argv = c("cdd.aa.gff2.gz", "refseq.to.hg38.gff3.gz", "appris.txt",
#         "cdd.all.igv.gff3", "cdd.isoform_collapse.igv.gff3",
#         "cdd.principal.domains.igv.gff3", "cdd.principal.sites.igv.gff3",
#         "cdd.all.Rdata", "cdd.isoform_collapse.Rdata", 
#         "cdd.principal.domains.Rdata", "cdd.principal.sites.Rdata");
argv = commandArgs(trailingOnly = TRUE);
cdd.aa = load.aa(argv[1]);
cdd.cDNA = aa.to.cDNA(cdd.aa);
map.cDNA_genome = load.cDNA.genome.map(argv[2]);
cdd.all = cDNA.to.genome(cdd.cDNA, map.cDNA_genome);
cdd.c = collapse.domains(cdd.all[grepl("^CDD\\:", cdd.all$db_xref) &
                                 cdd.all$type == "Region"]);
appris = sort(str_replace(read.table(argv[3], header = TRUE)$NP,"\\.\\d+$", ""));
cdd.p = collapse.domains(cdd.all[grepl("^CDD\\:", cdd.all$db_xref) & 
                                 cdd.all$type == "Region" &
                                 cdd.all$NP %in% appris]);
sites.p = collapse.sites(cdd.all[grepl("^CDD\\:", cdd.all$db_xref) &
                         cdd.all$type == "Site" &
                         cdd.all$NP %in% appris]);

cdd.all.browser = format.for.browser(cdd.all, c("type","name","seqnames","NP"));
cdd.c.browser = format.for.browser(cdd.c, c("gene","domain.name","seqnames"));
cdd.p.browser = format.for.browser(cdd.p, c("gene","domain.name","seqnames"));
sites.browser = format.for.browser(sites.p, c("gene", "type", "name", "note", "seqnames"))

export(cdd.all.browser, con = argv[4], format = "GFF3");
export(cdd.c.browser, con = argv[5], format = "GFF3");
export(cdd.p.browser, con = argv[6], format = "GFF3");
export(sites.browser, con = argv[7], format = "GFF3");

save(cdd.all, file = argv[8]);
save(cdd.c, file = argv[9]);
save(cdd.p, file = argv[10]);
save(sites.p, file = argv[11]);

#export(cdd.c, con = argv[8], format = "GFF3");

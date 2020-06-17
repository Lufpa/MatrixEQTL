# this code is borrowed from http://www.bios.unc.edu/research/genomic_software/Multi-Tissue-eQTL/code2.html
# to generate a combined matrix for all SNP_gene pairs tested in cis eQTL mappings runned in MAtrixeQTL
# for each tissue independently


### Read df.txt for the list of tissues
### and degrees of freedom of linear models
df = read.table("df.txt", stringsAsFactors=FALSE);
names(df) = c("tissue", "df");
#show(df)

### List vector for storing Matrix eQTL results
big.list = vector("list", nrow(df));

### Store gene and SNP names from the first tissue
### for matching with other tissues
genes = NULL;
snps = NULL;

### colClasses for faster reading of Matrix eQTL output
cc.file = NA;

### Loop over tissues
for(t1 in 1:nrow(df) ) {
  
  ### Get tissue name
  tissue = df$tissue[t1];
  
  ### Load Matrix eQTL output for the given tissue
  start.time = proc.time()[3];
  tbl = read.table(paste0("multitissue.results.cis.",tissue,"_voom_plate_sva_5PCs"),
                   header = T, stringsAsFactors=FALSE, colClasses=cc.file);
  end.time = proc.time()[3];
  cat(tissue, "loaded in", end.time - start.time, "sec.", nrow(tbl), "gene-SNP pairs.", "\n");

  ### set colClasses for faster loading of other results
  if(any(is.na(cc.file))) {
    cc.file = sapply(tbl, class);
  }

### Set gene and SNP names for matching
  if(is.null(snps))
    snps = unique(tbl$SNP);
  if(is.null(genes))
    genes = unique(tbl$gene);
    
  ### Match gene and SNP names from Matrix eQTL output
  ### to "snps" and "genes"
  gpos = match(tbl$gene, genes, nomatch=0L);
  spos = match(tbl$SNP,  snps,  nomatch=0L);
  
  ### Assign each gene-SNP pair a unique id for later matching with other tissues
  id = gpos + 2*spos*length(genes);
  
  ### Transform t-statistics into correlations
  r = tbl$t.stat / sqrt(df$df[t1] + tbl$t.stat^2);

  ### Record id's and correlations
  big.list[[t1]] = list(id = id, r = r);
  
  ### A bit of clean up to reduce memory requirements
  rm(tbl, gpos, spos, r, id, tissue, start.time, end.time);
  gc();
}

rm(t1, cc.file);

### Find the set of gene-SNP pairs
### present in results for all tissues
keep = rep(TRUE, length(big.list[[1]]$id));
for(t1 in 2:nrow(df)) {
  mch = match(big.list[[1]]$id, big.list[[t1]]$id, nomatch=0L);
  keep[ mch == 0] = FALSE;
  cat(df$tissue[t1], ", overlap size", sum(keep), "\n");
}
final.ids = big.list[[1]]$id[keep];
rm(keep, mch, t1);

### Create and fill in the matrix of z-scores
### Z-scores are calculated from correlations
big.matrix = matrix(NA_real_, nrow=length(final.ids), ncol=nrow(df));
fisher.transform = function(r){0.5*log((1+r)/(1-r))};
for(t1 in 1:nrow(df)) {
  mch = match(final.ids,  big.list[[t1]]$id);
  big.matrix[,t1] = fisher.transform(big.list[[t1]]$r[mch]) * sqrt(df$df[t1]-1);
  cat(t1,"\n");
}
stopifnot( !any(is.na(big.matrix)) )
rm(t1, mch);


### Save the big matrix
save(list="big.matrix", file="z-score.matrix.Rdata", compress=FALSE);

### Save gene names and SNP names for rows of big matrix
writeLines(text=genes[final.ids %% (length(genes)*2)], con="z-score.matrix.genes.txt")
writeLines(text=snps[final.ids %/% (length(genes)*2)], con="z-score.matrix.snps.txt")


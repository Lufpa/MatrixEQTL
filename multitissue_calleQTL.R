#taken from http://www.bios.unc.edu/research/genomic_software/Multi-Tissue-eQTL/code4.html

### Parameters
local.FDR.threshold = 0.05;
output.file.name = "MT-eQTLs_FDR5.txt";


### Load big matrix of z-scores
load(file="z-score.matrix.Rdata");
dim(big.matrix);

### Load gene names and SNP names
### matching the rows of big.matrix
gnames = readLines("z-score.matrix.genes.txt");
snames = readLines("z-score.matrix.snps.txt");

### Load tissue names
df = read.table("df.txt",stringsAsFactors=FALSE);
names(df) = c("tissue","df");
show(df);

### Load parameter estimates and pick the last one
load("paralist.Rdata");
param = tail(paralist,1)[[1]];

### Number of tissues
K = ncol(big.matrix);
m = nrow(big.matrix);

### The function for matrix power
mat.power = function(mat,pow) {
  e = eigen(mat);
  V = e$vectors;
  return( V %*% diag(e$values^pow) %*% t(V) );
}

### Matrix of possible tissue specificity profiles
Pmat = simplify2array(param$Psubs);


###
### Call eQTLs and save in a file
###

fid = file(description=output.file.name, open="wt");
writeLines(con=fid, paste0(
  "SNP\tgene\t",
   paste0("isEQTL.", df$tissue, collapse="\t"),
   "\t",
   paste0("marginalP.", df$tissue, collapse="\t")));


### Do calculations in slices of 10000 gene-SNP pairs
step1 = 1e4L;
cumdump = 0;
for( j in 1:ceiling(nrow(big.matrix)/step1)) {
  fr = step1*(j-1)+1;
  to = min(step1*j,nrow(big.matrix));
  X = big.matrix[fr:to, , drop=FALSE];
  
  ### likelihood for the slice
  prob = matrix(0, nrow(X), length(param$P));
  for( i in 1:length(param$Psubs)) {
    sigma_star = param$Delta + param$Sigma * tcrossprod(param$Psubs[[i]]);
    sigma_hfiv = mat.power(sigma_star,-0.5);
    sigma_dethfiv = (det(sigma_star))^(-0.5);
    w = (1/(2*pi)^(K/2)) * (param$P[i]*sigma_dethfiv);
    prob[,i] = exp(log(w) - colSums(tcrossprod(sigma_hfiv/sqrt(2),X)^2) ) ;
  }
  prob = prob / rowSums(prob);

  ### Select tests with eQTLs significant at local.FDR.threshold level
  keep = (prob[,1] <= local.FDR.threshold);
  if(any(keep)) {
    marginalProb = tcrossprod(prob[keep,,drop=FALSE], 1-Pmat);
    tissueSpecificity = t(Pmat)[
      apply(X=prob[keep,,drop=FALSE], MARGIN=1, FUN=which.max), ];
    
    dump = data.frame(
              snames[(fr:to)[keep]],
              gnames[(fr:to)[keep]],
              tissueSpecificity,
              marginalProb,
              row.names = NULL,
              check.rows = FALSE,
              check.names = FALSE,
              stringsAsFactors = FALSE);
    write.table(dump, file = fid, quote = FALSE,
                sep = "\t", row.names = FALSE, col.names = FALSE);
  }
  cumdump = cumdump + sum(keep);
  cat("Slice",j,"of",ceiling(nrow(big.matrix)/step1)," eQTLs recorded:",cumdump,"\n");
}
close(fid);



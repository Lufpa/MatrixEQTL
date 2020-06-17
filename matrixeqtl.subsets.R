# code for G eQTL mapping in matrixeQTL
# always activate conda R_MatrixEQTL environment before running the R script or any script

library("MatrixEQTL")
library("dplyr")
library("naniar")

args = commandArgs(trailingOnly=T)
cat("arguments: \n")
cat(args, "\n")

covfile  <- args[1]
analysis <- args[2]
vst      <- args[3]
nor      <- args[4]
PCs      <- args[5]
PCs      <- as.numeric(PCs)
file     <- args[6]
covfree  <- args[7]

cat ('covfree analysis :' , covfree, "\n")

### 
### FIRST: subset the genotype, phenotype, and covariate files according to the desired analysis
###

  # the cov files for the different subsets already contained customized svs, so the subsetting
  # is gonna be based on those pre-existing cov files so these analysis matches the diff. gene expression

  # load covariate file depending on the analysis to be run

    #load covariates  
    cvrt = SlicedData$new()
    cvrt$fileDelimiter="\t"
    cvrt$fileOmitCharacters="NA"
    cvrt$fileSkipRows=1
    cvrt$fileSkipColumns=1
    cvrt$fileSliceSize=2000
    cvrt$LoadFile(covfile)

    cat("number of samples to process: ", length(colnames(cvrt)), "\n")

   # c <- as.data.frame(colnames(cvrt))
   # c <- gsub("_", " ", c$`colnames(cvrt)`)

   # cat ("writing indtokeep file\n")
   # write.table(c, paste("indtokeep", file, ".txt", sep=""), col.names = F, row.names = F, quote=F)
    
  # load gene expression data and subset samples (this table is universal)
  # use matrixeQTL loading to make it faster
  
    #load gene expresion matrix
    gene = SlicedData$new()
    gene$fileDelimiter="\t"
    gene$fileOmitCharacters="NA"
    gene$fileSkipRows=1
    gene$fileSkipColumns=1
    gene$fileSliceSize=2000
    gene$LoadFile(vst)

    # subset table
    b<- gene$columnNames %in% cvrt$columnNames
    gene$ColumnSubsample(b)
    
    # quantile normalize while keeping rank # this means I will have to import VST data w/o normalize
    # i have previously generated a VST dataset for each subset, I will have to load that here
    # and normalize it here, or outside and import that.
    
    if ( nor == "yes"){
      cat ("vst counts are alerady normalized\n")
    } else {
      cat ("vst counts are being ranked normalized\n")
      genem <- as.matrix(gene)
      vst.norm.NA <- data.frame(matrix(nrow = nrow(genem), ncol=ncol(genem)))
      
      for( i in 1:nrow(genem)) {
        g <- genem[i,]
        g <- as.data.frame(replace_with_na_all(g, ~.x < 5.76))
        g <- t(apply(g,1, rank,ties.method="average", na.last="keep"))
        g <- g/(sum(!is.na(g))+1)
        gg <- qnorm(g)
        vst.norm.NA[i,] <- gg
      }
      
      colnames(vst.norm.NA) <- colnames(genem)
      rownames(vst.norm.NA) <- rownames(genem)
      gene <- SlicedData$new()
      gene$CreateFromMatrix(as.matrix(vst.norm.NA) )
      
	 cat ("writing vst normalized table\n")
    	 write.table(vst.norm.NA, paste("normVST", file, ".txt", sep=""), col.names = T, row.names=T, quote=F,sep="\t")
    	 rm (vst.norm.NA)
	    }

    
  # load genotype data and subset samples (this geno table is universal)   

    #load genotype matrix 
    snps = SlicedData$new()
    snps$fileDelimiter="\t"
    snps$fileOmitCharacters="NA"  #-1 missing data from vcftools matrix were converted to NA
    snps$fileSkipRows=1
    snps$fileSkipColumns=1
    snps$fileSliceSize=2000
    snps$LoadFile("/Genomics/ayroleslab2/lamaya/bigProject/Genotypes_feb2020/headbody_notsameindv/ctrl_headbody.final.forMatrixeQTL.Jun4.txt") #maybe replace this for an argument to be fed when running the script

    # subset table
    a<- snps$columnNames %in% cvrt$columnNames
    snps$ColumnSubsample(a)
        
	# subset covariate file based on the samples that overlapped with genotypes (misisng indv were removed during
	# QC for the snp file)
	 b <- cvrt$columnNames %in% snps$columnNames
	 cvrt$ColumnSubsample(b)
	
			# subset gene expression matrix to match the samples from geno and cvrt
			c <- gene$columnNames %in% snps$columnNames
		         gene$ColumnSubsample(c)
		
	cat ('after matching genotype with phenotype file, number of individuals in analysis: ', ncol(cvrt), "\n") 		
	c <- as.data.frame(colnames(cvrt))
   	c <- gsub("_", " ", c$`colnames(cvrt)`)
	cat ("writing indtokeep file\n")
	write.table(c, paste("indtokeep", file, ".txt", sep=""), col.names = F, row.names = F, quote=F)



    # filter SNPs based on MAF 5% given the subset of samples to be used
    # maf will have to be calculated again for each subset make sure to print how many snps are used
    
    cat ("caculating maf\n")
    maf.list = vector('list', length(snps))
    for(sl in 1:length(snps)) {
     slice = snps[[sl]];
     maf.list[[sl]] = rowMeans(slice,na.rm=TRUE)/2;
     maf.list[[sl]] = pmin(maf.list[[sl]],1-maf.list[[sl]]);
    }
    maf = unlist(maf.list)
  
    cat('SNPs before MAF filtering: ', nrow(snps),"\n")
    snps$RowReorder(maf>=0.05)
    cat('SNPs after MAF filtering: ', nrow(snps), "\n")
    rm(maf)
    write.table(rownames(snps), paste("snpstokeep", file, ".txt", sep=""), col.names = F, row.names = F, quote=F)

    
     
    #get GRM for this set of SNPs
    
	if (PCs == 0) {
		cat ("no GRM will be calculated because no PCs will be added\n")
		if (covfree == "yes") {
			cvrt$LoadFile()
			cat (" no covariate file will be used\n ")
			}
			
	} else {
		cat ("getting GRM and adding PCs to cov matrix\n")
    		outkinship = file
    		system(paste("~/bin/gcta_1.93.1beta/gcta64 --make-grm-bin --bfile /Genomics/ayroleslab2/lamaya/bigProject/Genotypes_feb2020/headbody_notsameindv/ctrl_headbody.final --autosome --keep ", paste("indtokeep",file,".txt",sep=""), " --extract ", paste("snpstokeep",file,".txt",sep="")," --out ", outkinship, sep=""))
    		system(paste("~/bin/gcta_1.93.1beta/gcta64 --grm ", outkinship, "--pca ", PCs, "--out ", outkinship)) 


		#system(paste("~/bin/ldak5.linux.fast --axes 50 --kinship-details NO --pca ", outkinship, " --grm ", outkinship, sep="" ))
		kinship <- read.table(paste(outkinship,".eigenvec",sep=""),h=F)
		kinshipid <- paste(kinship[,1],kinship[,2],sep="_")
		kinship <- kinship[, c(3:(PCs+2))] #extract the desired number of PCs
		kinship <- t(kinship)
		colnames(kinship) <- kinshipid
		rownames(kinship) <- c(paste(rep("PC", PCs), 1:PCs, sep=""))
    		kinship <- kinship[,order(match(colnames(kinship), cvrt$columnNames))]
                 
		#add first # PCA from kinship to cov file at the beginning of file because the last cov is the factor for GxE
		 if (covfree == "yes") {
                        
			cvrt$CreateFromMatrix(as.matrix(kinship))
                        cat (" only kinship PCs will be used as covariates\n ")
                        } else {

		cvrtm <- as.matrix(cvrt)
		cvrtm <- rbind(kinship, cvrtm)
		cvrt  <- SlicedData$new()
		cvrt$CreateFromMatrix(as.matrix(cvrtm) )
		cat (" PCs + covariates will be used as covariates\n ")
    	}
}

    # order files to have samples in the same order
    snps$ColumnSubsample(order(match(snps$columnNames,cvrt$columnNames)))
    all_equal(snps$columnNames, cvrt$columnNames)
    
    gene$ColumnSubsample(order(match(gene$columnNames, cvrt$columnNames)))
    all_equal(gene$columnNames, cvrt$columnNames)
    
    if (isTRUE(all_equal(snps$columnNames, gene$columnNames))){
      cat("all samples are in the right order\n")
    } else {
      cat("ERROR: something went wrong with the sample order in snps, gene, and/or cov\n")
    }
    
    
  #load snp position
  snpspos <- read.table("/Genomics/ayroleslab2/lamaya/bigProject/Genotypes_feb2020/headbody_notsameindv/ctrl_headbody.final.SNPposition.forMatrixeQTL.Jun4.txt", h=F, stringsAsFactors = FALSE)
  #subset table according to SNPs passing MAF filter
  snpspos <- subset(snpspos, snpspos$V1 %in% rownames(snps))
  
  
  #load gene position
  genepos <- read.table("~/scripts/MatrixEqtl/gene_pos_dmel.txt", h=T, stringsAsFactors = FALSE)



  #run MatrixeQTL in the analysis type (G or GxE)
  
  if (analysis == "G"){
    
    cat("matrixEQTL running in modelLINEAR with diet as covariate\n")
    G_eQTL = Matrix_eQTL_main(
      snps = snps,
      gene = gene,
      cvrt = cvrt,
      output_file_name.cis = paste("results.G.cis10Kb.", file, sep=""),
      output_file_name = paste("results.G.trans.", file, sep=""),
      pvOutputThreshold.cis = 1.5e-2,
      pvOutputThreshold = 3.5e-4,
      useModel = modelLINEAR, 
      errorCovariance = numeric(), 
      cisDist = 1e4,
      snpspos = snpspos, 
      genepos = genepos,
      verbose = TRUE,
      pvalue.hist = "qqplot",
      min.pv.by.genesnp = FALSE,
      noFDRsaveMemory = FALSE)
    
    saveRDS(G_eQTL, file= paste("matrixEQTL_resultsobject_G.", file, ".rds", sep="")) 
    
  
} else {
  
  cat("matrixEQTL running in modelLINEAR_CROSS with genotype-by-diet effect\n")
  G_eQTL = Matrix_eQTL_main(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name.cis = paste("results.GxE.cis10Kb.", file, sep=""),
    output_file_name = paste("results.GxE.trans.", file, sep=""),
    pvOutputThreshold.cis = 4e-5,
    pvOutputThreshold = 6e-8,
    useModel = modelLINEAR_CROSS,
    errorCovariance = numeric(),
    cisDist = 1e4,
    snpspos = snpspos,
    genepos = genepos,
    verbose = TRUE,
    pvalue.hist = "qqplot",
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE)
  
  saveRDS(G_eQTL, paste("matrixEQTL_resultsobject_GxE.", file, ".rds", sep="")) 
  
}


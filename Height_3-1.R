#Load qtl package for QTL mapping and snow package for multicore processing
library(qtl)
library(snow)
require(snow)

#Set working directory with data
setwd("/Users/jkta/Desktop/Projects/NAM CAM/data/")

#Data checking------------------------------------------------------------------

#Reads QTL data
col_sha <- read.cross(format = "csv", file = "col_sha_qtl_final2.csv",
                      genotypes = c("AA", "BB"))
class(col_sha)[1] <- "riself"

#Plot colxsha data to make sure there are no wonky phenotypes
summary(col_sha)
plot(col_sha)

#Estimates recombination frations between chromosomes to see if there are no linked
#markers between different chromosomes
col_sha <- est.rf(col_sha)
plot.rf(col_sha)

#Estimates a genetic map using MLE
newmap <- est.map(col_sha, verbose = TRUE, error.prob = 0.001)
plot.map(col_sha, newmap)
replace.map(col_sha, newmap)
plot.map(col_sha)

#Calculates error LODs with specific markers, to make sure there are no 
#genotyping errors
col_sha <- calc.errorlod(col_sha, error.prob = 0.001)
top.errorlod(col_sha, cutoff = 3)

#Interval mapping---------------------------------------------------------------

#Simulates genotypes between markers, and calculates their genotype 
#probabilities for use in interval mapping
col_sha <- sim.geno(col_sha, n.draws = 64, step = 1, error.prob = 0.001)
col_sha <- calc.genoprob(col_sha, step = 1, error.prob = 0.001)

#QTL model based on single-QTL analysis of Height 3 - Height 1------------------
h1h3_out_imp <- scanone(col_sha, pheno.col = 8, method = "imp", n.cluster = 4)
h1h3_perm_imp <- scanone(col_sha, pheno.col = 8, method = "imp", n.perm = 5000,
                         n.cluster = 4)
h1h3_perm95 <- summary(h1h3_perm_imp)[1]

#QTLs on Chr 1, 2, 4, and 5
plot(h1h3_out_imp, ylab = "LOD Score")
abline(h = h1h3_perm95, lty = 2)
summary(h1h3_out_imp, perms = h1h3_perm_imp, alpha = 0.05)

#QTL model based on putative QTLs; everything significant
h1h3_col_sha_qtl <- makeqtl(col_sha, chr = c(1, 2, 4, 5), 
                            pos = c(30.4, 39, 76.3, 13))

h1h3_col_sha_fq <- fitqtl(col_sha, pheno.col = 8, qtl = h1h3_col_sha_qtl, 
                          method = "imp", formula = y ~ Q1 + Q2 + Q3 + Q4)

summary(h1h3_col_sha_fq)

#No interactions
addint(col_sha, pheno.col = 8, qtl = h1h3_col_sha_qtl, method = "imp",
       formula = y ~ Q1 + Q2 + Q3 + Q4)

#Refining our QTL model
h1h3_col_sha_rq <- refineqtl(col_sha, pheno.col = 8, qtl = h1h3_col_sha_qtl, 
                             method = "imp", formula = y ~ Q1 + Q2 + Q3 + Q4)

h1h3_col_sha_fq1 <- fitqtl(col_sha, pheno.col = 8, qtl = h1h3_col_sha_rq,
                           method = "imp", formula = y ~ Q1 + Q2 + Q3 + Q4)

summary(h1h3_col_sha_fq1)

#Checking for additional QTLs; possibly another on Chr 5
h1h3_col_sha_aq <- addqtl(col_sha, pheno.col = 8, qtl = h1h3_col_sha_rq,
                          method = "imp", formula = y ~ Q1 + Q2 + Q3 + Q4)

summary(h1h3_col_sha_aq)

#Checking the fit of the additional QTL on Chr 5; everything significant
h1h3_col_sha_qtl1 <- makeqtl(col_sha, chr = c(1, 2, 4, 5, 5), 
                             pos = c(29, 41.7, 72.5, 13, 49))

h1h3_col_sha_fq2 <- fitqtl(col_sha, pheno.col = 8, qtl = h1h3_col_sha_qtl1,
                           method = "imp", formula = y ~ Q1 + Q2 + Q3 + 
                           Q4 + Q5)
summary(h1h3_col_sha_fq2)

#No interactions
addint(col_sha, pheno.col = 8, qtl = h1h3_col_sha_qtl1, method = "imp",
       formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5)

h1h3_col_sha_rq1 <- refineqtl(col_sha, pheno.col = 8, qtl = h1h3_col_sha_qtl1,
                              method = "imp", formula = y ~ Q1 + Q2 + Q3 +
                              Q4 + Q5)
h1h3_col_sha_fq3 <- fitqtl(col_sha, pheno.col = 8, qtl = h1h3_col_sha_rq1,
                           method = "imp", formula = y ~ Q1 + Q2 + Q3 +
                           Q4 + Q5)
summary(h1h3_col_sha_fq3)

#No additional QTLs
h1h3_col_sha_aq1 <- addqtl(col_sha, pheno.col = 8, qtl = h1h3_col_sha_qtl1,
                           method = "imp", formula = y ~ Q1 + Q2 + Q3 +
                           Q4 + Q5)
summary(h1h3_col_sha_aq1)

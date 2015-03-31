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

#QTL model based on single-QTL analysis of Height 3 - Height 2------------------

#QTL model based on single-QTL analysis
h2h3_out_imp <- scanone(col_sha, pheno.col = 7, method = "imp", n.cluster = 4)
h2h3_perm_imp <- scanone(col_sha, pheno.col = 7, method = "imp", n.perm = 5000,
                         n.cluster = 4)

h2h3_perm95 <- summary(h2h3_perm_imp)[1]

#Significant QTLs on Chr 1, 2, 4, 5
plot(h2h3_out_imp, ylab = "LOD Score")
abline(h = h2h3_perm95, lty = 2)
summary(h2h3_out_imp, perms = h2h3_perm_imp, alpha = 0.05)

#Checking the fit of our QTL; everything seems significant
h2h3_col_sha_qtl <- makeqtl(col_sha, chr = c(1, 2, 4, 5), 
                            pos = c(30.4, 40.8, 72.5, 26.7))
h2h3_col_sha_fq <- fitqtl(col_sha, pheno.col = 7, qtl = h2h3_col_sha_qtl,
                          method = "imp", formula = y ~ Q1 + Q2 + Q3 + Q4)
summary(h2h3_col_sha_fq)

#Checking for interactions in our 4 QTL model; seems like there is an 
#interaciton between QTLs on Chr 1 and 2
addint(col_sha, pheno.col = 7, qtl = h2h3_col_sha_qtl, method = "imp",
       formula = y ~ Q1 + Q2 + Q3 + Q4)

#Checking the fit with the new interaction; everything is significant
h2h3_col_sha_fq1 <- fitqtl(col_sha, pheno.col = 7, qtl = h2h3_col_sha_qtl,
                           method = "imp", formula = y ~ Q1 * Q2 + Q3 + Q4)
summary(h2h3_col_sha_fq1)

#Checking for additional QTLs; none present
h2h3_col_sha_aq <- addqtl(col_sha, pheno.col = 7, qtl = h2h3_col_sha_qtl,
                          method = "imp", formula = y ~ Q1 * Q2 + Q3 + Q4)
summary(h2h3_col_sha_aq)

#Refining the locations of our QTLs
h2h3_col_sha_rq <- refineqtl(col_sha, pheno.col = 7, qtl = h2h3_col_sha_qtl,
                             method = "imp", formula = y ~ Q1 * Q2 + Q3 + Q4)

#Checking the fit of our refined QTL model
h2h3_col_sha_fq2 <- fitqtl(col_sha, pheno.col = 7, qtl = h2h3_col_sha_rq,
                           method = "imp", formula = y ~ Q1 * Q2 + Q3 + Q4)
summary(h2h3_col_sha_fq2)

h2h3_Q1_Q2 <- find.marker(col_sha, chr = c(1, 2), pos = c(29, 40.8))
par(mfrow = c(1, 2))
plot.pxg(col_sha, marker = h2h3_Q1_Q2, pheno.col = 7)
effectplot(col_sha, mname1 = h2h3_Q1_Q2[1], mname2 = h2h3_Q1_Q2[2], 
           pheno.col = 7, add.legend = FALSE)

#QTL model based on 2D scan of Height 3 - Height 2------------------------------


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

#QTL mapping of Height 3--------------------------------------------------------

h3_out_imp <- scanone(col_sha, pheno.col = 4, method = "imp", n.cluster = 4)
h3_perm_imp <- scanone(col_sha, pheno.col = 4, method = "imp", n.perm = 5000, 
                       n.cluster = 4)
h3_perm95 <- summary(h3_perm_imp)[1]

plot(h3_out_imp, ylab = "LOD Score")
abline(h = h3_perm95, lty = 2)
summary(h3_out_imp, perms = h3_perm_imp, alpha = 0.05)

#QTL model for Height 3
h3_col_sha_qtl <- makeqtl(col_sha, chr = c(1, 2, 4, 5), 
                          pos = c(30.4, 39, 76.3, 12))
h3_col_sha_fq <- fitqtl(col_sha, pheno.col = 4, qtl = h3_col_sha_qtl, 
                        method = "imp", formula = y ~ Q1 + Q2 + Q3 + Q4)
summary(h3_col_sha_fq)

#Might be an interaction between QTL1 and QTL2
addint(col_sha, pheno.col = 5, qtl = h3_col_sha_qtl, method = "imp", 
       formula = y ~ Q1 + Q2 + Q3 + Q4)

#Making another QTL fit to see if the interaction is significant
h3_col_sha_fq1 <- fitqtl(col_sha, pheno.col = 4, qtl = h3_col_sha_qtl,
                         method = "imp", formula = y ~ Q1 * Q2 + Q3 + Q4)

#Interaction doesn't seem that significant
summary(h3_col_sha_fq1)

#Checking for additional QTLs
h3_col_sha_aq <- addqtl(col_sha, pheno.col = 4, qtl = h3_col_sha_qtl, 
                        method = "imp", formula = y ~ Q1 + Q2 + Q3 + Q4)

#Seems like there's another QTL on Chr 5
summary(h3_col_sha_aq)

#Make another QTL model with additional QTL on Chr 5 and checking the fit
h3_col_sha_qtl1 <- makeqtl(col_sha, chr = c(1, 2, 4, 5, 5), 
                           pos = c(30.4, 39, 76.3, 12, 49))

h3_col_sha_fq2 <- fitqtl(col_sha, pheno.col = 4, qtl = h3_col_sha_qtl1, 
                         method = "imp", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5) 

#All QTLs seem significant
summary(h3_col_sha_fq2)

#Doesn't seem like there's any interactions
addint(col_sha, pheno.col = 4, qtl = h3_col_sha_qtl1, method = "imp",
       formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5)

#Refining our QTL and checking the fit
h3_col_sha_rq <- refineqtl(col_sha, pheno.col = 4, qtl = h3_col_sha_qtl1, method = "imp",
                           formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5)

h3_col_sha_fq3 <- fitqtl(col_sha, pheno.col = 4, qtl = h3_col_sha_rq, 
                         method = "imp", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5)

summary(h3_col_sha_fq3)

#QTL analysis based on 2D scan of Height 3--------------------------------------

h3_s2_col_sha <- scantwo(col_sha, pheno.col = 4, method = "imp", verbose = TRUE,
                         n.cluster = 4)
h3_s2_col_sha_perm <- scantwo(col_sha, pheno.col = 4, method = "imp", 
                              n.perm = 5000, n.cluster = 4)

#The 2D scan seems to correspond with our single-QTL model
summary(h3_s2_col_sha_perm)
summary(h3_s2_col_sha, thresholds = c(5.56, 4.14, 3.47, 4.33, 2.4))

#Bayesian intervals for our QTLs------------------------------------------------

h3_bayes1 <- bayesint(h3_col_sha_rq, qtl.index = 1, prob = 0.95, 
                      expandtomarkers = TRUE)
h3_bayes2 <- bayesint(h3_col_sha_rq, qtl.index = 2, prob = 0.95, 
                      expandtomarkers = TRUE)
h3_bayes4 <- bayesint(h3_col_sha_rq, qtl.index = 3, prob = 0.95, 
                      expandtomarkers = TRUE)
h3_bayes51 <- bayesint(h3_col_sha_rq, qtl.index = 4, prob = 0.95, 
                       expandtomarkers = TRUE)
h3_bayes52 <- bayesint(h3_col_sha_rq, qtl.index = 5, prob = 0.95, 
                       expandtomarkers = TRUE)

#MQM analysis Height 3----------------------------------------------------------

#Setting every 5th marker as a cofactor and then analyzing for QTLs
auto_col_sha <- mqmautocofactors(aug_col_sha, 50)
h3_auto_mqm_col_sha <- mqmscan(aug_col_sha, auto_col_sha, pheno.col = 4, 
                               n.cluster = 4)

set_col_sha <- mqmsetcofactors(aug_col_sha, 5)
h3_set_mqm_col_sha <- mqmscan(aug_col_sha, set_col_sha, pheno.col = 4, 
                              n.cluster = 4)
summary(h3_auto_mqm_col_sha)
summary(h3_set_mqm_col_sha)

plot(h3_auto_mqm_col_sha, h3_set_mqm_col_sha, h3_out_imp, 
     col = c("blue", "green", "red"), lty = 1:3)

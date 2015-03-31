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

#QTL analysis on Height 2 - Height 1--------------------------------------------
h1h2_out_imp <- scanone(col_sha, pheno.col = 6, method = "imp")
h1h2_perm_imp <- scanone(col_sha, pheno.col = 6, method = "imp", n.perm = 5000,
                         n.cluster = 4)
h1h2_perm95 <- summary(h1h2_perm_imp)[1]

#QTLs on Chr 1 and 5
plot(h1h2_out_imp, ylab = "LOD Score")
abline(h = h1h2_perm95, lty = 2)
summary(h1h2_out_imp, perms = h1h2_perm_imp, alpha = 0.05)

#2 QTL model
h1h2_col_sha_qtl <- makeqtl(col_sha, chr = c(1, 5), pos = c(39.7, 12))
h1h2_col_sha_fq <- fitqtl(col_sha, pheno.col = 6, qtl = h1h2_col_sha_qtl,
                          method = "imp", formula = y ~ Q1 + Q2)

#All QTLs significant
summary(h1h2_col_sha_fq)

#Checking for interactions; doesn't seem like there are any
addint(col_sha, pheno.col = 6, qtl = h1h2_col_sha_qtl, method = "imp",
       formula = y ~ Q1 + Q2)

#Checking for additional QTLs; seems like there are some of Chr 4 and another on
#5
h1h2_col_sha_aq <- addqtl(col_sha, qtl = h1h2_col_sha_qtl, pheno.col = 6, 
                          method = "imp", formula = y ~ Q1 + Q2)
summary(h1h2_col_sha_aq)

#Checking the fit of our new QTL model
h1h2_col_sha_qtl1 <- makeqtl(col_sha, chr = c (1, 4, 5, 5), 
                             pos = c(39.7, 62, 12, 79.6))

#Everything seems significant
h1h2_col_sha_fq1 <- fitqtl(col_sha, pheno.col = 6, qtl = h1h2_col_sha_qtl1,
                           method = "imp", formula = y ~ Q1 + Q2 + Q3 + Q4)
summary(h1h2_col_sha_fq1)

#Checking for interactions - none
addint(col_sha, pheno.col = 6, qtl = h1h2_col_sha_qtl1, method = "imp",
       formula = y ~ Q1 + Q2 + Q3 + Q4)

#Refining the locaitons of our QTLs
h1h2_col_sha_rq <- refineqtl(col_sha, pheno.col = 6, qtl = h1h2_col_sha_qtl1,
                             method = "imp", formula = y ~ Q1 + Q2 + Q3 + Q4)

#Checking for additional QTLs
h1h2_col_sha_aq1 <- addqtl(col_sha, pheno.col = 6, qtl = h1h2_col_sha_rq,
                           method = "imp", formula = y ~ Q1 + Q2 + Q3 + Q4)
summary(h1h2_col_sha_aq1)

#2D scan based on QTL analysis of Height 2 - Height 1---------------------------
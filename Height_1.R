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

#Height 1 QTL analysis----------------------------------------------------------

#Doesn't seem like there are any significant QTLs for height 1, which is 
#interesting; but that might be because we start measuring flowering time only 
#only when they start bolting (so all the individuals are ~ 1 mm)

h1_out_imp <- scanone(col_sha, pheno.col = 2, method = "imp")
h1_out_perm <- scanone(col_sha, pheno.col = 2, method = "imp", n.perm = 2000, 
                       n.cluster = 4)
summary(h1_out_perm)
h1_perm95 <- summary(h1_out_perm)[1]

plot(h1_out_imp, ylab = "LOD Score", ylim = c(0, 2))
abline(h = h1_perm95, lty = 2)

#QTL analysis based on 2D scan of Height 1--------------------------------------
h1_s2_col_sha <- scantwo(col_sha, pheno.col = 2, method = "imp", 
                         verbose = TRUE, n.cluster = 12)
h1_s2_col_sha_perm <- scantwo(col_sha, pheno.col = 2, method = "imp", 
                              n.perm = 5000, n.cluster = 12)

#There are no significant QTL 
summary(h1_s2_col_sha_perm)
summary(h1_s2_col_sha, thresholds = c(485, 484, 426, 242, 241))

aug_col_sha <- mqmaugment(col_sha, minprob =  0.925, verbose = TRUE)
geno.image(aug_col_sha)

#MQM analysis Height 1----------------------------------------------------------

#Setting automatic cofactors (50 of them)
auto_col_sha <- mqmautocofactors(aug_col_sha, 50)
auto_h1_mqm_col_sha <- mqmscan(aug_col_sha, auto_col_sha, pheno.col = 2, 
                               n.cluster = 4)

#Setting every 5th marker as a cofactor and then analyzing for QTLs
set_col_sha <- mqmsetcofactors(aug_col_sha, 5)
set_h1_mqm_col_sha <- mqmscan(aug_col_sha, set_col_sha, pheno.col = 2, 
                              n.cluster = 4)
summary(auto_h1_mqm_col_sha)
summary(set_h1_mqm_col_sha)

par(mfrow = c(2, 1))
plot(auto_h1_mqm_col_sha, set_h1_mqm_col_sha, col = c("blue", "green"), 
     lty = 1:2)

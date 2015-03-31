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

#Building out model for Height 2------------------------------------------------

#Permutes the data to derive LOD scores; this is done by "freezing" the 
#genotypes and randomly assigning phenotypes. Threshold for background noise.
h2_out_imp <- scanone(col_sha, method = "imp", pheno.col = 3)
h2_perm_imp <- scanone(col_sha, method = "imp", pheno.col = 3, n.perm = 5000, 
                       verbose = TRUE, n.cluster = 4)

#Assign the 5% significance threshold so that we can visually inspect the LOD 
#score cutoff.
summary(h2_perm_imp)
h2_perm95 <- summary(h2_perm_imp)[1]

#Plotting the single-QTL analysis with LOD threshold; evidence of QTLs on 
#Chr 1, 4, and 5 (possibly 2 on 5)
plot(h2_out_imp, ylab = "LOD Score", bandcol = "gray90", main = "Imputation QTL on Height (Week 2)")
abline(h = h2_perm95, lty = 2)
summary(h2_out_imp, perms = h2_perm_imp, alpha = 0.05, pvalues = TRUE)

#Single QTL model---------------------------------------------------------------

#Building a single QTL model by adding highest LOD marker on Chr 5
h2_col_sha_qtl <- makeqtl(col_sha, chr = 5, pos = 12)
summary(h2_col_sha_qtl)

#Checking for other significant QTLs
h2_col_sha_aq <- addqtl(col_sha, pheno.col = 3, qtl = h2_col_sha_qtl, 
                        method = "imp")

#Statistically checks for other QTLs based on our permutation data; additional
#QTLs on Chr 1, 4 and another on 5
summary(h2_col_sha_aq, perms = h2_perm_imp, alpha = 0.05, pvalues = TRUE)
plot(h2_col_sha_aq, ylab = "LOD Score")

#Adding additional QTLs to our model with other putative QTL positions
h2_col_sha_qtl1 <- makeqtl(cross = col_sha, chr = c(1, 4, 5, 5), 
                           pos = c(37.7, 80.7, 12, 79.6))

#Refines the location of our QTLs using MLE
h2_col_sha_rq <- refineqtl(col_sha, qtl = h2_col_sha_qtl1, method = "imp", 
                           pheno.col = 3)

h2_col_sha_fq <- fitqtl(col_sha, qtl = h2_col_sha_qtl1, pheno.col = 3, 
                        method = "imp", formula = y ~ Q1 + Q2 + Q3 + Q4)

summary(h2_col_sha_fq)

#Checks for any interactions between QTLs in our 4 QTL model; possible 
#interaction between chromosome 4 and 5
addint(col_sha, qtl = h2_col_sha_qtl1, method = "imp", 
       pheno.col = 3)

#QTL model based on 2D scan on Height 2-----------------------------------------

#2D scan for Height 2
h2_s2_col_sha <- scantwo(col_sha, pheno.col = 3, method = "imp", 
                         verbose = TRUE, n.cluster = 4)
h2_s2_col_sha_perm <- scantwo(col_sha, pheno.col = 3, method = "imp", 
                              n.perm = 5000, n.cluster = 8)

#Generates thresholds for our 2D scan
summary(h2_s2_col_sha_perm)

#Determines putative QTLs and their locations based on our thresholds
summary(h2_s2_col_sha, perms = h2_s2_col_sha_perm, 
        thresholds = c(5.59, 4.15, 3.54, 4.38, 2.42))
summary(h2_col_sha_rq)

#Confirms our QTL model based on our single-QTL analysis
#QTLs on Chr 1, 4, and 5
plot(h2_s2_col_sha, main = "2D QTL Scan")
plot(h2_col_sha_rq)

#Effect plots-------------------------------------------------------------------

#Finding markers based on the positions of significant QTLs in 2D scan
h2_Q2_Q3 <- find.marker(col_sha, chr = c(4, 5), pos = c(62, 13))

#Phenotype and effect plots between markers; all paired markers seem additive
#except for significant markers on Chr 4 and 5
par(mfrow = c(1, 2))
plot.pxg(col_sha, marker = h2_Q2_Q3, pheno.col = 3)
effectplot(col_sha, mname1 = h2_Q2_Q3[1], mname2 = h2_Q2_Q3[2], pheno.col = 3, 
           add.legend = FALSE)

#Comparing putative QTLs between our single-QTL analysis and 2D scan
summary(h2_out_imp, perms = h2_perm_imp, alpha = 0.05, pvalues = TRUE)
summary(h2_s2_col_sha, perms = h2_s2_col_sha_perm, pvalues = TRUE, 
        thresholds = c(5.59, 4.15, 3.54, 4.38, 2.42))

#QTL model based on our 2D results
h2_s2_col_sha_qtl <- makeqtl(col_sha, chr = c(1, 4, 5, 5), 
                             pos = c(40, 81, 13, 80))
plot(h2_s2_col_sha_qtl)
h2_s2_col_sha_fq <- fitqtl(col_sha, qtl = h2_s2_col_sha_qtl, pheno.col = 3, 
                           formula = y ~ Q1 + Q2 * Q3 + Q4, method = "imp")

#QTLs seem significant except for the interaction between Chr 4 and 5
summary(h2_s2_col_sha_fq)

#Refine positions of our 2D QTL model
h2_s2_col_sha_rq <- refineqtl(col_sha, qtl = h2_s2_col_sha_qtl, pheno.col = 3,
                              method = "imp", formula = y ~ Q1 + Q2 * Q3 + Q4)
summary(h2_s2_col_sha_rq)

#Fit QTL model with our new positions; all QTLs are still significant, but the
#interaction between Chr 4 and 5 (between pos 12.2 and 12.6, respectively) are
#now significant
h2_s2_col_sha_fq1 <- fitqtl(col_sha, qtl = h2_s2_col_sha_rq, pheno.col = 3,
                            method = "imp", formula = y ~ Q1 + Q2 * Q3 + Q4)
summary(h2_s2_col_sha_fq1)

#Checks to see if there are additional QTLs that could be added to our model 
#based on a single-QTL scan
h2_s2_col_sha_aq <- addqtl(col_sha, qtl = h2_s2_col_sha_rq, pheno.col = 3, 
                           method = "imp", 
                           formula = y ~ Q1 + Q2 * Q3 + Q4)

#Additional putative QTLs on Chr 2 and 4; their LOD > 2
plot(h2_s2_col_sha_aq, ylab = "LOD Score")
summary(h2_s2_col_sha_aq)
plot(h2_s2_col_sha_rq)

#Bayesian intervals for our QTLs------------------------------------------------
h2_bayes1 <- bayesint(h2_s2_col_sha_rq, qtl.index = 1, prob = 0.95, 
                      expandtomarkers = TRUE)
h2_bayes2 <- bayesint(h2_s2_col_sha_rq, qtl.index = 2, prob = 0.95, 
                      expandtomarkers = TRUE)
h2_bayes3 <- bayesint(h2_s2_col_sha_rq, qtl.index = 3, prob = 0.95, 
                      expandtomarkers = TRUE)
h2_bayes4 <- bayesint(h2_s2_col_sha_rq, qtl.index = 4, prob = 0.95, 
                      expandtomarkers = TRUE)

#Extended QTL model on Height 2 (LOD scores only ~ 2)---------------------------
h2_s2_col_sha_qtl1 <- makeqtl(col_sha, chr = c(1, 2, 4, 4, 5, 5), 
                              pos = c(37.7, 36, 12.3, 80.7, 12, 77.3))
plot(h2_s2_col_sha_qtl1)

h2_s2_col_sha_fq2 <- fitqtl(col_sha, qtl = h2_s2_col_sha_qtl1, pheno.col = 3, 
                            formula = y ~ Q1 + Q2 + Q3 * Q5 + Q4 + Q6, 
                            method = "imp")

#All QTLs seem significant
summary(h2_s2_col_sha_fq2)

#Refining our new expanded model
h2_s2_col_sha_rq1 <- refineqtl(col_sha, qtl = h2_s2_col_sha_qtl1, 
                               pheno.col = 3, method = "imp", 
                               formula = y ~ Q1 + Q2 + Q3 * Q5 + Q4 + Q6)
summary(h2_s2_col_sha_rq1)

#Refined QTL model LOD increased by 0.32
h2_s2_col_sha_fq3 <- fitqtl(col_sha, qtl = h2_s2_col_sha_rq1, pheno.col = 3, 
                            formula = y ~ Q1 + Q2 + Q3 * Q5 + Q4 + Q6, 
                            method = "imp")

#Comparison of 4 QTL model vs 6 QTL model; 22 vs 30 LOD, and 54% vs 65% variance
#explained
summary(h2_s2_col_sha_fq3)
summary(h2_s2_col_sha_fq1)

#Possible interaction between both QTLs on Chr 5
addint(col_sha, qtl = h2_s2_col_sha_rq1, 
       formula = y ~ Q1 + Q2 + Q3 * Q5 + Q4 + Q6, 
       pheno.col = 3, method = "imp")

#Checking for additional QTLs; no significant QTLs
h2_s2_col_sha_aq1 <- addqtl(col_sha, qtl = h2_s2_col_sha_rq1, pheno.col = 3, 
                            method = "imp", 
                            formula = y ~ Q1 + Q2 + Q3 * Q5 + Q4 + Q6)
plot(h2_s2_col_sha_aq1, ylab = "LOD Score")
summary(h2_s2_col_sha_aq1)

#MQM model on Height 2----------------------------------------------------------

#Augments the data (basically imputes the data), by adding additional invidiuals
aug_col_sha <- mqmaugment(col_sha, minprob =  0.925, verbose = TRUE)
geno.image(aug_col_sha)

#Comparison of single-QTL analysis scan (on original data) and mqmscan on 
#augmented data
one_h2_col_sha <- scanone(col_sha, pheno.col = 3, method = "imp")
mqm_h2_col_sha <- mqmscan(aug_col_sha, pheno.col = 3, n.clusters = 4)
plot(one_mqm_col_sha, scan_mqm_col_sha, col = c("red", "blue"), lty = 1:2)
abline(h = h2_perm95, lty = 2)
summary(scan_mqm_col_sha)

#Finding markers that have significant LOD scores in the mqmscan, and then 
#setting them as cofactors to find additional QTLs
h2_mqm_markers <- find.marker(aug_col_sha, chr = c(1, 4, 5), pos = c(40, 80, 10))
h2_mqm_mts <- find.markerindex(aug_col_sha, name = h2_mqm_markers)
h2_mqm_cofactors <- mqmsetcofactors(aug_col_sha, cofactors = h2_mqm_mts)
h2_mqm_col_sha_co1 <- mqmscan(aug_col_sha, cofactors = h2_mqm_cofactors, 
                              pheno.col = 3, n.cluster = 4)
plot(h2_mqm_col_sha_co1)
summary(h2_mqm_col_sha_co1)

#Doing backwards method of finding QTLs-----------------------------------------

#Setting automatic cofactors (50 of them)
auto_col_sha <- mqmautocofactors(aug_col_sha, 50)
h2_auto_mqm_col_sha <- mqmscan(aug_col_sha, auto_col_sha, pheno.col = 3, 
                               n.cluster = 4)

#Setting every 5th marker as a cofactor and then analyzing for QTLs
set_col_sha <- mqmsetcofactors(aug_col_sha, 5)
h2_set_mqm_col_sha <- mqmscan(aug_col_sha, set_col_sha, pheno.col = 3, 
                              n.cluster = 4)
summary(h2_auto_mqm_col_sha)
summary(h2_set_mqm_col_sha)

plot(h2_auto_mqm_col_sha, h2_set_mqm_col_sha, h2_out_imp, 
     col = c("blue", "green", "red"), lty = 1:3)

#Checking putative QTL locations in mqm backwards models and forward selection 
#QTL model
par(mfrow = c(2, 2))
plot(mqmgetmodel(h2_auto_mqm_col_sha))
plot(mqmgetmodel(h2_set_mqm_col_sha))
plot(h2_s2_col_sha_rq1)
plot(h2_s2_col_sha_rq)

#Setting permutations to find QTL significance thresholds
mqm.perm <- mqmpermutation(aug_col_sha, scanfunction = mqmscan, 
                           cofactors = h2_mqm_cofactors, pheno.col = 3, 
                           batchsize = 25, n.perm = 2000)
mqm.perm.process <- mqmprocesspermutation(mqm.perm)
summary(mqm.perm.process)
mqmplot.permutations(mqm.perm, legend = FALSE)

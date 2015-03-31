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

#Starting interval mapping------------------------------------------------------

#Simulates genotypes between markers, and calculates their genotype 
#probabilities for use in interval mapping
col_sha <- sim.geno(col_sha, n.draws = 64, step = 1, error.prob = 0.001)
col_sha <- calc.genoprob(col_sha, step = 1, error.prob = 0.001)

#Building a QTL model for bolt days---------------------------------------------

#Starting out with an imputation model
bd_out_imp <- scanone(col_sha, method = "imp", pheno.col = 5, n.cluster = 4)
bd_perm_imp <- scanone(col_sha, method = "imp", pheno.col = 5, n.perm = 2000, 
                       verbose = TRUE, n.cluster = 4)
summary(bd_perm_imp)
bd_perm_95 <- summary(bd_perm_imp)[1]

#Seems like there are QTLs on Chr 1, 4, and 5 (possibly 2 on Chr 5)
plot(bd_out_imp, ylab = "LOD score", main = "Imputation Analysis on Bolt Days")
abline(h = bd_perm_95, lty = 2)

#Again this verifies the location of QTLs on Chr 1, 4, and 5
summary(bd_out_imp, perms = bd_perm_imp, alpha = 0.05)

#Building a QTL model using the significant QTL's from out single-QTL analysis
bd_col_sha_qtl <- makeqtl(col_sha, chr = c(1, 4, 5), pos = c(66.02, 6.34, 12))
plot(bd_col_sha_qtl)

#Fitting our QTL model; all putative QTLs are significant
bd_col_sha_fq <- fitqtl(col_sha, pheno.col = 5, qtl = bd_col_sha_qtl, 
                        method = "imp", formula = y ~ Q1 + Q2 + Q3)
summary(bd_col_sha_fq)

bd_col_sha_rq <- refineqtl(col_sha, qtl = bd_col_sha_qtl, pheno.col = 5,
                           method = "imp", formula = y ~ Q1 + Q2 + Q3)

#QTL 1 moved a bit
par(mfrow = c(1, 2))
plot(bd_col_sha_qtl)
plot(bd_col_sha_rq)

bd_col_sha_fq1 <- fitqtl(col_sha, pheno.col = 5, qtl = bd_col_sha_rq, 
                         method = "imp", formula = y ~ Q1 + Q2 + Q3)
summary(bd_col_sha_fq1)

#Checking for interactions between QTLs in our model; seems like QTL2 (Chr 4) 
#interacts with QTL3 (Chr 5)
addint(col_sha, qtl = bd_col_sha_rq, pheno.col = 5, method = "imp",
       formula = y ~ Q1 + Q2 + Q3)

#Adding interactions to our model, all QTLs and their interactions are 
#significant
bd_col_sha_fq2 <- fitqtl(col_sha, pheno.col = 5, qtl = bd_col_sha_rq, 
                         method = "imp", formula = y ~ Q1 + Q2 + Q2 * Q3)
summary(bd_col_sha_fq2)

#Checking for additional QTLs; none really significant
bd_col_sha_aq <- addqtl(col_sha, qtl = bd_col_sha_rq, pheno.col = 5, 
                        method = "imp", formula = y ~ Q1 + Q2 + Q2 * Q3)

#Building a QTL model based on 2D scan of bolt days-----------------------------

#2D QTL scan on bolt days
bd_s2_col_sha <- scantwo(col_sha, pheno.col = 5, method = "imp")
bd_s2_col_sha_perm <- scantwo(col_sha, pheno.col = 5, method = "imp",
                              n.perm = 5000)

summary(bd_s2_col_sha_perm)

#Seems like there are only significant QTLs on Chr 4 and 5
summary(bd_s2_col_sha, perms = bd_s2_col_sha_perm, 
        thresholds = c(6.59, 5.38, 3.84, 4.69, 3.36))
plot(bd_s2_col_sha)

bd_s2_col_sha_qtl <- makeqtl(col_sha, chr = c(4, 5), pos = c(3, 15))
bd_s2_col_sha_fq <- fitqtl(col_sha, pheno.col = 5, qtl = bd_s2_col_sha_qtl, 
                           method = "imp", formula = y ~ Q1 * Q2)

#These are definitely QTLs because they account for ~60% variance
summary(bd_s2_col_sha_fq)

#Refining locations of bolt days QTL
bd_s2_col_sha_rq <- refineqtl(col_sha, pheno.col = 5, qtl = bd_s2_col_sha_qtl, 
                              method = "imp", formula = y ~ Q1 * Q2)
summary(bd_s2_col_sha_rq)

#Checking for additional QTLs; there is another putative QTL on Chr 1
bd_s2_col_sha_aq <- addqtl(col_sha, pheno.col = 5, qtl = bd_s2_col_sha_qtl,
                           method = "imp", formula = y ~ Q1 * Q2)
summary(bd_s2_col_sha_aq)

#Fitting our new 3 QTL model
bd_s2_col_sha_qtl2 <- makeqtl(col_sha, chr = c(1, 4, 5), pos = c(59.16, 4, 14))
bd_s2_col_sha_fq1 <- fitqtl(col_sha, pheno.col = 5, qtl = bd_s2_col_sha_qtl2,
                            method = "imp", formula = y ~ Q1 + Q2 * Q3)
summary(bd_s2_col_sha_fq1)

#Checking for additional QTLs; doesn't seem like there are any
bd_s2_col_sha_aq1 <- addqtl(col_sha, pheno.col = 5, qtl = bd_s2_col_sha_qtl2,
                            method = "imp", formula = y ~ Q1 + Q2 * Q3)
summary(bd_s2_col_sha_aq1)

#Checking for interactions between our 3 QTL model; no interactions
addint(col_sha, pheno.col = 5, qtl = bd_s2_col_sha_qtl2, method = "imp",
       formula = y ~ Q1 + Q2 * Q3)

#Refining our model
bd_s2_col_sha_rq1 <- refineqtl(col_sha, pheno.col = 5, qtl = bd_s2_col_sha_qtl2,
                               method = "imp", formula = y ~ Q1 + Q2 * Q3)

bd_s2_col_sha_fq2 <- fitqtl(col_sha, pheno.col = 5, qtl = bd_s2_col_sha_rq1,
                            method = "imp", formula = y ~ Q1 + Q2 * Q3)

summary(bd_s2_col_sha_fq2)

#Building confidence intervals for our QTLs to find candidate genes-------------
bd_bayes1 <- bayesint(bd_s2_col_sha_rq1, prob = 0.95, qtl.index = 1,
                      expandtomarkers = TRUE)
bd_bayes4 <- bayesint(bd_s2_col_sha_rq1, prob = 0.95, qtl.index = 2,
                      expandtomarkers = TRUE)
bd_bayes5 <- bayesint(bd_s2_col_sha_rq1, prob = 0.95, qtl.index = 3,
                      expandtomarkers = TRUE)

bayesint(bd_out_imp, chr = 5, prob = 0.95, expandtomarkers = TRUE)

bd_bayes1
bd_bayes4
bd_bayes5

#Bolt days effect plots---------------------------------------------------------
bd_Q2_Q3 <- find.marker(col_sha, chr = c(4, 5), pos = c(3, 15))
par(mfrow = c(1, 2))
plot.pxg(col_sha, marker = bd_Q2_Q3, pheno.col = 5)
effectplot(col_sha, mname1 = bd_Q2_Q3[1], mname2 = bd_Q2_Q3[2], pheno.col = 5, 
           add.legend = FALSE)

#MQM model of bolt days QTL-----------------------------------------------------
aug_col_sha <- mqmaugment(col_sha, minprob =  0.925, verbose = TRUE)
geno.image(aug_col_sha)

#Comparison of single-QTL analysis scan (on original data) and mqmscan on 
#augmented data
one_bd_mqm_col_sha <- scanone(col_sha, pheno.col = 5, method = "imp")
scan_bd_mqm_col_sha <- mqmscan(aug_col_sha, pheno.col = 5, n.clusters = 4)
plot(one_bd_mqm_col_sha, scan_bd_mqm_col_sha, col = c("red", "blue"), lty = 1:2)
abline(h = perm95)
summary(scan_bd_mqm_col_sha)

#Finding markers that have significant LOD scores in the mqmscan, and then 
#setting them as cofactors to find additional QTLs
markers_bd_mqm <- find.marker(aug_col_sha, chr = c(1, 4, 5), pos = c(65, 5, 10))
mqm_bd_mts <- find.markerindex(aug_col_sha, name = markers_bd_mqm)
mqm_bd_cofactors <- mqmsetcofactors(aug_col_sha, cofactors = mqm_bd_mts)
co1_bd_mqm_col_sha <- mqmscan(aug_col_sha, cofactors = mqm_bd_cofactors, 
                              pheno.col = 5, n.cluster = 4)
plot(co1_bd_mqm_col_sha)
summary(co1_bd_mqm_col_sha)

#Bolt days backwards method of finding QTLs-------------------------------------

#Setting automatic cofactors (50 of them)
auto_bd_col_sha <- mqmautocofactors(aug_col_sha, 50)
auto_bd_mqm_col_sha <- mqmscan(aug_col_sha, auto_bd_col_sha, pheno.col = 5, 
                               n.cluster = 4)

#Setting every 5th marker as a cofactor and then analyzing for QTLs
set_bd_col_sha <- mqmsetcofactors(aug_col_sha, 5)
set_bd_mqm_col_sha <- mqmscan(aug_col_sha, set_bd_col_sha, pheno.col = 5, 
                              n.cluster = 4)
summary(auto_bd_mqm_col_sha)
summary(set_bd_mqm_col_sha)
summary(co1_bd_mqm_col_sha)

par(mfrow = c(2, 1))
plot(auto_bd_mqm_col_sha, set_bd_mqm_col_sha, co1_bd_mqm_col_sha, 
     col = c("blue", "green", "red"), lty = 1:3)

#Checking putative QTL locations in mqm backwards models and forward selection 
#QTL model
par(mfrow = c(2, 2))
plot(mqmgetmodel(auto_bd_mqm_col_sha))
plot(mqmgetmodel(set_bd_mqm_col_sha))
plot(bd_col_sha_qtl2)
summary(set_bd_mqm_col_sha)

#It seems like choosing the starting qtl map is arbitrary because the methods
#produce similar maps; using setcofactors seem to add extraneous QTLs
bd_mqm_col_sha_qtl <- makeqtl(col_sha, chr = c(1, 2, 3, 4, 5), 
                              pos = c(75, 5, 15, 5, 10))
plot(bd_mqm_col_sha_qtl)
bd_mqm_col_sha_fq <- fitqtl(col_sha, pheno.col = 5, qtl = bd_mqm_col_sha_qtl, 
                            method = "imp", 
                            formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5)
summary(bd_mqm_col_sha_fq)
addint(col_sha, pheno.col = 5, qtl = bd_mqm_col_sha_qtl, 
       method = "imp", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5)
bd_mqm_col_sha_fq1 <- fitqtl(col_sha, pheno.col = 5, qtl = bd_mqm_col_sha_qtl,
                             method = "imp", 
                             formula = y ~ Q1 + Q2 * Q3 + Q4 * Q5)
summary(bd_mqm_col_sha_fq1)
bd_mqm_col_sha_rq <- refineqtl(col_sha, pheno.col = 5, qtl = bd_mqm_col_sha_qtl,
                               method = "imp",
                               formula = y ~ Q1 + Q2 * Q3 + Q4 * Q5)
plot(bd_mqm_col_sha_rq)
bd_mqm_col_sha_fq2 <- fitqtl(col_sha, pheno.col = 5, qtl = bd_mqm_col_sha_rq,
                             method = "imp",
                             formula = y ~ Q1 + Q2 * Q3 + Q4 * Q5)
summary(bd_mqm_col_sha_fq2)

bd_mqm_col_sha_aq <- addqtl(col_sha, pheno.col = 5, qtl = bd_mqm_col_sha_rq, 
                            method = "imp", 
                            formula = y ~ Q1 + Q2 * Q3 + Q4 * Q5)
summary(bd_mqm_col_sha_aq)
plot(bd_mqm_col_sha_aq)

par(mfrow = c(1, 2))
plot(bd_col_sha_qtl2)
summary(bd_col_sha_qtl2)
summary(bd_mqm_col_sha_rq)
plot(bd_mqm_col_sha_rq)
#Load qtl package for QTL mapping and snow package for multicore processing
library(qtl)
library(snow)
require(snow)

#Set working directory with data
setwd("/Users/jkta/Desktop/data/")

#Data checking------------------------------------------------------------------

#Reads QTL data
col_sha <- read.cross(format = "csv", file = "col_sha_qtl_final2.csv",
                      genotypes = c("AA", "BB"))

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
col_sha <- calc.genoprob(col_sha, step = 1, error.prob = 0.01)

#Uses a bunch of different methods to analyze QTLs
out.mrh2 <- scanone(col_sha, method = "mr", pheno.col = 3)
out.emh2 <- scanone(col_sha, method = "em", pheno.col = 3)
out.hkh2 <- scanone(col_sha, method = "hk", pheno.col = 3)
out.ehkh2 <- scanone(col_sha, method = "ehk", pheno.col = 3)

par(mfrow = c(2, 2))
plot(out.mrh2, main = "Marker Regression", ylab = "LOD Score", col = "red")
plot(out.emh2, main = "EM Method", ylab = "LOD Score", col = "blue")
plot(out.hkh2, main = "Haley-Knott Regression", ylab = "LOD Score", 
     col = "green")
plot(out.ehkh2, main = "Extended Haley-Knott Regression", ylab = "LOD Score", 
     col = "black")
plot(out.hkh2 - out.emh2, ylim = c(-0.5, 1), ylab = "LOD[HK] - LOD[EM]")
plot(out.ehkh2 - out.emh2, ylim = c(-0.5, 1), ylab = "LOD[EHK] - LOD[EM]")
plot(out.ehkh2 - out.hkh2, ylim = c(-0.5, 1), ylab = "LOD[EHK] - LOD[HK]")

col_sha <- sim.geno(col_sha, step = 1, n.draws = 64, error.prob = 0.001)
out.imph2 <- scanone(col_sha, method = "imp", pheno.col = 3, n.cluster = 4)

#Permutes the data to derive LOD scores; this is done by "freezing" the 
#genotypes and randomly assigning phenotypes. Threshold for background noise.
perm.imph2 <- scanone(col_sha, method = "imp", pheno.col = 3, n.perm = 5000, 
                      verbose = TRUE, n.cluster = 4)

#Assign the 5% significance threshold so that we can visually inspect the LOD 
#score cutoff.
summary(perm.imph2)
perm95 <- summary(perm.imph2)[1]

#Plotting the single-QTL analysis with LOD threshold; evidence of QTLs on 
#Chr 1, 4, and 5 (possibly 2 on 5)
plot(out.imph2, ylab = "LOD Score")
abline(h = perm95, lty = 2)
summary(out.imph2, perms = perm.imph2, alpha = 0.05, pvalues = TRUE)

#Building QTL model for Height 2------------------------------------------------

#Building a single QTL model by adding highest LOD marker on Chr 5
qtl_col_sha <- makeqtl(col_sha, chr = 5, pos = 12.3, what = "draws")
summary(qtl_col_sha)

#Checking for other significant QTLs
col_sha_qtl_model <- addqtl(col_sha, qtl = qtl_col_sha, method = "imp", 
                            pheno.col = 3)

#Statistically checks for other QTLs based on our permutation data; additional
#QTLs on Chr 1, 4 and another on 5
summary(col_sha_qtl_model, perms = perm.imph2, alpha = 0.05, pvalues = TRUE)
plot(col_sha_qtl_model, ylab = "LOD Score")


#Adding additional QTLs to our model with other putative QTL positions
qtl1_col_sha <- makeqtl(cross = col_sha, chr = c(1, 4, 5, 5), 
                pos = c(37.7, 81, 12.3, 79.6), 
                what = "draws")

#Refines the location of our QTLs using MLE
rqtl_col_sha <- refineqtl(col_sha, qtl = qtl1_col_sha, method = "imp", 
                          pheno.col = 3)
summary(fitqtl(col_sha, qtl = rqtl_col_sha, method = "imp", pheno.col = 3))

#Checks for any interactions between QTLs in our 4 QTL model
addint(col_sha, qtl = rqtl_col_sha, qtl.only = TRUE, method = "imp", 
       pheno.col = 3)

plot(rqtl_col_sha)

#2D scan for QTLs; pairwise comparison for each interval location between 
#chromosomes
scantwo_col_sha <- scantwo(col_sha, pheno.col = 3, method = "imp", 
                           verbose = TRUE, n.cluster = 4)

#Generates thresholds for our 2D scan
scantwo_col_sha_perm <- scantwo(col_sha, pheno.col = 3, method = "imp", 
                                n.perm = 1000, n.cluster = 8)
summary(scantwo_col_sha_perm)

#Determines putative QTLs and their locations based on our thresholds
summary(scantwo_col_sha, perms = scantwo_col_sha_perm, thresholds = c(5.54, 4.1, 3.45, 4.43, 2.48))
summary(rqtl_col_sha)

#Confirms our QTL model based on our single-QTL analysis
#QTLs on Chr 1, 4, and 5
plot(scantwo_col_sha, main = "2D QTL Scan")

#Finding markers based on the positions of significant QTLs in 2D scan
mar14 <- find.marker(col_sha, chr = c(1, 4), pos = c(40, 81))
mar15 <- find.marker(col_sha, chr = c(1, 5), pos = c(38, 13))
mar45 <- find.marker(col_sha, chr = c(4, 5), pos = c(62, 13))
mar55 <- find.marker(col_sha, chr = c(5, 5), pos = c(11, 80))

#Phenotype and effect plots between markers; all paired markers seem additive
#except for significant markers on Chr 4 and 5
par(mfrow = c(1, 2))
plot.pxg(col_sha, marker = mar14, pheno.col = 3)
effectplot(col_sha, mname1 = mar14[1], mname2 = mar14[2], pheno.col = 3, 
           add.legend = FALSE)
plot.pxg(col_sha, marker = mar15, pheno.col = 3)
effectplot(col_sha, mname1 = mar15[1], mname2 = mar15[2], pheno.col = 3,
           add.legend = FALSE)
plot.pxg(col_sha, marker = mar45, pheno.col = 3)
effectplot(col_sha, mname1 = mar45[1], mname2 = mar45[2], pheno.col = 3,
           add.legend = FALSE)
plot.pxg(col_sha, marker = mar55, pheno.col = 3)
effectplot(col_sha, mname1 = mar55[1], mname2 = mar55[2], pheno.col = 3,
           add.legend = FALSE)

#Comparing putative QTLs between our single-QTL analysis and 2D scan
summary(out.imph2, perms= perm.imph2, alpha = 0.05, pvalues = TRUE)
summary(scantwo_col_sha, perms = scantwo_col_sha_perm, thresholds = c(5.54, 4.1, 3.45, 4.43, 2.48))

#QTL model based on our 2D results
col_sha <- sim.geno(col_sha, step = 0.1, n.draws = 128, error.prob = 0.001)
qtl_col_sha_2D <- makeqtl(col_sha, chr = c(1, 4, 5, 5), pos = c(40, 81, 12, 80))
plot(qtl_col_sha_2D)
qtl_col_sha_2D_fq <- fitqtl(col_sha, qtl = qtl_col_sha_2D, pheno.col = 3, 
                            formula = y ~ Q1 + Q2 * Q3 + Q4, method = "imp")

#QTLs seem significant except for the interaction between Chr 4 and 5
summary(qtl_col_sha_2D_fq)

#Refine positions of our 2D QTL model
qtl_col_sha_2D_ref <- refineqtl(col_sha, qtl = qtl_col_sha_2D, pheno.col = 3,
                                method = "imp", formula = y ~ Q1 + Q2 * Q3 + Q4)
summary(qtl_col_sha_2D_ref)

#Fit QTL model with our new positions; all QTLs are still significant, but the
#interaction between Chr 4 and 5 (between pos 12.2 and 12.6, respectively) are
#now significant
qtl_col_sha_2D_fq1 <- fitqtl(col_sha, qtl = qtl_col_sha_2D_ref, pheno.col = 3,
                          method = "imp", formula = y ~ Q1 + Q2 * Q3 + Q4)
summary(qtl_col_sha_2D_fq1)

#Checks to see if there are additional QTLs that could be added to our model 
#based on a single-QTL scan
qtl_col_sha_2D_aq <- addqtl(col_sha, qtl = qtl_col_sha_2D_ref, pheno.col = 3, 
                            method = "imp", 
                            formula = y ~ Q1 + Q2 * Q3 + Q4)

#Potentially additional putative QTLs on Chr 2 and 4; their LOD > 2
max(qtl_col_sha_2D_aq)
plot(qtl_col_sha_2D_aq, ylab = "LOD Score")
summary(qtl_col_sha_2D_aq)

#New QTL model based on these additional QTLs
qtl_col_sha_2D_rv1 <- makeqtl(col_sha, chr = c(1, 2, 4, 4, 5, 5), 
                              pos = c(37.2, 37.2, 12.2, 62, 12.6, 77.3))
plot(qtl_col_sha_2D_rv1)

qtl_col_sha_2D_fq2 <- fitqtl(col_sha, qtl = qtl_col_sha_2D_rv1, pheno.col = 3, 
                             formula = y ~ Q1 + Q2 + Q3 * Q5 + Q4 + Q6, 
                             method = "imp", get.ests = TRUE)

#All QTLs seem significant
summary(qtl_col_sha_2D_fq2)

#Refining our new expanded model
qtl_col_sha_2D_ref2 <- refineqtl(col_sha, qtl = qtl_col_sha_2D_rv1, 
                                 pheno.col = 3, method = "imp", 
                                 formula = y ~ Q1 + Q2 + Q3 * Q5 + Q4 + Q6)
summary(qtl_col_sha_2D_ref2)

#Refined QTL model LOD increased by 0.8
qtl_col_sha_2D_fq3 <- fitqtl(col_sha, qtl = qtl_col_sha_2D_ref2, pheno.col = 3, 
                             formula = y ~ Q1 + Q2 + Q3 * Q5 + Q4 + Q6, 
                             method = "imp", get.ests = TRUE)

#Comparison of 4 QTL model vs 6 QTL model; 22 vs 30 LOD, and 54% vs 65% variance
#explained
summary(qtl_col_sha_2D_fq3)
summary(qtl_col_sha_2D_fq)

#Possible interaction between both QTLs on Chr 5
addint(col_sha, qtl = qtl_col_sha_2D_ref2, 
       formula = y ~ Q1 + Q2 + Q3 * Q5 + Q4 + Q6, 
       pheno.col = 3, method = "imp")

#The included QTL interaction on Chr 5 don't seem significant in our model
qtl_col_sha_2D_fq4 <- fitqtl(col_sha, qtl = qtl_col_sha_2D_ref2, pheno.col = 3, 
                             formula = y ~ Q1 + Q2 + Q3 * Q5 + Q4 + Q5 * Q6, 
                             method = "imp", get.ests = TRUE)
summary(qtl_col_sha_2D_fq4)

#Checking for additional QTLs; putative QTLs on Chr 1 and 5; LOD > 1
qtl_col_sha_2D_aq1 <- addqtl(col_sha, qtl = qtl_col_sha_2D_ref2, pheno.col = 3, 
                            method = "imp", 
                            formula = y ~ Q1 + Q2 + Q3 * Q5 + Q4 + Q6)
max(qtl_col_sha_2D_aq1)
plot(qtl_col_sha_2D_aq1, ylab = "LOD Score")
summary(qtl_col_sha_2D_aq1)
summary(qtl_col_sha_2D_ref2)

#Adding the 2 QTLs to our model
qtl_col_sha_2D_rv2 <- makeqtl(col_sha, chr = c(1, 1, 2, 4, 4, 5, 5, 5), 
                              pos = c(25.5, 37.2, 39.6, 12.2, 
                                      81.5, 12.6, 55.5, 81.1))
plot(qtl_col_sha_2D_rv2)

#Only the additional QTL on 1 is significant, but the additional QTL on Chr 1 
#might be an artifact because it is close to the other QTL on Chr 1 (<12 cM)
qtl_col_sha_2D_fq5 <- fitqtl(col_sha, qtl = qtl_col_sha_2D_rv2, pheno.col = 3, 
                             formula = y ~ Q1 + Q2 + Q3 + Q4 * Q6 + Q5 + Q7 + Q8, 
                             method = "imp", get.ests = TRUE)
summary(qtl_col_sha_2D_fq5)

#Didn't run because time consuming
qtl_col_sha_2D_ap <- addpair(col_sha, qtl = qtl_col_sha_2D_ref2, pheno.col = 3, 
                             method =  "imp",  
                             formula = y ~ Q1 + Q2 + Q3 * Q5 + Q4 + Q6)

#MQM Code-----------------------------------------------------------------------

#Augments the data (basically imputes the data), by adding additional invidiuals
aug_col_sha <- mqmaugment(col_sha, minprob =  0.925, verbose = TRUE)
geno.image(aug_col_sha)

#Comparison of single-QTL analysis scan (on original data) and mqmscan on 
#augmented data
one_mqm_col_sha <- scanone(col_sha, pheno.col = 3, method = "imp")
scan_mqm_col_sha <- mqmscan(aug_col_sha, pheno.col = 3, n.clusters = 4)
plot(one_mqm_col_sha, scan_mqm_col_sha, col = c("red", "blue"), lty = 1:2)
abline(h = perm95)
summary(scan_mqm_col_sha)

#Finding markers that have significant LOD scores in the mqmscan, and then 
#setting them as cofactors to find additional QTLs
markers_mqm <- find.marker(aug_col_sha, chr = c(1, 4, 5), pos = c(40, 80, 10))
mqm_mts <- find.markerindex(aug_col_sha, name = markers_mqm)
mqm_cofactors <- mqmsetcofactors(aug_col_sha, cofactors = mqm_mts)
mqm_col_sha_co1 <- mqmscan(aug_col_sha, cofactors = mqm_cofactors, 
                           multicore = 4, pheno.col = 3)
plot(mqm_col_sha_co1)
summary(mqm_col_sha_co1)

#Doing backwards method of finding QTLs-----------------------------------------

#Setting automatic cofactors (50 of them)
auto_col_sha <- mqmautocofactors(aug_col_sha, 50)
auto_mqm_col_sha <- mqmscan(aug_col_sha, auto_col_sha, pheno.col = 3, 
                            multicore = 4)

#Setting every 5th marker as a cofactor and then analyzing for QTLs
set_col_sha <- mqmsetcofactors(aug_col_sha, 5)
set_mqm_col_sha <- mqmscan(aug_col_sha, set_col_sha, pheno.col = 3, 
                           multicore = 4)
summary(auto_mqm_col_sha)
summary(set_mqm_col_sha)
par(mfrow = c(2, 1))
plot(auto_mqm_col_sha, set_mqm_col_sha, out.imph2, 
     col = c("blue", "green", "red"), lty = 1:3)

#Checking putative QTL locations in mqm backwards models and forward selection 
#QTL model
par(mfrow = c(2, 2))
plot(mqmgetmodel(auto_mqm_col_sha))
plot(mqmgetmodel(set_mqm_col_sha))
plot(qtl_col_sha_2D_ref2)

#Setting permutations to find QTL significance thresholds
mqm.perm <- mqmpermutation(aug_col_sha, scanfunction = mqmscan, 
                           cofactors = set_col_sha, pheno.col = 3, 
                           batchsize = 25, n.perm = 10)
mqm.perm.process <- mqmprocesspermutation(mqm.perm)
summary(mqm.perm.process)
mqmplot.permutations(mqm.perm, legend = FALSE)


#Building a QTL model for bolt days---------------------------------------------

#Starting out with an imputation model
out.impbd <- scanone(col_sha, method = "imp", pheno.col = 5, n.cluster = 4)
perm.impbd <- scanone(col_sha, method = "imp", pheno.col = 5, n.perm = 2000, 
                      verbose = TRUE, n.cluster = 4)
summary(perm.impbd)
permbd <- summary(perm.impbd)[1]

#Seems like there are QTLs on Chr 1, 4, and 5 (possibly 2 on Chr 5)
plot(out.impbd, ylab = "LOD score")
abline(h = permbd, lty = 2)

#Again this verifies the location of QTLs on Chr 1, 4, and 5
summary(out.impbd, perms = perm.impbd, alpha = 0.05)

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
plot(bd_col_sha_qtl)
plot(bd_col_sha_rq)

bd_col_sha_fq1 <- fitqtl(col_sha, pheno.col = 5, qtl = bd_qtl_col_sha_rq, 
                        method = "imp", formula = y ~ Q1 + Q2 + Q3)
summary(bd_col_sha_fq1)

#Checking for interactions between QTLs in our model; seems like QTL1 (Chr 1) 
#and QTL2 (Chr 4) interact, as well as QTL2 (Chr 4) with QTL3 (Chr 5)
addint(col_sha, qtl = bd_qtl_col_sha_rq, pheno.col = 5, method = "imp",
       formula = y ~ Q1 + Q2 + Q3)

#Adding interactions to our model, all QTLs and their interactions are 
#significant
bd_col_sha_fq2 <- fitqtl(col_sha, pheno.col = 5, qtl = bd_qtl_col_sha_rq, 
                         method = "imp", formula = y ~ Q1 * Q2 + Q2 * Q3)
summary(bd_col_sha_fq2)

#Checking for additional QTLs
bd_col_sha_aq <- addqtl(col_sha, qtl = bd_qtl_col_sha_rq, pheno.col = 5, 
                        method = "imp", formula = y ~ Q1 * Q2 + Q2 * Q3)

#Seems like there are additional QTLs on Chr 1, 4, and 5
summary(bd_col_sha_aq)

#Making another model with additional 3 QTLs
bd_col_sha_qtl1 <- makeqtl(col_sha, chr = c(1, 1, 3, 4, 5, 5), 
                           pos = c(39.7, 76.2, 4, 4, 13, 70))
plot(bd_col_sha_qtl1)

#Seems like all the QTLs and the interactions seem significant
bd_col_sha_fq3 <- fitqtl(col_sha, pheno.col = 5, qtl = bd_col_sha_qtl1, 
                         method = "imp", 
                         formula = y ~ Q1 * Q4 + Q2 + Q3 + Q4 * Q5 + Q6)
summary(bd_col_sha_fq3)

#Refining the QTL locations of our new model
bd_col_sha_rq1 <- refineqtl(col_sha, qtl = bd_col_sha_qtl1, pheno.col = 5,
                            method = "imp", 
                            formula = y ~ Q1 * Q4 + Q2 + Q3 + Q4 * Q5 + Q6)
plot(bd_col_sha_rq1)

#Checking for interactions in our 6 QTL model; interactions between Chr 1 and 3,
#Chr 1 and 5, Chr 4 with Chr 1, and Chr 4 with Chr 5
addint(col_sha, qtl = bd_col_sha_rq1, pheno.col = 5, method = "imp",
       formula = y ~ Q1 * Q4 + Q2 + Q3 + Q4 * Q5 + Q6)

#Adding the interactions and checking for the fit of the model; everything is 
#significant again
bd_col_sha_fq4 <- fitqtl(col_sha, pheno.col = 5, qtl = bd_col_sha_rq1, 
                         method = "imp", 
                         formula = y ~ Q1 * Q4 + Q2 + Q3 + Q4 * Q5 + 
                                       Q6 + Q2 * Q4 + Q6 * Q4)
summary(bd_col_sha_fq4)

#Refining our model
bd_col_sha_rq2 <- refineqtl(col_sha, qtl = bd_col_sha_rq1, pheno.col = 5, 
                            method = "imp", 
                            formula = y ~ Q1 * Q4 + Q2 + Q3 + Q4 * Q5 + 
                                          Q6 + Q2 * Q4 + Q6 * Q4)
plot(bd_col_sha_rq2)

#Checking for additional QTLs; additional QTLs on Chr 1 and 5, but the one on 
#Chr 5 pos 20.5 might be an artifact because it is close to QTL 4 
bd_col_sha_aq1 <- addqtl(col_sha, qtl = bd_col_sha_rq2, pheno.col = 5, 
                         method = "imp", 
                         formula = y ~ Q1 * Q4 + Q2 + Q3 + Q4 * Q5 + 
                                       Q6 + Q2 * Q4 + Q6 * Q4)
summary(bd_col_sha_aq1)

summary(bd_col_sha_rq2)

bd_col_sha_qtl2 <- makeqtl(col_sha, chr = c(1, 1, 3, 4, 5, 5),
                           pos = c(39.7, 78, 17, 4, 10.6, 70))
par(mfrow = c(1, 2))
plot(bd_col_sha_rq2)
plot(bd_col_sha_qtl2)
bd_col_sha_fq4 <- fitqtl(col_sha, pheno.col = 5, qtl = bd_col_sha_qtl2, 
                         method = "imp", 
                         formula = y ~ Q1 * Q4 + Q2 + Q3 + Q4 * Q5 + Q6 + 
                                   Q2 * Q4 + Q4 * Q6)
summary(bd_col_sha_fq4)
bd_col_sha_rq3 <- refineqtl(col_sha, pheno.col = 5, qtl = bd_col_sha_qtl2, 
                            method = "imp", 
                            formula = y ~ Q1 * Q4 + Q2 + Q3 + Q4 * Q5 + Q6 + 
                                      Q2 * Q4 + Q4 * Q6)
plot(bd_col_sha_rq3)

bd_col_sha_aq2 <- addqtl(col_sha, qtl = bd_col_sha_qtl2, pheno.col = 5, 
                         method = "imp", 
                         formula = y ~ Q1 * Q4 + Q2 + Q3 + Q4 * Q5 + Q6 + 
                                   Q2 * Q4 + Q4 * Q6)
summary(bd_col_sha_aq2)
addint(col_sha, qtl = bd_col_sha_qtl2, pheno.col = 5, method = "imp",
       formula = y ~ Q1 * Q4 + Q2 + Q3 + Q4 * Q5 + Q6 + Q2 * Q4 + Q4 * Q6)
plot(bd_col_sha_qtl2)

bd_col_sha_fq5 <- fitqtl(col_sha, pheno.col = 5, qtl = bd_col_sha_qtl2, 
                         method = "imp", 
                         formula = y ~ Q1 * Q3 + Q1 * Q4 + Q1 * Q6 + Q2 + Q3 + 
                                   Q4 * Q5 + Q6 + Q2 * Q4 + Q4 * Q6)
summary(bd_col_sha_fq5)

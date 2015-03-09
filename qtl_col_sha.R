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
plot(h2_out_imp, ylab = "LOD Score")
abline(h = h2_perm95, lty = 2)
summary(h2_out_imp, perms = h2_perm_imp, alpha = 0.05, pvalues = TRUE)

#Building QTL model for Height 2------------------------------------------------

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
one_mqm_col_sha <- scanone(col_sha, pheno.col = 3, method = "imp")
scan_mqm_col_sha <- mqmscan(aug_col_sha, pheno.col = 3, n.clusters = 4)
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

#MQM model of bolt days QTL--------------------------------------------------
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

#Bolt days effect plots---------------------------------------------------------
bd_Q2_Q3 <- find.marker(col_sha, chr = c(4, 5), pos = c(3, 15))
par(mfrow = c(1, 2))
plot.pxg(col_sha, marker = bd_Q2_Q3, pheno.col = 5)
effectplot(col_sha, mname1 = bd_Q2_Q3[1], mname2 = bd_Q2_Q3[2], pheno.col = 5, 
           add.legend = FALSE)

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


#QTL model based on 2D scan of Height 3 - Height 2------------------------------


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

#Genetic maps of all QTL models for all phenotypes
par(mfrow = c(2, 3))
plot(h2_s2_col_sha_rq, main = "Height 2 QTL Map")
plot(h3_col_sha_rq, main = "Height 3 QTL Map")
plot(bd_s2_col_sha_rq, main = "Bolt Days QTL Map")
plot(h1h2_col_sha_rq, main = "Height 2 - Height 1 QTL Map")
plot(h2h3_col_sha_rq, main = "Height 3 - Height 2 QTL Map")
plot(h1h3_col_sha_rq1, main = "Height 3 - Height 1 QTL Map")

library(qtl)
setwd("/Users/jkta/Desktop/data/")

col_sha <- read.cross(format = "csv", file = "col_sha_qtl_final2.csv",
                     genotypes = c("AA", "BB"))

summary(col_sha)
plot(col_sha)

col_sha <- est.rf(col_sha)
plot.rf(col_sha)

newmap <- est.map(col_sha, verbose = TRUE, error.prob = 0.001)
plot.map(col_sha, newmap)
replace.map(col_sha, newmap)
plot.map(col_sha)

col_sha <- calc.errorlod(col_sha, error.prob = 0.001)
top.errorlod(col_sha, cutoff = 3)

col_sha <- sim.geno(col_sha, n.draws = 64, step = 1, error.prob = 0.001)
col_sha <- calc.genoprob(col_sha, step = 1, error.prob = 0.01)

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
perm.imph2 <- scanone(col_sha, method = "imp", pheno.col = 3, n.perm = 5000, 
                      verbose = TRUE, n.cluster = 4)
summary(perm.imph2)
perm95 <- summary(perm.imph2)[1]

plot(out.imph2, ylab = "LOD Score")
abline(h = perm95, lty = 2)

summary(out.imph2, perms= perm.imph2, alpha = 0.05, pvalues = TRUE)

qtl_col_sha <- makeqtl(col_sha, chr = 5, pos = 12.3, what = "draws")
col_sha_qtl_model <- addqtl(col_sha, qtl = qtl_col_sha, method = "imp", 
                            pheno.col = 3)
summary(col_sha_qtl_model, perms = perm.imph2, alpha = 0.05, pvalues = TRUE)
plot(col_sha_qtl_model, ylab = "LOD Score")
abline(h = perm95, lty = 2)

qtl1_col_sha <- makeqtl(cross = col_sha, chr = c(1, 2, 4, 5, 5), 
                pos = c(37.7, 37.0, 80.7, 79.6, 12.3), 
                what = "draws")

rqtl_col_sha <- refineqtl(col_sha, qtl = qtl1_col_sha, method = "imp", 
                          pheno.col = 3)

options(width = 64)
rqtl_col_sha

summary(fitqtl(col_sha, qtl = rqtl_col_sha, method = "imp", pheno.col = 3))
addint(col_sha, qtl = rqtl_col_sha, qtl.only = TRUE, method = "imp", 
       pheno.col = 3)

plot(rqtl_col_sha)

scantwo_col_sha <- scantwo(col_sha, pheno.col = 3, method = "imp", 
                           verbose = TRUE, n.cluster = 4)
scantwo_col_sha_perm <- scantwo(col_sha, pheno.col = 3, method = "imp", 
                                n.perm = 1000, n.cluster = 8)

summary(scantwo_col_sha_perm)
summary(scantwo_col_sha, perms = scantwo_col_sha_perm, thresholds = c(5.54, 4.1, 3.45, 4.43, 2.48))
summary(rqtl_col_sha)
#Interactions on chromosome 4 with 5, but also an additive QTL too
plot(scantwo_col_sha, main = "2D QTL Scan")

mar14 <- find.marker(col_sha, chr = c(1, 4), pos = c(40, 81))
mar15 <- find.marker(col_sha, chr = c(1, 5), pos = c(38, 13))
mar45 <- find.marker(col_sha, chr = c(4, 5), pos = c(62, 13))
mar55 <- find.marker(col_sha, chr = c(5, 5), pos = c(11, 80))

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

summary(out.imph2, perms= perm.imph2, alpha = 0.05, pvalues = TRUE)
summary(scantwo_col_sha, perms = scantwo_col_sha_perm, thresholds = c(5.54, 4.1, 3.45, 4.43, 2.48))
summary(rqtl_col_sha)

col_sha <- sim.geno(col_sha, step = 0.1, n.draws = 128, error.prob = 0.001)
qtl_col_sha_2D <- makeqtl(col_sha, chr = c(1, 2, 4, 4, 5, 5), pos = c(40, 41, 11, 80, 13, 80))
plot(qtl_col_sha_2D)

qtl_col_sha_2D_fq <- fitqtl(col_sha, qtl = qtl_col_sha_2D, pheno.col = 3, 
                            formula = y ~ Q1 + Q2 + Q3 * Q4 + Q5, method = "imp")
summary(qtl_col_sha_2D_fq)

qtl_col_sha_2D_fq2 <- fitqtl(col_sha, qtl = qtl_col_sha_2D, pheno.col = 3, 
                            formula = y ~ Q1 + Q2 + Q3 * Q4 + Q5, 
                            method = "imp", 
                            dropone = FALSE, get.ests = TRUE)
summary(qtl_col_sha_2D_fq2)

rqtl_col_sha_2D <- refineqtl(col_sha, qtl = qtl_col_sha_2D, pheno.col = 3,
                             method = "imp", 
                             formula = y ~ Q1 + Q2 + Q3 * Q4 + Q5)

qtl_col_sha_2D_fq3 <- fitqtl(col_sha, qtl = rqtl_col_sha_2D, pheno.col = 3, 
                             method = "imp", 
                             formula = y ~ Q1 + Q2 + Q3 * Q4 + Q5)

qtl_col_sha_2D_fq4 <- fitqtl(col_sha, qtl = qtl_col_sha_2D, pheno.col = 3, 
                             method = "imp", 
                             formula = y ~ Q1 + Q2 + Q3 * Q5 + Q4 + Q6)
summary(qtl_col_sha_2D_fq)
summary(qtl_col_sha_2D_fq3)
summary(qtl_col_sha_2D_fq4)

rqtl2_col_sha_2D <- refineqtl(col_sha, qtl = qtl_col_sha_2D, pheno.col = 3, 
                              method = "imp", 
                              formula = y ~ Q1 + Q2 + Q3 * Q5 + Q4 + Q6)
summary(rqtl2_col_sha_2D)

final_qtl_col_sha_2D <- fitqtl(col_sha, qtl = rqtl2_col_sha_2D, pheno.col = 3, 
                               method = "imp", 
                               formula = y ~ Q1 + Q2 + Q3 * Q5 + Q4 + Q6)
summary(final_qtl_col_sha_2D)

plotLodProfile(rqtl2_col_sha_2D, ylab = "Profile LOD Score")

addint(col_sha, qtl = rqtl2_col_sha_2D, 
       formula = y ~ Q1 + Q2 + Q3 * Q5 + Q4 + Q6, 
       pheno.col = 3, method = "imp")

qtl_col_sha_2D_aq <- addqtl(col_sha, qtl = rqtl2_col_sha_2D, pheno.col = 3, 
                            method = "imp", 
                            formula = y ~ Q1 + Q2 + Q3 * Q5 + Q4 + Q6)
max(qtl_col_sha_2D_aq)
plot(qtl_col_sha_2D_aq, ylab = "LOD Score")

#Didn't run because time consuming
qtl_col_sha_2D_ap <- addpair(col_sha, qtl = rqtl_col_sha_2D, pheno.col = 3,  
                             formula = y ~ Q1 + Q2 * Q3 + Q4)
#MQM Code

#Augments the data (basically imputes the data)
aug_col_sha <- mqmaugment(col_sha, minprob =  0.925, verbose = TRUE)
geno.image(aug_col_sha)

mqm_one_col_sha <- scanone(col_sha, pheno.col = 3, method = "imp")
mqm_col_sha <- mqmscan(aug_col_sha, pheno.col = 3, n.clusters = 4)
plot(mqm_one_col_sha, mqm_col_sha, col = c("red", "blue"), lty = 1:2)

max(mqm_col_sha)
mqm_markers <- mqmextractmarkers(mqm_col_sha)
find.marker(aug_col_sha, chr = 5, pos = 10)
mqm_mts <- find.markerindex(aug_col_sha, "Chr5_2572017")
mqm_cofactors <- mqmsetcofactors(aug_col_sha, cofactors = mqm_mts)
mqm_col_sha_co1 <- mqmscan(aug_col_sha, cofactors = mqm_cofactors, pheno.col = 3)
plot(mqm_col_sha_co1)
summary(mqm_col_sha_co1)

mqm_mts2 <- c(mqm_mts, 
              find.markerindex(aug_col_sha, find.marker(aug_col_sha, 1, 40)))
mqm_cofactors2 <- mqmsetcofactors(aug_col_sha, cofactors = mqm_mts2)
mqm_col_sha_co2 <- mqmscan(aug_col_sha, cofactors = mqm_mts2, pheno.col = 3, 
                           multicore = 4)
plot(mqm_col_sha_co2, ylim = c(-100, 100))
summary(mqm_col_sha_co2)

plot(mqm_col_sha_co1, mqm_col_sha_co2, col = c("blue", "green"), lty = 1:2, 
     ylim = c(-50, 100))

mqm_col_sha_aco <- mqmautocofactors(aug_col_sha, 50)
mqm_col_sha_auto <- mqmscan(aug_col_sha, mqm_col_sha_aco, pheno.col = 3, 
                            multicore = 4)
mqm_col_sha_sco <- mqmsetcofactors(aug_col_sha, 5)
mqm_col_sha_back <- mqmscan(aug_col_sha, mqm_col_sha_sco, pheno.col = 3, 
                            multicore = 4)
par(mfrow = c(2, 1))
plot(mqm_col_sha_auto, mqm_col_sha_back)


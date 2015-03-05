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

plot(out.imph2, out.emh2, col = c("blue", "red"), ylab = "LOD Score", lty = 1:2)
abline(h = perm95, lty = 2)

summary(out.imph2, perms= perm.imph2, alpha = 0.05, pvalues = TRUE)

qtl_col_sha <- makeqtl(col_sha, chr = 5, pos = 12.3, what = "draws")
col_sha_qtl_model <- addqtl(col_sha, qtl = qtl_col_sha, method = "imp", 
                            pheno.col = 3)
summary(col_sha_qtl_model, perms = perm.imph2, alpha = 0.1, pvalues = TRUE)
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
summary(scantwo_col_sha, perms = scantwo_col_sha_perm, thresholds = c(5.54, 4.1, Inf, 4.43, 2.48))
summary(rqtl_col_sha)
plot(scantwo_col_sha, lower = "fv1", main = "2D QTL Scan")

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

#MQM Code
#Augments the data (basically imputes the data)
mqm_col_sha <- mqmaugment(col_sha, minprob =  0.925, verbose = TRUE)
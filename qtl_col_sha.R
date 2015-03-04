library(qtl)
setwd("/Users/jkta/Desktop/data/")

col_sha <- read.cross(format = "csv", file = "col_sha_qtl_final2.csv",
                     genotypes = c("AA", "BB"))

summary(col_sha)
plot(col_sha)

col_sha <- est.rf(col_sha)
plot.rf(col_sha)

newmap <- est.map(col_sha, verbose = TRUE, error.prob = 0.01)
plot.map(col_sha, newmap)
replace.map(col_sha, newmap)
plot.map(col_sha)

col_sha <- calc.errorlod(col_sha, error.prob = 0.01)
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
plot(out.hkh2, main = "Haley-Knott Regression", ylab = "LOD Score", col = "green")
plot(out.ehkh2, main = "Extended Haley-Knott Regression", ylab = "LOD Score", col = "black")
plot(out.hkh2 - out.emh2, ylim = c(-0.5, 1), ylab = "LOD[HK] - LOD[EM]")
plot(out.ehkh2 - out.emh2, ylim = c(-0.5, 1), ylab = "LOD[EHK] - LOD[EM]")
plot(out.ehkh2 - out.hkh2, ylim = c(-0.5, 1), ylab = "LOD[EHK] - LOD[HK]")

col_sha <- sim.geno(col_sha, step = 1, n.draws = 64, error.prob = 0.001)
out.imph2 <- scanone(col_sha, method = "imp", pheno.col = 3)
perm.imph2 <- scanone(col_sha, method = "imp", pheno.col = 3, n.perm = 5000, verbose = TRUE)
perm95 <- summary(perm.imph2)[1]

plot(out.imph2, out.emh2, col = c("blue", "red"), ylab = "LOD Score", lty = 1:2)
abline(h = perm95, lty = 2)

summary(out.imph2, perms= perm.imph2, alpha = 0.05, pvalues = TRUE)

qtl_col_sha <- makeqtl(col_sha, chr = 5, pos = 39.5, what = "draws")
col_sha_qtl_model <- addqtl(col_sha, qtl = qtl_col_sha, method = "imp", pheno.col = 3)
summary(col_sha_qtl_model, perms = perm.imph2, alpha = 0.1, pvalues = TRUE)
plot(col_sha_qtl_model, ylab = "LOD Score")
abline(h = perm95, lty = 2)



#MQM Code
#Augments the data (basically imputes the data)
mqm_col_sha <- mqmaugment(col_sha, minprob =  0.925, verbose = TRUE)
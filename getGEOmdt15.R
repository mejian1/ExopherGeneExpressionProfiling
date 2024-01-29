library(GEOquery)
library(limma)
gse9720 <- getGEO("GSE9720", AnnotGPL = FALSE)[[1]]
gse9720$condition <- sub("-.*$", "", gse9720$title)
pData(gse9720) <- pData(gse9720)[, c("title", "geo_accession", "condition")]
gse9720 <- gse9720[, order(gse9720$condition)]

eset <- exprs(gse9720)

fData(gse9720) <- fData(gse9720)[, c("ORF", "ID")]
exprs(gse9720) <- normalizeBetweenArrays(log2(exprs(gse9720)), method="quantile")


#  exprs(gse14308) <- normalizeBetweenArrays(log2(exprs(gse14308)+1), method="quantile")
#  
ess <- list(GSE9720Norm=gse9720)

preloadedDir <- tempdir()

save(ess, file=file.path(preloadedDir, "GSE9720Norm.rda"))

MA <- normalizeWithinArrays(RG)
MA <- normalizeWithinArrays(RG, method="robustspline")
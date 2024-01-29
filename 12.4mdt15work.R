library(GEOquery)
library(limma)
gse9720 <- getGEO("GSE9720", AnnotGPL = TRUE)[[1]]
gse9720$condition <- sub("-.*$", "", gse9720$title)
pData(gse9720) <- pData(gse9720)[, c("title", "geo_accession", "condition")]
gse9720 <- gse9720[, order(gse9720$condition)]

fData(gse9720) <- fData(gse9720)[, c("ORF", "ID")]
df<-data.frame(gse9720)

nrow(df)
ncol(df)

# now transpose
df <- t(df)
df <- as.data.frame(df)

# that worked; now to make it an RGList
RGLista <- new("RGList", df)

# analyze our RG
RGListb <- RG.MA(RGLista)
RGListc <- backgroundCorrect(RGListb, method="none", offset=50)
MA <- MA.RG(RGListc, bc.method="subtract", offset=50)
RG <- RG.MA(MA)
MA <- normalizeWithinArrays(MA, method ="none")


gs <- factor(sml)
groups <- make.names(c("control","exptal"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)
# gset <- gset[complete.cases(exprs(gset)), ]
v <- vooma(gset, design, plot=T)
v$genes <- fData(gset)
fit  <- lmFit(v)
cts <- c(paste(groups[1],"-",groups[2],sep=""))
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)
tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","ORF"))
write.table(tT, file=stdout(), row.names=F, sep="\t")
tT2 <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "P-adj value distribution")
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05, lfc=0)
vennDiagram(dT, circle.col=palette())
t.good <- which(!is.na(fit2$F)) # filter out bad probes
qqt(fit2$t[t.good], fit2$df.total[t.good], main="Moderated t statistic")
# volcano plot (log P-value vs log fold change)
colnames(fit2) # list contrast names
ct <- 1        # choose contrast of interest
# Please note that the code provided to generate graphs serves as a guidance to
# the users. It does not replicate the exact GEO2R web display due to multitude
# of graphical options.
# The following will produce basic volcano plot using limma function:
volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20,
            highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)))

# MD plot (log fold change vs mean log expression)
# highlight statistically significant (p-adj < 0.05) probes
plotMD(fit2, column=ct, status=dT[,ct], legend=F, pch=20, cex=1)
abline(h=0)

# General expression data analysis
ex <- exprs(gset)

# box-and-whisker plot
ord <- order(gs)  # order samples by group
palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
          "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))
par(mar=c(7,4,2,1))
title <- paste ("GSE9720", "/", annotation(gset), sep ="")
boxplot(ex[,ord], boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=gs[ord])
legend("topleft", groups, fill=palette(), bty="n")

# expression value distribution
par(mar=c(4,4,2,1))
title <- paste ("GSE9720", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, group=gs, main=title, legend ="topright")


testMA <- new("MAList", u)
list(preprocessing)

# RGListd <- normalizeWithinArrays(RGListc, layout, method="printtiploess", weights=1,
                   bc.method="subtract", offset=0) 

# gset <- getGEO("GSE9720", GSEMatrix =TRUE, AnnotGPL=FALSE)
# if (length(gset) > 1) idx <- grep("GPL5859", attr(gset, "names")) else idx <- 1
# gset <- gset[[idx]]
# fvarLabels(gset) <- make.names(fvarLabels(gset))
# gsms <- "01111"
# sml <- strsplit(gsms, split="")[[1]]
# ex <- exprs(gset)
# qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
# LogC <- (qx[5] > 100) ||
#  (qx[6]-qx[1] > 50 && qx[2] > 0)
# if (LogC) { ex[which(ex <= 0)] <- NaN
# exprs(gset) <- log2(ex) }
# A <- subset MA to extract A values and do the same for M then plot!
# A <- c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5)
view(dim(MA))
plotMA(MA)
plotMA(RGListd, array = 1, xlab = "A", ylab = "M", main = colnames(RGListd), status = NULL,
       zero.weights = TRUE)
spottypes <- readSpotTypes()


MA <- normalizeWithinArrays(RGListb)


 normalizeWithinArrays(df)
?view(exprs(gse9720))
MA.RG(exprs(gse9720))
read.maimages(exprs(gse9720))
exprs(gse14308) <- normalizeBetweenArrays(log2(exprs(gse14308)+1), method="quantile")
#  
ess <- list(GSE9720Norm=gse9720)

preloadedDir <- tempdir()

save(ess, file=file.path(preloadedDir, "GSE9720Norm.rda"))
view(ess)
# remove NaNs !!
ess <- na.omit(ess)
view(ess)
x

# try removing NaNs again
ESS <- ess[!is.na(ess)]
view(ESS)

# 3rd try
subset(ESS, !is.nan(ESS[1]))

# 4th try
Filter(function(x) !any(is.nan(x), ESS))

# 5th
purrr::keep(ESS, ~all(!is.na(x)))

# 6th
apply(ESS, 1, function(x) x[is.nan(x)])

for(i in ESS)
  for()
    
    a <- c(0.1, NaN, 0.3, 0.4, 0.5)
view(a)
b <- MA.RG(a)

MA <- normalizeWithinArrays(ESS, na.rm=TRUE)
MA <- normalizeWithinArrays(RGLista)
MA <- MA.RG(ESS)
MA <- normalizeWithinArrays(RG, method="robustspline")

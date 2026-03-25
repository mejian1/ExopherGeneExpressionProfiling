# ============================================================
# Differential Gene Expression Analysis - Microarray Dataset
# GSE9720: "The Mediator Subunit MDT-15 Confers Metabolic
# Adaptation to Ingested Material"
# DOI: 10.1371/journal.pgen.1000021
# Platform: GPL5859 (Affymetrix C. elegans Genome Array)
# ============================================================
# Requirements: BiocManager::install(c("GEOquery","limma","Biobase"))

library(Biobase)
library(GEOquery)
library(limma)

# ----------------------------------------------------------------
# 1. Load series and platform data from GEO
# ----------------------------------------------------------------
gset <- getGEO("GSE9720", GSEMatrix = TRUE, AnnotGPL = FALSE)
if (length(gset) > 1) idx <- grep("GPL5859", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# Sanitize feature variable labels so they are valid R names
fvarLabels(gset) <- make.names(fvarLabels(gset))

# ----------------------------------------------------------------
# 2. Assign sample groups
# GSE9720 has 5 samples: GSM246680 (N2 control), GSM246681 (fat-6
# RNAi), GSM246682 (fat-7 RNAi), GSM246683 (fat-6/fat-7 RNAi),
# GSM246684 (mdt-15 RNAi).
# gsms encodes one character per sample:
#   0 = control (G0), 1 = mdt-15 RNAi (G1), X = exclude
# Update this string to reflect the sample order in your download.
# ----------------------------------------------------------------
gsms <- "0XXX1"          # 1 control (GSM246680), exclude fat-6/7 RNAi (GSM246681-3), 1 mdt-15 RNAi (GSM246684)
sml  <- strsplit(gsms, split = "")[[1]]

# Exclude samples marked as "X"
sel  <- which(sml != "X")
sml  <- sml[sel]
gset <- gset[, sel]

# ----------------------------------------------------------------
# 3. Log2-transform expression values (if not already on log scale)
# ----------------------------------------------------------------
ex  <- exprs(gset)
qx  <- as.numeric(quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))
LogC <- (qx[5] > 100) ||
        (qx[6] - qx[1] > 50 && qx[2] > 0) ||
        (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) {
  ex[ex <= 0] <- NaN
  exprs(gset)  <- log2(ex)
}

# ----------------------------------------------------------------
# 4. Build design matrix and fit linear model with limma
# ----------------------------------------------------------------
sml    <- paste0("G", sml)       # e.g. "G0", "G1"
fl     <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)

fit <- lmFit(gset, design)

# Contrast: mdt-15 RNAi (G1) vs control (G0)
cont.matrix <- makeContrasts(G1 - G0, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, proportion = 0.01)

# ----------------------------------------------------------------
# 5. Extract top differentially expressed genes
# ----------------------------------------------------------------
tT <- topTable(fit2, adjust = "fdr", sort.by = "B", number = 250)

# Select informative columns (adjust column names to what is available)
available_cols <- intersect(c("ID", "adj.P.Val", "P.Value", "t", "B", "logFC", "ORF", "Gene.Symbol"),
                            colnames(tT))
tT <- subset(tT, select = available_cols)
write.table(tT, file = "GSE9720_mdt15_DGE_top250.tsv", row.names = FALSE, sep = "\t")
cat("Top 250 DEGs written to GSE9720_mdt15_DGE_top250.tsv\n")

# Full table (all genes)
tT_all <- topTable(fit2, adjust = "fdr", sort.by = "B", number = Inf)
write.table(tT_all, file = "GSE9720_mdt15_DGE_all_genes.tsv", row.names = FALSE, sep = "\t")
cat("All-gene table written to GSE9720_mdt15_DGE_all_genes.tsv\n")

# ----------------------------------------------------------------
# 6. Significant DEG subsets (adj.P.Val < 0.05, |logFC| >= 1)
# ----------------------------------------------------------------
dT      <- decideTests(fit2, adjust.method = "fdr", p.value = 0.05, lfc = 0)
sig_up  <- tT_all[tT_all$adj.P.Val < 0.05 & tT_all$logFC > 0, ]
sig_dn  <- tT_all[tT_all$adj.P.Val < 0.05 & tT_all$logFC < 0, ]
sig_fc2 <- tT_all[tT_all$adj.P.Val < 0.05 & abs(tT_all$logFC) >= 1, ]

write.table(sig_up,  file = "GSE9720_mdt15_DEGs_up.tsv",        row.names = FALSE, sep = "\t")
write.table(sig_dn,  file = "GSE9720_mdt15_DEGs_down.tsv",      row.names = FALSE, sep = "\t")
write.table(sig_fc2, file = "GSE9720_mdt15_DEGs_P0.05_FC2.tsv", row.names = FALSE, sep = "\t")
cat(sprintf("Up-regulated DEGs: %d\n",   nrow(sig_up)))
cat(sprintf("Down-regulated DEGs: %d\n", nrow(sig_dn)))
cat(sprintf("DEGs |logFC|>=1:    %d\n",  nrow(sig_fc2)))

# ----------------------------------------------------------------
# 7. Diagnostic Plots
# ----------------------------------------------------------------

# -- P-value histogram --
hist(tT_all$adj.P.Val,
     col = "grey", border = "white",
     xlab = "Adjusted P-value (FDR)",
     ylab = "Number of genes",
     main = "GSE9720 – Adjusted P-value Distribution")

# -- Venn diagram of up/down regulated genes --
vennDiagram(dT, circle.col = palette())

# -- QQ plot of moderated t-statistics --
t_good <- which(!is.na(fit2$F))
qqt(fit2$t[t_good], fit2$df.total[t_good],
    main = "Moderated t-statistic QQ Plot (GSE9720)")

# -- Volcano plot --
ct <- 1  # first (and only) contrast
volcanoplot(fit2, coef = ct,
            main  = colnames(fit2)[ct],
            pch   = 20,
            highlight = length(which(dT[, ct] != 0)),
            names = rep("+", nrow(fit2)))

# -- MA (mean-difference) plot --
plotMD(fit2, column = ct, status = dT[, ct],
       legend = FALSE, pch = 20, cex = 1,
       main   = "MA Plot – MDT-15 RNAi vs Control (GSE9720)")
abline(h = 0)

# ----------------------------------------------------------------
# 8. Box-and-whisker plot of expression across samples
# ----------------------------------------------------------------
ex_ord <- exprs(gset)[, order(fl)]
fl_ord <- fl[order(fl)]
palette(c("#1B9E77", "#7570B3"))
par(mar = c(7, 4, 2, 1))
boxplot(ex_ord,
        boxwex = 0.6, notch = TRUE,
        main    = paste0("GSE9720 / ", annotation(gset)),
        outline = FALSE, las = 2,
        col     = fl_ord)
legend("topleft", levels(fl), fill = palette(), bty = "n")

# -- Expression value density plot --
par(mar = c(4, 4, 2, 1))
plotDensities(exprs(gset), group = fl,
              main   = paste0("GSE9720 / ", annotation(gset), " Value Distribution"),
              legend = "topright")

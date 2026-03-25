# ============================================================
# Differential Gene Expression Analysis - Microarray Dataset
# GSE27677: "A Fasting-Responsive Signaling Pathway that
# Extends Life Span in C. elegans"
# DOI: 10.1016/j.celrep.2012.12.018
# Platform: GPL200 (Affymetrix C. elegans Genome Array)
# ============================================================
# Requirements: BiocManager::install(c("GEOquery","limma","Biobase"))

library(Biobase)
library(GEOquery)
library(limma)

# ----------------------------------------------------------------
# 1. Load series and platform data from GEO
# ----------------------------------------------------------------
gset <- getGEO("GSE27677", GSEMatrix = TRUE, AnnotGPL = TRUE)
if (length(gset) > 1) idx <- grep("GPL200", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# Sanitize feature variable labels
fvarLabels(gset) <- make.names(fvarLabels(gset))

# ----------------------------------------------------------------
# 2. Assign sample groups
# gsms string: one character per sample, X = exclude
# Groups:
#   0 = control (fed/untreated)
#   1 = fasting 3 hr
#   2 = fasting 6 hr
#   3 = fed (refeeding control)
#   4 = fasting (long-term)
# ----------------------------------------------------------------
gsms <- "0000001122XXXXXXXXXXXX333444XXXXXXXXXXXXXXXXXX"
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
sml    <- paste0("G", sml)
fl     <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)

fit <- lmFit(gset, design)

# Contrasts:
#   G4-G0: long-term fasting vs control (primary comparison)
#   G1-G0: 3 hr fasting vs control
#   G2-G1: 6 hr fasting vs 3 hr fasting (sequential time-point)
#   G3-G2: refed vs 6 hr fasting      (sequential time-point)
#   G4-G3: long-term fasting vs refed  (sequential time-point)
cont.matrix <- makeContrasts(
  G4 - G0,   # long-term fasting vs control
  G1 - G0,   # 3 hr fasting vs control
  G2 - G1,   # 6 hr vs 3 hr fasting
  G3 - G2,   # refed vs 6 hr fasting
  G4 - G3,   # long-term fasting vs refed
  levels = design
)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, proportion = 0.01)

# ----------------------------------------------------------------
# 5. Extract top differentially expressed genes (primary contrast:
#    long-term fasting vs control, G4-G0)
# ----------------------------------------------------------------
tT <- topTable(fit2, adjust = "fdr", sort.by = "B", number = 250,
               coef = "G4 - G0")

# Select informative columns
available_cols <- intersect(
  c("ID", "adj.P.Val", "P.Value", "F", "t", "B", "logFC",
    "Gene.symbol", "Gene.title", "Gene.Symbol"),
  colnames(tT)
)
tT <- subset(tT, select = available_cols)
write.table(tT, file = "GSE27677_fasting_DGE_top250.tsv", row.names = FALSE, sep = "\t")
cat("Top 250 DEGs written to GSE27677_fasting_DGE_top250.tsv\n")

# Full table (all genes, primary contrast)
tT_all <- topTable(fit2, adjust = "fdr", sort.by = "B", number = Inf,
                   coef = "G4 - G0")
write.table(tT_all, file = "GSE27677_fasting_DGE_all_genes.tsv", row.names = FALSE, sep = "\t")
cat("All-gene table written to GSE27677_fasting_DGE_all_genes.tsv\n")

# ----------------------------------------------------------------
# 6. Significant DEG subsets (adj.P.Val < 0.05, |logFC| >= 1)
# ----------------------------------------------------------------
sig_up  <- tT_all[!is.na(tT_all$adj.P.Val) & tT_all$adj.P.Val < 0.05 & tT_all$logFC > 0, ]
sig_dn  <- tT_all[!is.na(tT_all$adj.P.Val) & tT_all$adj.P.Val < 0.05 & tT_all$logFC < 0, ]
sig_fc2 <- tT_all[!is.na(tT_all$adj.P.Val) & tT_all$adj.P.Val < 0.05 & abs(tT_all$logFC) >= 1, ]

write.table(sig_up,  file = "GSE27677_fasting_DEGs_up.tsv",        row.names = FALSE, sep = "\t")
write.table(sig_dn,  file = "GSE27677_fasting_DEGs_down.tsv",      row.names = FALSE, sep = "\t")
write.table(sig_fc2, file = "GSE27677_fasting_DEGs_P0.05_FC2.tsv", row.names = FALSE, sep = "\t")
cat(sprintf("Up-regulated DEGs   (fasting vs control): %d\n", nrow(sig_up)))
cat(sprintf("Down-regulated DEGs (fasting vs control): %d\n", nrow(sig_dn)))
cat(sprintf("DEGs |logFC|>=1     (fasting vs control): %d\n", nrow(sig_fc2)))

# ----------------------------------------------------------------
# 7. Diagnostic Plots
# ----------------------------------------------------------------
dT <- decideTests(fit2, adjust.method = "fdr", p.value = 0.05, lfc = 0)

# -- P-value histogram (primary contrast) --
hist(tT_all$adj.P.Val,
     col = "grey", border = "white",
     xlab = "Adjusted P-value (FDR)",
     ylab = "Number of genes",
     main = "GSE27677 – Adjusted P-value Distribution (Fasting vs Control)")

# -- Venn diagram across contrasts --
vennDiagram(dT, circle.col = palette())

# -- QQ plot of moderated t-statistics --
t_good <- which(!is.na(fit2$F))
qqt(fit2$t[t_good], fit2$df.total[t_good],
    main = "Moderated t-statistic QQ Plot (GSE27677)")

# -- Volcano plot (primary contrast: G4-G0) --
ct <- 1  # first contrast: G4-G0
volcanoplot(fit2, coef = ct,
            main      = colnames(fit2)[ct],
            pch       = 20,
            highlight = length(which(dT[, ct] != 0)),
            names     = rep("+", nrow(fit2)))

# -- MA (mean-difference) plot --
plotMD(fit2, column = ct, status = dT[, ct],
       legend = FALSE, pch = 20, cex = 1,
       main   = "MA Plot – Fasting vs Control (GSE27677)")
abline(h = 0)

# ----------------------------------------------------------------
# 8. Box-and-whisker and density plots across samples
# ----------------------------------------------------------------
labels_map <- c(G0 = "control", G1 = "fasting 3 hr",
                G2 = "fasting 6 hr", G3 = "refeeding control",
                G4 = "fasting (long-term)")

ex_ord <- exprs(gset)[, order(fl)]
fl_ord <- fl[order(fl)]
palette(c("#dfeaf4", "#f4dfdf", "#f2cb98", "#b7d9b7", "#d7a3d0"))
par(mar = c(7, 4, 2, 1))
boxplot(ex_ord,
        boxwex = 0.6, notch = TRUE,
        main    = paste0("GSE27677 / ", annotation(gset), " selected samples"),
        outline = FALSE, las = 2,
        col     = fl_ord)
legend("topleft", labels_map[levels(fl)], fill = palette(), bty = "n")

par(mar = c(4, 4, 2, 1))
plotDensities(exprs(gset), group = fl,
              main   = paste0("GSE27677 / ", annotation(gset), " Value Distribution"),
              legend = "topright")

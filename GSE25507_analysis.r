# Version info: R 4.2.2, Biobase 2.58.0, GEOquery 2.66.0, limma 3.54.0, DESeq2 1.38.3
################################################################
#libraries needed
library(GEOquery)
library(limma)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)


# QUESTION 1 - downloading dataset
# load series and platform data from GEO
gset <- getGEO("GSE25507", GSEMatrix = TRUE, AnnotGPL = TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
# dataset only has GPL570 platform only, the if condition is for universally acceptable code
gset <- gset[[idx]]

# QUESTION 2 - preprocessing and EDA

# make proper column names to match toptable
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- paste0(
  "00000000000000000000000000000000000000000000000000",
  "00000000000000111111111111111111111111111111111111",
  "1111111111111111111111111111111111111111111111"
)
sml <- strsplit(gsms, split = "")[[1]]

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("control", "autism"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~ group + 0, gset)
colnames(design) <- levels(gs)

# General expression data analysis
ex <- normalizeBetweenArrays(exprs(gset)) # Normalise data
ex[which(ex <= 0)] <- NaN

summary(ex)
# box-and-whisker plot
ord <- order(gs) # order samples by group
palette(c(
  "#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
  "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"
))
par(mar = c(7, 4, 2, 1))
title <- paste("GSE25507", "/", annotation(gset), sep = "")
boxplot(ex[, ord], boxwex = 0.6, notch = T, main = title, outline = FALSE, las = 2, col = gs[ord])
legend("topleft", groups, fill = palette(), bty = "n")

# expression value distribution
par(mar = c(4, 4, 2, 1))
title <- paste("GSE25507", "/", annotation(gset), " value distribution", sep = "")
plotDensities(ex, group = gs, main = title, legend = "topright")

# generate pdata
cat("pdata with log transformation")
pdata <- pData(gset)
attributes(pdata)

# generate fdata
cat("fdata with log transformation")
fdata <- fData(gset)
attributes(fdata)

# QUESTION 3 - change after log transformation
exprs(gset) <- log2(ex) # log2 transform
ex <- normalizeBetweenArrays(exprs(gset)) # Normalise data
summary(ex)

# box-and-whisker plot with log transformation
ord <- order(gs) # order samples by group
palette(c(
  "#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
  "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"
))
par(mar = c(7, 4, 2, 1))
title <- paste("GSE25507", "/", annotation(gset), "with log transformation", sep = "")
boxplot(ex[, ord], boxwex = 0.6, notch = T, main = title, outline = FALSE, las = 2, col = gs[ord])
legend("topleft", groups, fill = palette(), bty = "n")

# expression value distribution with log transformation
par(mar = c(4, 4, 2, 1))
title <- paste("GSE25507", "/", annotation(gset), " value distribution with log transformation", sep = "")
plotDensities(ex, group = gs, main = title, legend = "topright")

# generate pdata with log transformation
cat("pdata with log transformation")
pdata <- pData(gset)
attributes(pdata)

# generate fdata with log transformation
cat("fdata with log transformation")
fdata <- fData(gset)
attributes(fdata)

# QUESTION 4 - Perform sdifferential expression analysis using simple t-test with using the limma library
ctrl <- which(pData(gset)$group == "control") # control group
atsm <- which(pData(gset)$group == "autism") # treatment group

# perform t-test
ttest <- apply(ex, 1, function(x) t.test(x[ctrl], x[atsm]))

# calculate log fold change
logfc <- apply(ex, 1, function(x) mean(x[atsm]) - mean(x[ctrl]))

# calculate p-values
pvals <- sapply(ttest, function(x) x$p.value)

# correct p-values using Holm's method
adj_pvals <- p.adjust(pvals, method = "holm")

# create a data frame of results
results <- data.frame(logFC = logfc, PValue = pvals, Adj_PValue = adj_pvals)

# plot volcano plot
ggplot(results, aes(x = logFC, y = -log10(PValue))) +
  geom_point(aes(color = Adj_PValue < 0.15), size = 1) +
  scale_color_manual(values = c("black", "red")) +
  theme_classic() +
  ggtitle("Volcano Plot") +
  xlab("Log2 fold change") +
  ylab("-Log10 p-value")

# QUESTION-5 perform differential expression analysis using the limma packag.
fit <- lmFit(gset, design) # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- c(paste(groups[1], "-", groups[2], sep = ""))
cont.matrix <- makeContrasts(contrasts = cts, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust = "holm", sort.by = "B", number = 250)

tT <- subset(tT, select = c("ID", "adj.P.Val", "P.Value", "t", "B", "logFC", "Gene.symbol", "Gene.title"))
write.table(tT, file = stdout(), row.names = F, sep = "\t")
tT2 <- topTable(fit2, adjust = "holm", sort.by = "B", number = Inf)

# summarize test results as "up", "down" or "not expressed"
dT <- decideTests(fit2, adjust.method = "holm", p.value = 0.15)

# volcano plot (log P-value vs log fold change)
colnames(fit2) # list contrast names
ct <- 1 # choose contrast of interest
volcanoplot(fit2,
  coef = ct, main = colnames(fit2)[ct], pch = 20,
  highlight = length(which(dT[, ct] != 0)), names = rep("+", nrow(fit2))
)

# Build histogram of P-values for all genes to visualise control test results
hist(tT2$adj.P.Val,
  col = "grey", border = "white", xlab = "P-adj",
  ylab = "Number of genes", main = "P-adj value distribution"
)

# Venn diagram of results
vennDiagram(dT, circle.col = palette())

# create Q-Q plot for t-statistic
t.good <- which(!is.na(fit2$F)) # filter out bad probes
qqt(fit2$t[t.good], fit2$df.total[t.good], main = "Moderated t statistic")

# MD plot (log fold change vs mean log expression)
# highlight statistically significant (p-adj < 0.15) probes
plotMD(fit2, column = ct, status = dT[, ct], legend = F, pch = 20, cex = 1)
abline(h = 0)

#QUESTION 7
# Perform enrichment analysis using the set of genes obtained using GSEA
# Set up input data (here, using a vector of example genes)
gene_list <- tT$Gene.symbol

# gene set enrichment analysis using the Gene Ontology Biological Processes database
enrich_result <- enrichGO(
  gene = gene_list,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pvalueCutoff = 0.15,
  qvalueCutoff = 0.2
)

# View the results
enrich_result@result

#QUESTION 8
# Subset the top 50 pathways
top_pathways <- head(enrich_result@result, 50)
print("Top 50 pathways:\n")
top_pathways$Description#QUESTION 9

# Plot the top 50 enriched pathways using ggplot2
dot_plot <- ggplot(top_pathways, aes(x = -log10(pvalue), y = Description))
dot_plot + geom_point(size = 3) + xlab("-log10(p-value)") + ylab("Pathway") # dot plot called and initialised

bar_plot <- ggplot(top_pathways, aes(y = reorder(Description, -log10(pvalue)), x = -log10(pvalue))) +
  geom_bar(stat = "identity", fill = "#8F00FF") +
  ylab("Pathway") +
  xlab("-log10(p-value)") +
  ggtitle("Enrichment Scores") +
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5))
bar_plot # bar plot called

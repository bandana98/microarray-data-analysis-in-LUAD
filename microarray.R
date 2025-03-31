# GSE32863

# install.packages(c("BiocManager", "tidyverse"))
# BiocManager::install(c("GEOquery", "limma", "WGCNA", "clusterProfiler", "org.Hs.eg.db"))

# Load libraries
library(GEOquery)
library(limma)
library(WGCNA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(pheatmap)

setwd("D:/Band_research/microarray/LUAD")

# Download from GEO database 
gse_id <- "GSE32863"  # Replace with actual GEO dataset ID
options('download.file.method.GEOquery' = 'libcurl')  # Use HTTP instead of FTP
gse <- getGEO(gse_id, GSEMatrix = TRUE, AnnotGPL = TRUE)
# Extract expression matrix
exprSet <- exprs(gse[[1]])
dim(exprSet)  
head(exprSet)  

########-------------------------  Annotate Probe IDs with Gene Symbols --------------------
# Extract annotation table
anno <- fData(gse[[1]])

# Map probe IDs to gene symbols
exprSet <- as.data.frame(exprSet)

exprSet$GeneSymbol <- anno$Symbol

# Remove probes without gene symbols
exprSet <- exprSet[!is.na(exprSet$GeneSymbol) & exprSet$GeneSymbol != "", ]

# Aggregate duplicate gene symbols (if any) by mean expression
exprSet <- exprSet %>%
  group_by(GeneSymbol) %>%
  summarise(across(everything(), mean, na.rm = TRUE))  

rownames(exprSet) <- exprSet$GeneSymbol

exprSet$GeneSymbol <- NULL

# Normalize the Data 
boxplot(exprSet[, 1:10], main = "Raw Expression Data", las = 2)

# Log2 transformation
exprSet <- log2(exprSet + 1)

###--------------- DEGs analysis --------------------------

# Extract sample information
pdata <- pData(gse[[1]])

# sample groups 
group <- factor(pdata$`characteristics_ch1.9`)  
group <- factor(pdata$`characteristics_ch1.9`, levels = c("tissue: Normal lung", "tissue: Lung tumor"))
table(group)
design <- model.matrix(~ 0 + group)
colnames(design) <- c("Normal", "Tumor")
contrast_matrix <- makeContrasts(Tumor_vs_Normal = Tumor - Normal, levels = design)
fit <- lmFit(exprSet, design)
fit <- contrasts.fit(fit, contrast_matrix)  # Apply contrast
fit <- eBayes(fit)
deg_results <- topTable(fit, coef = "Tumor_vs_Normal", adjust = "BH", number = Inf)

deg_results <- deg_results[deg_results$adj.P.Val < 0.05 & abs(deg_results$logFC) > 1, ]
head(deg_results) 
write.csv(deg_results, "DEGs_GSE32863.csv")

#----------------------------------- heat map     ----------------------------------------
top_degs <- deg_results[order(deg_results$adj.P.Val), ][1:50, ]
top_expr <- exprSet[rownames(top_degs), ]

# Define annotation for samples (Tumor vs Normal)
annotation_col <- data.frame(Group = group)  # Ensure 'group' is your sample classification vector
rownames(annotation_col) <- colnames(top_expr)

# Generate heatmap
pheatmap(top_expr,
         scale = "row",  # Normalize by row (gene expression)
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "ward.D2",
         annotation_col = annotation_col,
         fontsize = 12,  # Increase font size
         fontsize_row = 10,  # Adjust row font
         fontsize_col = 10,  # Adjust column font
         cellwidth = 12,  # Adjust cell width
         cellheight = 10,  # Adjust cell height
         color = colorRampPalette(c("blue", "white", "red"))(50),  # Improved color scale
         main = "Heatmap of Top 50 DEGs")

# Save high-resolution image (600 dpi)
ggsave("heatmap_high_res.png", width = 8, height = 6, dpi = 600)

#----------------------------------- Volcano plot ----------------------------------------
# Ensure logFC is numeric
deg_results$logFC <- as.numeric(deg_results$logFC)
deg_results$adj.P.Val <- as.numeric(deg_results$adj.P.Val)

# significance thresholds
deg_results$Significance <- "Not Significant"
deg_results$Significance[deg_results$logFC >= 1 & deg_results$adj.P.Val < 0.05] <- "Upregulated"
deg_results$Significance[deg_results$logFC <= -1 & deg_results$adj.P.Val < 0.05] <- "Downregulated"

volcano_plot <- ggplot(deg_results, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "green", "Not Significant" = "gray")) +
  labs(title = "Volcano Plot of DEGs", x = "Log2 Fold Change", y = "-Log10 Adjusted P-Value") +
  theme_minimal()

ggsave("volcano_plot.png", plot = volcano_plot, dpi = 600)


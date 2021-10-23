## Clone project from GitHub here
## Go to RStudio Cloud 
## Select under blue New Project menu: New Project from Git Repository
## Paste in URL for BIOL119 repos: https://github.com/tgirke/BIOL119.git 

## Install required packages
if(!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GEOquery", "limma", "clusterProfiler", "org.Hs.eg.db", "enrichplot", "pheatmap", "RColorBrewer"))

## Import data from GEO
library(GEOquery) 
gds <- getGEO("GDS2778") # Download GDS2778 data from GEO
eset <- GDS2eSet(gds)
pData(eset)[ ,1:2] # Inspect experimental design 
exprs(eset)[1:4, ] # Inspect expression data that are RMA normalized and on log2 scale

## DEG analysis with Limma
library(limma) # Loads limma library.
targets <- pData(eset)
design <- model.matrix(~ -1+factor(as.character(targets$agent)))
colnames(design) <- c("Treatment", "Control") 
fit <- lmFit(eset, design) # Fits a linear model for each gene based on the given series of arrays.
contrast.matrix <- makeContrasts(Control-Treatment, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix) # Computes estimated coefficients and standard errors for a given set of contrasts.
fit2 <- eBayes(fit2)
deg_df <- topTable(fit2, coef=1, adjust="fdr", sort.by="B", number=Inf)
deg_df <- deg_df[!is.na(deg_df$adj.P.Val),]

## Over-represention analysis (ORA) with clusterProfiler
cutoff <- 0.05 # Cutoff to use for filtering on adjusted p-value (FDR)
ids <- deg_df[deg_df$adj.P.Val <= cutoff, "Gene.ID"]
ids <- ids[!grepl("/", ids)] # Removes ambiguous mappings
ids <- ids[grepl("", ids)]
library(clusterProfiler); library(org.Hs.eg.db); library(enrichplot)
ego_mf <- enrichGO(gene=ids, OrgDb=org.Hs.eg.db, ont="MF", pAdjustMethod="BH", pvalueCutoff=0.1, readable=TRUE)
ego_bp <- enrichGO(gene=ids, OrgDb=org.Hs.eg.db, ont="BP", pAdjustMethod="BH", pvalueCutoff=0.1, readable=TRUE)

## Visualize BP and MF ontology results
dim(ego_mf)
head(ego_mf)
barplot(ego_mf, showCategory=10)
goplot(ego_mf)
dim(ego_bp)
head(ego_bp)
barplot(ego_bp, showCategory=10)

## Clustering
library(pheatmap); library("RColorBrewer")
affy_ids <- row.names(deg_df[deg_df$adj.P.Val <= cutoff, ])
deg_ma <- exprs(eset)[affy_ids, ] 
pheatmap(deg_ma, color=brewer.pal(9, "Blues"))





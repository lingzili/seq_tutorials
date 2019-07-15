library("tidyverse")
library("limma")
library("DESeq2")
library("edgeR")
library("pheatmap")
library("RColorBrewer")

# Prepare and load data into DESeq2 ---------------------------------------

# Read in count data
data <- read.delim("~/seq_tutorials/data/GSE76268_rMatrix.txt", row.names = 1, comment.char = "#", stringsAsFactors = FALSE)
View(data)

# Modify the header line to remove the extensions followed by sample ID
colnames(data) <- gsub(".bam", "", colnames(data))
colnames(data)

# Create a data frame of experimental design
sample <- c(
  "SRR3048055", "SRR3048059", "SRR3048060", "SRR3048070", "SRR3048074",
  "SRR3048075", "SRR3048076", "SRR3048083", "SRS1219077", "SRS1219078",
  "SRS1219082", "SRS1219084", "SRS1219085", "SRS1219086", "SRS1219089"
)

condition <- c(
  "alpha", "alpha", "alpha", "beta", "beta",
  "beta", "beta", "beta", "beta", "beta",
  "beta", "alpha", "alpha", "alpha", "alpha"
)

meta <- data.frame(sample, condition, row.names = 1, stringsAsFactors = FALSE)

meta

# Check that sample names match in both files
all(colnames(data) %in% rownames(meta))
all(colnames(data) == rownames(meta))

# Load data into DESeq2
padjCutoff <- 0.1 # If this is not 0.1, we have to tweak DESeq's filtering.

dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~condition)

colData(dds)$condition <- factor(colData(dds)$condition, levels = c("beta", "alpha"))

View(counts(dds))

## Total number of raw counts per sample
colSums(counts(dds))

# Normalize samples
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized = TRUE)

## Total number of normalized counts per sample
colSums(counts(dds, normalized = T))

write.table(normalized_counts, file = "data/normalized_counts.txt", sep = "\t", quote = F, col.names = NA)

## Unsupervised clustering
### PCA
### Transform counts for data visualization
rld <- rlog(dds, blind = TRUE)

### Plot PCA of PC1 and PC2
plotPCA(rld, intgroup = c("condition"))

# Input is a matrix of log transformed values
### Extract the rlog matrix from the object
## assay() is function from the "SummarizedExperiment" package that was loaded when you loaded DESeq2
rld_mat <- assay(rld)
pca <- prcomp(t(rld_mat))

# Create data frame with metadata and PC3 and PC4 values for input to ggplot
df <- cbind(meta, pca$x)
ggplot(df) + geom_point(aes(x = PC3, y = PC4, color = condition))

## Correlation Heatmap

### Compute pairwise correlation values
rld_cor <- cor(rld_mat) ## cor() is a base R function

head(rld_cor) ## check the output of cor(), make note of the rownames and colnames

### Plot heatmap
pheatmap(rld_cor)

# Change colors
heat.colors <- brewer.pal(6, "Blues")
pheatmap(rld_cor,
  color = heat.colors, border_color = NA, fontsize = 10,
  fontsize_row = 10, height = 20
)

# Estimate biological variance and visualize
dds <- estimateDispersions(dds)

## Plot dispersion estimates
plotDispEsts(dds)

## Define contrasts, extract results table, and shrink the log2 fold changes

contrast_oe <- c("sample", "condition", "beta")

res_tableOE_unshrunken <- results(dds, contrast = contrast_oe, alpha = 0.05)

res_tableOE <- lfcShrink(dds, contrast = contrast_oe, res = res_tableOE_unshrunken)


# Differential expression analysis can be done using several available tools, each with differing methodologies.
## Differential gene expression using DESeq2
## We considered genes with an FDR<.1 as significantly differentially expressed.

## Correlation Heatmap
## Spearman correlation heatmap between biological replicates

## Heatmap with pre-defined genes With 50 most variable genes

## Volcano plot

## GOTerm analysis with GOstats

## KEGG pathway with GAGE and clusterProfiler

## Gene discovery rate (detected genes vs. reads)
## See <Gene Discovery Rate in RNA-seq data>(http://www.nxn.se/valent/gene-discovery-rate-in-rna-seq-data)

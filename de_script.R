## Setup
### Bioconductor and CRAN libraries used
library(tidyverse)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)
library(DEGreport)

data <- read.table("data/Mov10_full_counts.txt", header = T, row.names = 1)

meta <- read.table("meta/Mov10_full_meta.txt", header = T, row.names = 1)

class(data)
class(meta)

view(data)
view(meta)

all(colnames(data) %in% rownames(meta))
all(colnames(data) == rownames(meta)) 

dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ sampletype)

dds <- estimateSizeFactors(dds)

normalized_counts <- counts(dds, normalized = TRUE)

write.table(normalized_counts, file="data/normalized_counts.txt", sep="\t", quote=F, col.names=NA)

rld <- rlog(dds, blind = TRUE)
rld_mat <- assay(rld)

plotPCA(rld,intgroup="sampletype")

pca <- prcomp(t(rld_mat))


variance <- pca$sdev^2 / sum(pca$sdev^2)
PC = paste0("PC", 1:length(variance))
variance


dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ sampletype)

dds <- DESeq(dds)

resultsNames(dds)

plotDispEsts(dds)

#For MOV10_OE
contrast_oe <- c("sampletype","MOV10_overexpression","control")
res_tableOE_unshrunken <- results(dds, contrast = contrast_oe, alpha = 0.05)
res_table_OE <- lfcShrink(dds, contrast = contrast_oe,res = res_tableOE_unshrunken)

# Define contrast
contrast_oe <- c("sampletype", "MOV10_overexpression", "control")

resultsNames(dds)

#[1] "Intercept"                                  
#[2] "sampletype_MOV10_knockdown_vs_control"     
#[3] "sampletype_MOV10_overexpression_vs_control"

# Extract results
res_tableOE_unshrunken <- results(dds, contrast = contrast_oe, alpha = 0.05)

# Shrink log fold changes (using apeglm method)
res_table_OE <- lfcShrink(dds, coef = 3, type = "apeglm")

# DGE workshop code -not working
# res_tableOE <- lfcShrink(dds, contrast=contrast_oe, res=res_tableOE_unshrunken)

class(res_table_OE)

mcols(res_table_OE,use.names = T)

res_table_OE %>% data.frame() %>% view()

#For Mov10_Knockdown

dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ sampletype)

dds <- DESeq(dds)

resultsNames(dds)

contrast_kd <- c("sampletype", "MOV10_knockdown", "control")

res_tableKD <- results(dds,contrast = contrast_kd, alpha = 0.05)

res_tableKD <- lfcShrink(dds, coef = 2, type="apeglm")

## Summarize results
summary(res_table_OE)

### Set thresholds (How I will set thresholds for any data, which criteria do I need to focus on?)
padj.cutoff <- 0.05
lfc.cutoff <- 0.58

res_table_OE_tb <- res_table_OE %>% 
  data.frame() %>% 
  rownames_to_column(var="gene") %>% 
  as_tibble()


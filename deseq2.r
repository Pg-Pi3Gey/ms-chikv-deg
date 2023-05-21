if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("apeglm")
BiocManager::install("vsn")
BiocManager::install('EnhancedVolcano')
BiocManager::install("org.Hs.eg.db")
BiocManager::install("VennDetail")
BiocManager::install("clusterProfiler")
BiocManager::install("enrichplot")

install.packages("ggplot2")
install.packages("ggnewscale")
install.packages("tidyverse")

library(org.Hs.eg.db)
library(DESeq2)
library(apeglm)
library(vsn)
library(VennDetail)
library(EnhancedVolcano)
library(ggplot2)
library(ggnewscale)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(tidyverse)

## Simple function for plotting a Volcano plot, returns a ggplot object
volcano.plot <- function(res, datasetName) {
  return(EnhancedVolcano(res, x = 'log2FoldChange', y = 'padj',
                         lab=rownames(res),
                         title = datasetName,
                         subtitle = bquote(italic('FDR <= 0.05 and absolute FC >= 1.5')),
                         # Change text and icon sizes
                         labSize = 3, pointSize = 1.5, axisLabSize=10, titleLabSize=12,
                         subtitleLabSize=8, captionLabSize=10,
                         # Disable legend
                         legendPosition = "none",
                         # Set cutoffs
                         pCutoff = 0.05, FCcutoff = 1.5))
}

setwd("~/DESEQ2_R")

### Step 1: import the featurecounts results into R & prepare countdata ###
countdata <- read.table("featurecounts_output.txt", header=TRUE, row.names=1)
# Remove first five columns (chr, start, end, strand, length)
countdata <- countdata[ ,6:ncol(countdata)]
# Remove rudementary filename and rename accoring to Sample name
colnames(countdata) <- gsub("......results.hisat2.", "", colnames(countdata))
colnames(countdata) <- gsub(".hisat.sorted.bam", "", colnames(countdata))
# Convert to matrix
countdata <- as.matrix(countdata)

# read in sample info
indexed <- read.csv('female_sample.csv', header = TRUE)
sampleinfo <- as.matrix(indexed[, -1])
row.names(sampleinfo) <- indexed[, 1]
# check if sample info is matching with count matrix
all(rownames(sampleinfo) == colnames(countdata))

###  Step 2: Analysis with DESeq2 ###
## DESEQ2 With Days of infection as design ##
dds_doi <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = sampleinfo,
                              design = ~ doi)
# pre-filtering: removing rows with low gene counts
# keeping rows that have at least 10 reads total
keep <- rowSums(counts(dds_doi)) >= 10
dds_doi <- dds_doi[keep,]
# set the factor level
dds_doi$doi <- relevel(dds_doi$doi, ref = "Control")
# Run the DESeq pipeline
dds_doi <- DESeq(dds_doi)

## contrast classes ##
resultsNames(dds_doi)

# Contrast class: D0or1 #
resD0or1 <- results(dds_doi, name="doi_D0or1_vs_Control", alpha=0.05)
resD0or1 <- lfcShrink(dds_doi, coef="doi_D0or1_vs_Control", type="apeglm", res=resD0or1)
resD0or1 <- resD0or1[order(resD0or1$pvalue),]
summary(resD0or1)

#visualize D0or1 Results
volcano.plot(res = resD0or1, datasetName = "Day 0or1 v/s Control")

#annotate Gene Symbol and entrezID & export D0or1 Results
resD0or1$symbol <- mapIds(org.Hs.eg.db, keys = rownames(resD0or1), 
                       keytype = "ENSEMBL", column = "SYMBOL")
resD0or1$entrez <- mapIds(org.Hs.eg.db, keys = rownames(resD0or1), 
                       keytype = "ENSEMBL", column = "ENTREZID")
resD0or1_export <- subset(resD0or1, padj < 0.05)
resD0or1_export <- subset(resD0or1_export, abs(log2FoldChange) >= 1.5)
write.csv(resD0or1_export, 
          file="~/DESEQ2_R/RESULTS/DESEQ2/DEG_D0or1.csv")
summary(resD0or1_export)

# Contrast class: D2 #
resD2 <- results(dds_doi, name="doi_D2_vs_Control", alpha=0.05)
resD2 <- lfcShrink(dds_doi, coef="doi_D2_vs_Control", type="apeglm", res=resD2)
resD2 <- resD2[order(resD2$pvalue),]
summary(resD2)

#visualize D2 Results
volcano.plot(res = resD2, datasetName = "Day 2 v/s Control")

#annotate Gene Symbol and entrezID & export D2 Results
resD2$symbol <- mapIds(org.Hs.eg.db, keys = rownames(resD2), 
                          keytype = "ENSEMBL", column = "SYMBOL")
resD2$entrez <- mapIds(org.Hs.eg.db, keys = rownames(resD2), 
                          keytype = "ENSEMBL", column = "ENTREZID")
resD2_export <- subset(resD2, padj < 0.05)
resD2_export <- subset(resD2_export, abs(log2FoldChange) >= 1.5)
write.csv(resD2_export, 
          file="~/DESEQ2_R/RESULTS/DESEQ2/DEG_D2.csv")
summary(resD2_export)

# Contrast class: D3or4 #
resD3or4 <- results(dds_doi, name="doi_D3or4_vs_Control", alpha=0.05)
resD3or4 <- lfcShrink(dds_doi, coef="doi_D3or4_vs_Control", type="apeglm", res=resD3or4)
resD3or4 <- resD3or4[order(resD3or4$pvalue),]
summary(resD3or4)

#visualize D3or4  Results
volcano.plot(res = resD3or4, datasetName = "Day 3or4 v/s Control")

#annotate Gene Symbol and entrezID & export D3or4 Results
resD3or4$symbol <- mapIds(org.Hs.eg.db, keys = rownames(resD3or4), 
                       keytype = "ENSEMBL", column = "SYMBOL")
resD3or4$entrez <- mapIds(org.Hs.eg.db, keys = rownames(resD3or4), 
                       keytype = "ENSEMBL", column = "ENTREZID")
resD3or4_export <- subset(resD3or4, padj < 0.05)
resD3or4_export <- subset(resD3or4_export, abs(log2FoldChange) >= 1.5)
write.csv(resD3or4_export, 
          file="~/DESEQ2_R/RESULTS/DESEQ2/DEG_D3or4.csv")
summary(resD3or4_export)

# Contrast class: D20 #
resD20 <- results(dds_doi, name="doi_D20_vs_Control", alpha=0.05)
resD20 <- lfcShrink(dds_doi, coef="doi_D20_vs_Control", type="apeglm", res=resD20)
resD20 <- resD20[order(resD20$pvalue),]
summary(resD20)

#visualize D20  Results
volcano.plot(res = resD20, datasetName = "Day 20 v/s Control")

#annotate Gene Symbol and entrezID & export D20 Results
resD20$symbol <- mapIds(org.Hs.eg.db, keys = rownames(resD20), 
                          keytype = "ENSEMBL", column = "SYMBOL")
resD20$entrez <- mapIds(org.Hs.eg.db, keys = rownames(resD20), 
                          keytype = "ENSEMBL", column = "ENTREZID")
resD20_export <- subset(resD20, padj < 0.05)
resD20_export <- subset(resD20_export, abs(log2FoldChange) >= 1.5)
write.csv(resD20_export, 
          file="~/DESEQ2_R/RESULTS/DESEQ2/DEG_D20.csv")
summary(resD20_export)

## venn diagram
# Data preparation
resD0or1.degs <- subset(resD0or1, abs(log2FoldChange) >= 1.5)
resD0or1.degs <- row.names(subset(resD0or1.degs, padj < 0.05))
resD2.degs <- subset(resD2, abs(log2FoldChange) >= 1.5)
resD2.degs <- row.names(subset(resD2.degs, padj < 0.05))
resD3or4.degs <- subset(resD3or4, abs(log2FoldChange) >= 1.5)
resD3or4.degs <- row.names(subset(resD3or4.degs, padj < 0.05))
resD20.degs <- subset(resD20, abs(log2FoldChange) >= 1.5)
resD20.degs <- row.names(subset(resD20.degs, padj < 0.05))

ven <- venndetail(list(D0or1 = resD0or1.degs, D20 = resD20.degs, 
                       D2 = resD2.degs, D3or4 = resD3or4.degs)) 
plot(ven)
detail(ven)
#extract genes of interest from venn Diagram & save into csv file
venn_shared_genes <- getSet(ven, subset = "Shared")
venn_shared_genes$symbol <- mapIds(org.Hs.eg.db, keys = venn_shared_genes$Detail, 
                        keytype = "ENSEMBL", column = "SYMBOL")
venn_shared_genes$entrez <- mapIds(org.Hs.eg.db, keys = venn_shared_genes$Detail, 
                        keytype = "ENSEMBL", column = "ENTREZID")
write.csv(venn_shared_genes, file = "~/DESEQ2_R/RESULTS/venn_shared_genes.csv")

ven_D0or1_D2_D3or4 <- getSet(ven, subset = "D0or1_D2_D3or4")
ven_D0or1_D2_D3or4$symbol <- mapIds(org.Hs.eg.db, keys = ven_D0or1_D2_D3or4$Detail, 
                                   keytype = "ENSEMBL", column = "SYMBOL")
ven_D0or1_D2_D3or4$entrez <- mapIds(org.Hs.eg.db, keys = ven_D0or1_D2_D3or4$Detail, 
                                   keytype = "ENSEMBL", column = "ENTREZID")
write.csv(ven_D0or1_D2_D3or4, file = "~/DESEQ2_R/RESULTS/ven_D0or1_D2_D3or4.csv")

ven_D20 <- getSet(ven, subset = "D20")
ven_D20$symbol <- mapIds(org.Hs.eg.db, keys = ven_D20$Detail, 
                                    keytype = "ENSEMBL", column = "SYMBOL")
ven_D20$entrez <- mapIds(org.Hs.eg.db, keys = ven_D20$Detail, 
                                    keytype = "ENSEMBL", column = "ENTREZID")
write.csv(ven_D20, file = "~/DESEQ2_R/RESULTS/ven_D20.csv")

## DESEQ2 With Condition as design ##
dds_c <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = sampleinfo,
                              design = ~ condition)
# pre-filtering: removing rows with low gene counts
# keeping rows that have at least 10 reads total
keep <- rowSums(counts(dds_c)) >= 10
dds_c <- dds_c[keep,]
# set the factor level
dds_c$condition <- relevel(dds_c$condition, ref = "CTRL")
# Run the DESeq pipeline
dds_c <- DESeq(dds_c)
res_c <-results(dds_c, name="condition_INF_vs_CTRL", alpha=0.05)
res_c <- lfcShrink(dds_c, coef="condition_INF_vs_CTRL", type="apeglm", res=res_c)
res_c <- res_c[order(res_c$pvalue),]
summary(res_c)

#visualize  Results
volcano.plot(res = res_c, datasetName = "Infected v/s Control")
# Count data transformations with Regularized log transformation
vsd_c <- vst(dds_c, blind=TRUE)
plotPCA(vsd_c, intgroup="doi")

#annotate Gene Symbol and entrezID & export Results
res_c$symbol <- mapIds(org.Hs.eg.db, keys = rownames(res_c), 
                        keytype = "ENSEMBL", column = "SYMBOL")
res_c$entrez <- mapIds(org.Hs.eg.db, keys = rownames(res_c), 
                        keytype = "ENSEMBL", column = "ENTREZID")
res_export <- subset(res_c, padj < 0.05)
res_export <- subset(res_export, abs(log2FoldChange) >= 1.5)
write.csv(res_export, 
          file="~/DESEQ2_R/RESULTS/DESEQ2/DEG_all.csv")
summary(res_export)

### Step 3: Enrichment Analysis ###
## From DOI as design ##
#D0or1
# Extract significant results
sigOE_d1 <- subset(resD0or1, padj < 0.05)
sigOE_d1 <- subset(sigOE_d1, abs(log2FoldChange) >= 1.5)
sigOE_d1_genes <- as.character(row.names(sigOE_d1))
sigOE_d1_genes = sort(sigOE_d1_genes, decreasing = TRUE)
# GO over-representation analysis
ego_d1 <- enrichGO(gene = sigOE_d1_genes, 
                universe = as.character(row.names(resD0or1)),
                keyType = "ENSEMBL",
                OrgDb = org.Hs.eg.db, 
                ont = "MF", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05,
                readable = TRUE)
## Output results from GO analysis to a table
cluster_summary_d1 <- data.frame(ego_d1)
write.csv(cluster_summary_d1, "~/DESEQ2_R/RESULTS/clusterprofiler/clusterProfiler_D1.csv")
# barplot 
barplot(ego_d1, 
        drop = TRUE, 
        showCategory = 3, 
        title = "Top GO Molecular Functions For D0or1",
        font.size = 10)
# Cnetplot
#To color genes by log2 fold changes, we need to extract the log2 fold changes from our results table creating a named vector
OE_foldchanges_d1 <- sigOE_d1$log2FoldChange
names(OE_foldchanges_d1) <- as.character(row.names(sigOE_d1))

cnetplot(ego_d1, 
         categorySize="pvalue", 
         showCategory = 3, 
         color.params = list(foldChange = OE_foldchanges_d1), 
         vertex.label.font=6)

#D2
# Extract significant results
sigOE_d2 <- subset(resD2, padj < 0.05)
sigOE_d2 <- subset(sigOE_d2, abs(log2FoldChange) >= 1.5)
sigOE_d2_genes <- as.character(row.names(sigOE_d2))
sigOE_d2_genes = sort(sigOE_d2_genes, decreasing = TRUE)
# GO over-representation analysis
ego_d2 <- enrichGO(gene = sigOE_d2_genes, 
                universe = as.character(row.names(resD2)),
                keyType = "ENSEMBL",
                OrgDb = org.Hs.eg.db, 
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05,
                readable = TRUE)
## Output results from GO analysis to a table
cluster_summary_d2 <- data.frame(ego_d2)
write.csv(cluster_summary_d2, "~/DESEQ2_R/RESULTS/clusterprofiler/clusterProfiler_D2.csv")
# barplot 
barplot(ego_d2, 
        drop = TRUE, 
        showCategory = 3, 
        title = "Top 3 GO Biological Process For D2",
        font.size = 10)
# Cnetplot
#To color genes by log2 fold changes, we need to extract the log2 fold changes from our results table creating a named vector
OE_foldchanges_d2 <- sigOE_d2$log2FoldChange
names(OE_foldchanges_d2) <- as.character(row.names(sigOE_d2))

cnetplot(ego_d2, 
         categorySize="pvalue", 
         showCategory = 3, 
         color.params = list(foldChange = OE_foldchanges_d2), 
         vertex.label.font=6)

#D3
# Extract significant results
sigOE_d3or4 <- subset(resD3or4, padj < 0.05)
sigOE_d3or4 <- subset(sigOE_d3or4, abs(log2FoldChange) >= 1.5)
sigOE_d3or4_genes <- as.character(row.names(sigOE_d3or4))
sigOE_d3or4_genes = sort(sigOE_d3or4_genes, decreasing = TRUE)
# GO over-representation analysis
ego_d3or4 <- enrichGO(gene = sigOE_d3or4_genes, 
                   universe = as.character(row.names(resD3or4)),
                   keyType = "ENSEMBL",
                   OrgDb = org.Hs.eg.db, 
                   ont = "MF", 
                   pAdjustMethod = "BH", 
                   qvalueCutoff = 0.05,
                   readable = TRUE)
## Output results from GO analysis to a table
cluster_summary_d3or4 <- data.frame(ego_d3or4)
write.csv(cluster_summary_d3or4, "~/DESEQ2_R/RESULTS/clusterprofiler/clusterProfiler_d3or4.csv")
# barplot 
barplot(ego_d3or4, 
        drop = TRUE, 
        showCategory = 3, 
        title = "Top 3 GO Molecular Function For D3or4",
        font.size = 10)
# Cnetplot
#To color genes by log2 fold changes, we need to extract the log2 fold changes from our results table creating a named vector
OE_foldchanges_d3or4 <- sigOE_d3or4$log2FoldChange
names(OE_foldchanges_d3or4) <- as.character(row.names(sigOE_d3or4))

cnetplot(ego_d3or4, 
         categorySize="pvalue", 
         showCategory = 3, 
         color.params = list(foldChange = OE_foldchanges_d3or4), 
         vertex.label.font=6)


#D20
# Extract significant results
sigOE_d20 <- subset(resD20, padj < 0.05)
sigOE_d20_genes <- as.character(row.names(sigOE_d20))
sigOE_d20_genes = sort(sigOE_d20_genes, decreasing = TRUE)
# GO over-representation analysis
ego_d20 <- enrichGO(gene = sigOE_d20_genes, 
                      universe = as.character(row.names(resD20)),
                      keyType = "ENSEMBL",
                      OrgDb = org.Hs.eg.db, 
                      ont = "BP", 
                      pAdjustMethod = "BH", 
                      qvalueCutoff = 0.05,
                      readable = TRUE)
## Output results from GO analysis to a table
cluster_summary_d20 <- data.frame(ego_d20)
write.csv(cluster_summary_d20, "~/DESEQ2_R/RESULTS/clusterprofiler/clusterProfiler_d20.csv")
# barplot 
barplot(ego_d20, 
        drop = TRUE, 
        showCategory = 3, 
        title = "Top 3 GO Biological Process For D20",
        font.size = 10)
# Cnetplot
#To color genes by log2 fold changes, we need to extract the log2 fold changes from our results table creating a named vector
OE_foldchanges_d20 <- sigOE_d20$log2FoldChange
names(OE_foldchanges_d20) <- as.character(row.names(sigOE_d20))

cnetplot(ego_d20, 
         categorySize="pvalue", 
         showCategory = 3, 
         color.params = list(foldChange = OE_foldchanges_d20), 
         vertex.label.font=6)
## From Condition as design ##
## Create background dataset for hypergeometric testing using all genes tested for significance in the results                 
allOE_genes <- as.character(row.names(res_c))
## Extract significant results
sigOE <- subset(res_c, padj < 0.05)

sigOE_genes <- as.character(row.names(sigOE))
sigOE_genes = sort(sigOE_genes, decreasing = TRUE)

## GO over-representation analysis
ego <- enrichGO(gene = sigOE_genes, 
                universe = allOE_genes,
                keyType = "ENSEMBL",
                OrgDb = org.Hs.eg.db, 
                ont = "MF", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05,
                readable = TRUE)
## Output results from GO analysis to a table
cluster_summary <- data.frame(ego)
write.csv(cluster_summary, "results/clusterProfiler_res_c.csv")

# barplot 
barplot(ego, 
        drop = TRUE, 
        showCategory = 3, 
        title = "Top 3 GO Molecular Function for All Sample",
        font.size = 10)
# Cnetplot
#To color genes by log2 fold changes, we need to extract the log2 fold changes from our results table creating a named vector
OE_foldchanges <- sigOE$log2FoldChange
names(OE_foldchanges) <- as.character(row.names(sigOE))

cnetplot(ego, 
         categorySize="pvalue", 
         showCategory = 3, 
         color.params = list(foldChange = OE_foldchanges), 
         vertex.label.font=6)

#gseGO clusterprofiler
gse <- gseGO(geneList=sigOE_genes, 
             ont ="ALL", 
             keyType = "ENSEMBL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "none")

#### DESEQ2 Males
## Simple function for plotting a Volcano plot, returns a ggplot object
setwd("~/DESEQ2_R")

### Step 1: import the featurecounts results into R & prepare countdata ###
countdata_male <- read.table("featurecounts_males.txt", header=TRUE, row.names=1)
# Remove first five columns (chr, start, end, strand, length)
countdata_male <- countdata_male[ ,6:ncol(countdata_male)]
# Remove rudementary filename and rename accoring to Sample name
colnames(countdata_male) <- gsub("......results.hisat2.", "", colnames(countdata_male))
colnames(countdata_male) <- gsub(".hisat.sorted.bam", "", colnames(countdata_male))
# Convert to matrix
countdata_male <- as.matrix(countdata_male)
# read in sample info
indexed_male <- read.csv('male_sample.csv', header = TRUE)
sampleinfo_male <- as.matrix(indexed_male[, -1])
row.names(sampleinfo_male) <- indexed_male[, 1]
# check if sample info is matching with count matrix
all(rownames(sampleinfo_male) == colnames(countdata_male))
# run DESEQ2
dds_m <- DESeqDataSetFromMatrix(countData = countdata_male,
                                colData = sampleinfo_male,
                                design = ~ condition)
# pre-filtering: removing rows with low gene counts
# keeping rows that have at least 10 reads total
keep <- rowSums(counts(dds_m)) >= 10
dds_m <- dds_m[keep,]
# set the factor level
dds_m$condition <- relevel(dds_m$condition, ref = "CTRL")
# Run the DESeq pipeline
dds_m <- DESeq(dds_m)
res_m <-results(dds_m, name="condition_INF_vs_CTRL", alpha=0.05)
res_m <- lfcShrink(dds_m, coef="condition_INF_vs_CTRL", type="apeglm", res=res_m)
res_m <- res_m[order(res_m$pvalue),]
summary(res_m)

#visualize  Results
volcano.plot(res = res_m, datasetName = "Infected v/s Control")
# Count data transformations with Regularized log transformation
vsd_m <- vst(dds_m, blind=TRUE)
plotPCA(vsd_m, intgroup="doi")
#annotate Gene Symbol and entrezID & export Results
res_m$symbol <- mapIds(org.Hs.eg.db, keys = rownames(res_m), 
                       keytype = "ENSEMBL", column = "SYMBOL")
res_m$entrez <- mapIds(org.Hs.eg.db, keys = rownames(res_m), 
                       keytype = "ENSEMBL", column = "ENTREZID")
resm_export <- subset(res_m, padj < 0.05)
resm_export <- subset(resm_export, abs(log2FoldChange) >= 1.5)
write.csv(resm_export, 
          file="~/DESEQ2_R/RESULTS/DESEQ2/DEG_all_males.csv")


## venn diagram
# Data preparation
res_c.degs <- subset(res_c, abs(log2FoldChange) >= 1.5)
res_c.degs <- row.names(subset(res_c.degs, padj < 0.05))
res_m.degs <- subset(res_m, abs(log2FoldChange) >= 1.5)
res_m.degs <- row.names(subset(res_m.degs, padj < 0.05))

venn_gender <- venndetail(list(Females = res_c.degs, Males = res_m.degs)) 
plot(venn_gender)
detail(venn_gender)
#extract genes of interest from venn Diagram & save into csv file
venn_gender_shared_genes <- getSet(venn_gender, subset = "Shared")
venn_gender_shared_genes$symbol <- mapIds(org.Hs.eg.db, keys = venn_gender_shared_genes$Detail, 
                                   keytype = "ENSEMBL", column = "SYMBOL")
venn_gender_shared_genes$entrez <- mapIds(org.Hs.eg.db, keys = venn_gender_shared_genes$Detail, 
                                   keytype = "ENSEMBL", column = "ENTREZID")
write.csv(venn_gender_shared_genes, file = "~/DESEQ2_R/RESULTS/venn_gender_shared_genes.csv")

venn_gender_female_genes <- getSet(venn_gender, subset = "Females")
venn_gender_female_genes$symbol <- mapIds(org.Hs.eg.db, keys = venn_gender_female_genes$Detail, 
                                          keytype = "ENSEMBL", column = "SYMBOL")
venn_gender_female_genes$entrez <- mapIds(org.Hs.eg.db, keys = venn_gender_female_genes$Detail, 
                                          keytype = "ENSEMBL", column = "ENTREZID")
write.csv(venn_gender_female_genes, file = "~/DESEQ2_R/RESULTS/venn_gender_female_genes.csv")
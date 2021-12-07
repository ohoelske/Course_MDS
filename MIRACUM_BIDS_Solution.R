###################
# RNA Seq Analyse #

################
# Introduction #

# We ask you here to perform the basic steps of a RNA-Seq analysis. We already provide you with the aligned reads to the human referencr genome, i.e.
# you're starting point are the count data for each sample. Together with the count data we provide you with the annotation table linking sample ids to actuall
# biological conditions and "real world samples". The data set is from a real scientific project we were working on here in Freiburg.
# The data sets contains count data from a single cell line, namely MiaPaCa2, and two biological conditions. The aim of the study is to identify
# the influence of an KLF7 over-expression. Therefore, the cell line was modified to over-express KLF7. As control the untreated cell line was used.
# The over-expressing condition is labeled "*_KLF7" and the control "*_LV". LV stand here for "Leervektor".
# Your task is now to identify the significantly altered genes (in terms of gene expression) between the treated samples (KLF7) and the untreated control (LV).
# Therefore a differentially expressed genes (DEG) analysis has to be performed.
# In the following script we provide you with some hints and already completed commands to do so.

# The following steps are covered

# 1. the count data has to be imported
# 2. linking samples to conditions
# 3. generate a meta table continaing all measured genes and all available gene identifiers
# 4. calculate a principal component analysis (PCA) to get a first overview of the data set
# 5. perform the actual DEG analysis. How many genes are significantly altered up to an adjusted p-value < 0.05; How many of them are up or down regulated respectively?
# 6. visualize and export your results
# 7. generate heatmpas of the highly variable genes

####################
# useful commands #

# head()
# ?functionName or help()
# install.packages()
# accessing columns by name, e.g. result$padj
# filter data.frame, e.g. res[res$test > 5,] the output contains all columns of res but only the rows fulfilling that the value within the "test" column contains a value bigger than 5
# dim()
# length()
# is.na() / !is.na()

################
# useful links #

# links / tutorials / vignettes
# https://www.bioconductor.org/install/
# https://bioconductor.org/packages/release/bioc/html/DESeq2.html
# https://cran.r-project.org/web/packages/FactoMineR/index.html
# https://ggplot2.tidyverse.org ; https://cran.r-project.org/web/packages/ggplot2/index.html
# https://cran.r-project.org/web/packages/ggrepel/vignettes/ggrepel.html
# https://bioconductor.org/packages/release/bioc/html/genefilter.html

#############
# Libraries #

# used packages
# if not installed please install with either bioconductor BiocManager::install() (see useful links) or install.packages()

library(DESeq2)
library(openxlsx)
library(pheatmap)
library(org.Hs.eg.db)
library(ggplot2)
library(genefilter)
library(FactoMineR)
library(ggrepel)
library(apeglm)
library(ggrepel)
library(genefilter)

# Define Working Directory
setwd("~/Projekte/MIRACUM/BIDS/MIRACUM_BIDS_Bioinformatik_Systembiologie_RNA_Seq") # change according to were you put the data

############################################
###                                      ###
###               FUNCTIONS              ###
###                                      ###
############################################

# used to convert different gene identifiers and build the geneMat data.frame
# containing all the different identifiers per gene contained in the data set
# all conversions are based on the org.Hs.eg.db database
# execute the following function before you proceed further

ensembl2entrez <- function(ensembl)
{
  entrez <- mget(as.character(ensembl), org.Hs.egENSEMBL2EG, ifnotfound=NA)
  entrez <- lapply(entrez, function(i) return(i[1]))
  return(unlist(entrez))
}

entrez2ensembl <- function(entrez)
{
  esbl <- mget(as.character(entrez), org.Hs.egENSEMBL, ifnotfound=NA)
  esbl <- lapply(esbl, function(i) return(i[1]))
  return(unlist(esbl))
}

entrez2symbol <- function(entrez)
{
  symbol <- mget(as.character(entrez), org.Hs.egSYMBOL, ifnotfound=NA)
  symbol <- unlist(lapply(symbol, function(i) return(i[1])))
  return(symbol)
}

entrez2genename <- function(entrez)
{
  symbol <- mget(as.character(entrez), org.Hs.egGENENAME, ifnotfound=NA)
  symbol <- unlist(lapply(symbol, function(i) return(i[1])))
  return(symbol)
}

getGeneMat <- function(ensIDs){
  geneMat <- data.frame(ENSEMBL=ensIDs)
  geneMat$ENTREZ <- ensembl2entrez(geneMat$ENSEMBL)
  idxNA <- !is.na(geneMat$ENTREZ)
  sym <- entrez2symbol(na.omit(geneMat$ENTREZ))
  genename <- entrez2genename(na.omit(geneMat$ENTREZ))
  geneMat$Symbol <- NA
  geneMat$Symbol[idxNA] <- sym
  geneMat$Genename <- NA
  geneMat$Genename[idxNA] <- genename
  rownames(geneMat) <- geneMat$ENSEMBL
  return(geneMat)
}

########
# MAIN #

##### 
# Provide

## Importing and creating the count matrix together with the gene matrix
## Create Count Matrix

# Load /import count data
targets <- read.xlsx(xlsxFile = "targets.xlsx")
countFiles <- list.files(pattern = "ReadsPerGene.out.tab")
names(countFiles) <- gsub("_ReadsPerGene.out.tab", "", countFiles)
countList <- lapply(countFiles, read.delim, skip = 4, header = FALSE, stringsAsFactors = FALSE)
names(countList) <- gsub("_ReadsPerGene.out.tab", "", targets$ID)
rownames(targets) <- targets$ID

# obtain gene identifier
ensIDs <- countList[[1]]$V1

# Make count table
countList <- lapply(countList, function(i) return(i[,2]))
countMat <- do.call(cbind, countList)
rownames(countMat) <- ensIDs

# Create geneMat from all genes contained in the experiment (measured)
geneMat <- getGeneMat(ensIDs = ensIDs)

# save count matrix, geneMat and targets; export count matrix to TXT file

write.table(countMat, "counts.txt", sep = "\t", quote = FALSE)
save(countMat, geneMat, targets, file = 'countData.RData')

## Building DESeq2 model as preparation to calculate differentially expressed genes
# We are interested in the comparison between KLF7 over-expressing cells and cells containing only the empty construct
# Define Comparisons, Design, etc.

targets$Group <- factor(targets$Group, levels = c("MiaPaCa2_LV", "MiaPaCa2_KLF7"))

# build DESeq2 model
dds <- DESeqDataSetFromMatrix(countData = countMat,
                              colData = targets,
                              design = ~ Group)

# add gene meta data
mcols(dds) <- DataFrame(mcols(dds), geneMat)

# filtering low expressed genes
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]

# calculate and get normalized count matrix
dds <- estimateSizeFactors(dds)
dds.counts <- counts(dds, normalized=TRUE)

# get normalized rlogs
rld <- rlogTransformation(dds, blind=TRUE)
rlds <- assay(rld)

#####
# homework

# Plot PCA
# calculate a principal component analysis (PCA) of the expression set with the samples as the individuals and visualize the results in a "pretty manner"
# you could use the FactoMineR R package for calculation and the ggplot2 package for plotting
# save the resulting pca as a pdf. The ggsave function from the ggplot2 package could be helpful.
# #"quick and dirty"
# DESeq2::plotPCA(rld, intgroup = 'Group')

# PCA
pca <- PCA(t(rlds), graph = F, scale.unit = F)
scores <- as.data.frame(pca$ind$coord)
scores$Condition <- targets$Group
scores$labels <- targets$Sample

pp <- ggplot(data = scores, aes(x=Dim.1,y=Dim.2,color=Condition,label=labels)) + geom_point(size = 4) + theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
pp <- pp + labs(title = 'PCA', x = paste('PC1','(',round(pca$eig[1,2],2),'%)'), y = paste('PC2','(',round(pca$eig[2,2],2),'%)'))
pp
ggsave(pp, filename = "PCA.pdf", device = "pdf", dpi = "retina")

## differential analysis
# perform the actual differential expression analysis
# use the "Wald" test and fitType "mean"
# use also the lfcShrink function with type "apeglm" for more precise log fold-change (logFC, LFC) values
# How many DEGs to we have in total (padj < 0.05)
# How many are up-regulated and how many are down-regulated?

dds <- DESeq(dds, test = "Wald", fitType = "mean", parallel = T)
res <- results(object = dds, alpha = 0.05, parallel = T)
# get the output names of the result
resultsNames(dds)
# Re-estimating logFCs
resLFC <- lfcShrink(dds = dds, coef="Group_MiaPaCa2_KLF7_vs_MiaPaCa2_LV", type="apeglm", parallel = T, res = res)
# get a summary of the results
summary(resLFC)

# out of 24469 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 243, 0.99%
# LFC < 0 (down)     : 257, 1.1%
# outliers [1]       : 10, 0.041%
# low counts [2]     : 2372, 9.7%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# export results
# export the results of the analysis including the gene meta information from the geneMat matrix
# think about which identifier are common between your DEG result and the geneMat. the funtion merge() could be useful
# write the results to an excel table; in the R package openxlsx the function write.xlsx can be used
resLFC$Ensembl_ID <- rownames(resLFC)
results.deseq2 <- as.data.frame(resLFC)
results.deseq2 <- merge(results.deseq2, geneMat, by.x = "Ensembl_ID", by.y = "ENSEMBL")

write.xlsx(results.deseq2, file = "MiaPaCa2_KLF7vsLV_DESeq2.xlsx", asTable = F, firstRow = T, headerStyle = createStyle(textDecoration = 'bold'), keepNA = F, row.names = F)

# create an additional excel table reporting only the significant DEGs
# sort them from highest significance to lowest; hint: function order()
# basically filter your results based on the padj column, be careful with NA values in the padj column; helpful na.omit(), is.na or !is.na()
results.deseq2.tmp <- results.deseq2[!is.na(results.deseq2$padj),]
results.deseq2.sig <- results.deseq2.tmp[results.deseq2.tmp$padj < 0.05,]
results.deseq2.sig.sort <- results.deseq2.sig[order(results.deseq2.sig$padj, decreasing = F),]
write.xlsx(results.deseq2.sig.sort, file = "MiaPaCa2_KLF7vsLV_sig_DEGs_DESeq2.xlsx", asTable = F, firstRow = T, headerStyle = createStyle(textDecoration = 'bold'), keepNA = F, row.names = F, overwrite = T)


## volcano plot
# create a volcano plot out of the DEG results and add meaningful axes names and plot title
# use again padj, i.e. the adjusted or FDR corrected p-value, < 0.05 as significance cutoff for highlighting significant genes but use all genes for plotting
# reminder: volcano plot: plot logFC against significance 
# to improve readability you can transform the p-values to a log10 (function log10()) scale to get rid of the small numbers
# the simplest way of creating the plot is to generate a new data.frame
# tmp <- data.frame(ColumnName1 = values, ColumnName2 = values, etc.); e.g. Symbol as gene identifier, logFC from the log2FoldChange column and pValues from the padj column
# changing ggplot labels/title etc. labs could be useful
# store plot as pdf
# optional: add the most significant genes as labels into the plot; decide for a meaningful cutoff so that you don't have a too crowded plot; together with a p-value cutoff an additional logFC cutoff might help
# useful package is ggrepel with the function geom_label_repel

df <- data.frame(ID=results.deseq2$Symbol,
                 logFC = results.deseq2$log2FoldChange,
                 adj.P.Val = -log10(results.deseq2$padj),
                 Condition = "KLF7_LV")

df$Significance <- "adj. p-value => 0.05"
df$Significance[df$adj.P.Val >= -log10(0.05)] <- "adj. p-value < 0.05"
df$Significance <- factor(x = df$Significance, levels = c("adj. p-value < 0.05","adj. p-value => 0.05"))

# plot
p <- ggplot(df, aes(x=logFC, y=adj.P.Val, color=Significance, label=ID)) + geom_point(size=0.75) + scale_color_manual(values=c('red','black'))
p <- p + labs( title="DEGs", x = "log2FC", y = "-log10(adj. p-value)") + theme_bw()
p <- p + geom_label_repel(aes(logFC, adj.P.Val, label=ifelse((adj.P.Val>=-log10(1e-5) & abs(logFC) > 2), as.character(ID), '')), point.padding=0.5, segment.color= "grey", segment.size= 0.5, size=4)
p
ggsave(p, filename = "Volcano_MiaPaCa2_KLF7vsLV.pdf", device = "pdf", dpi = "retina")

## identify and visualized the genes with the highest variance
## Another interesting way the dig into the data is the visualize highly variable genes as a heatmap
## identifying highly variable genes can be done with the genefilter package and the function varFilter();
## possible var.func is IQR (interquartile range) with filtering by quantiles to filter for a var.cutoff of 0.99, i.e. only to top 1 % o the genes with the highest IQR are kept
## the resulting genes should be visualized as normalized expressions (contained in the rlds variable) and separately as z-scores; the function genescale could be used to calculate z-scores out of the rlds object;
## reminder: z-scores; zero mean and a standard deviation of 1 per gene;
## second reminder: genes are contained in rows ;)
## for plotting the pheatmap package could be used
## usually a colormap from green over black to red is used for normalized expression and a colorscale from blue to yellow is just for z-scores
## optional: color the two groups, i.e. KLF7 over-expression and LV (empty vector)
## the parameter annotation_col in the pheatmap function is useful for this part
## save the plot as pdf; check the pheatmap parameters ;)

tmp <- varFilter(eset = rlds, var.func=IQR, var.cutoff=0.99, filterByQuantile=TRUE)
zscore <- genescale(tmp, axis = 1, method = "Z")

# column annotation
column.annotation <- data.frame(Condition = targets$Group)
rownames(column.annotation) <- targets$ID

# normalized expression
pheatmap(tmp, show_rownames = F, color = colorRampPalette(c("green","black","red"))(50), annotation_col = column.annotation, filename = "Highly_variable_genes_normalized_expression.pdf")

# z-scoress
pheatmap(zscore, show_rownames = F, color = colorRampPalette(c("blue","yellow"))(50), annotation_col = column.annotation, filename = "Highly_variable_genes_z_score.pdf")



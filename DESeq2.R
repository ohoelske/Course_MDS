
###### Excluding low/ none expressed genes
keep <- rowSums(countMat > 0) >= 3 

# BUILD DESEQ DATA
dds <- DESeqDataSetFromMatrix(countData = countMat[keep, ],
                              colData = ann,
                              design = ~ condition)
mcols(dds) <- DataFrame(mcols(dds), geneMat[keep, ])

# GET COUNT
dds <- estimateSizeFactors(dds)
dds.counts <- counts(dds, normalized=TRUE)

# DIFFERENTIAL ANALYSIS
dds <- DESeq(dds, test = "Wald", fitType = "mean")
result_dss = results(dds)


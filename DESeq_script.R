library(DESeq2)
library(EnhancedVolcano)

# load data
counts_data <- read.csv('counts_combined.csv', row.names = "X")
column_data <- read.csv('column_data.csv', row.names = "X")

# DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = counts_data,
                       colData = column_data,
                       design = ~ status)

# optional filtering of genes with low counts
keep <- rowSums(counts(dds)) >= 0
dds <- dds[keep,]

# set factor level e.g. status: original vs perturbed cells
dds$status <- relevel(dds$status, ref = "original")

# run DESeq
dds <- DESeq(dds)

# results
res <- results(dds)
summary(res)

# plotting
plotMA(res)

# volcano plot
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                pointSize = 2.0,
                labSize = 2.0)



# BiocManager::install("DESeq2")
# BiocManager::install("apeglm")
# BiocManager::install("EnhancedVolcano")

library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)

dir <- "C:/Users/User/Desktop/BI/Practicum/NO_project/215213_Counts"
ff <- list.files(path = dir, pattern = "*ReadsPerGene.out.tab$", full.names = TRUE ) 
counts.files <- lapply(ff, read.table, skip = 4)
counts <- as.data.frame(sapply(counts.files, function(x) x[ , 2]))
ff <- gsub("C:/Users/User/Desktop/BI/Practicum/NO_project/215213_counts/", "", ff)
ff <- gsub("ReadsPerGene.out.tab", "", ff)
colnames(counts) <- ff 
row.names(counts) <- counts.files[[1]]$V1
View(counts)
write.csv(counts, "C:/Users/User/Desktop/BI/Practicum/NO_project/215213_counts/counts_215213.csv")

matrix <- read.table(file = "C:/Users/User/Desktop/BI/Practicum/NO_project/215213_counts/counts_215213.csv", header = TRUE, sep = ",")

#variant 1: to sum technical replicates
matrix$SRR9703414_15 <- matrix$SRR9703414 + matrix$SRR9703415
matrix$SRR9703416_17 <- matrix$SRR9703416 + matrix$SRR9703417
matrix$SRR9703418_19 <- matrix$SRR9703418 + matrix$SRR9703419
matrix$SRR9703420_21 <- matrix$SRR9703420 + matrix$SRR9703421
matrix$SRR9703422_23 <- matrix$SRR9703422 + matrix$SRR9703423
matrix$SRR9703424_25 <- matrix$SRR9703424 + matrix$SRR9703425
matrix$SRR9703426_27 <- matrix$SRR9703426 + matrix$SRR9703427
matrix$SRR9703428_29 <- matrix$SRR9703428 + matrix$SRR9703429
matrix$SRR9703430_31 <- matrix$SRR9703430 + matrix$SRR9703431
matrix$SRR9703432_33 <- matrix$SRR9703432 + matrix$SRR9703433
matrix <- matrix[ ,-c(14:33)]

#variant 2: to find out average between technical replicates
matrix$SRR9703414_15 <- round((matrix$SRR9703414 + matrix$SRR9703415)/2, 0)
matrix$SRR9703416_17 <- round((matrix$SRR9703416 + matrix$SRR9703417)/2, 0)
matrix$SRR9703418_19 <- round((matrix$SRR9703418 + matrix$SRR9703419)/2, 0)
matrix$SRR9703420_21 <- round((matrix$SRR9703420 + matrix$SRR9703421)/2, 0)
matrix$SRR9703422_23 <- round((matrix$SRR9703422 + matrix$SRR9703423)/2, 0)
matrix$SRR9703424_25 <- round((matrix$SRR9703424 + matrix$SRR9703425)/2, 0)
matrix$SRR9703426_27 <- round((matrix$SRR9703426 + matrix$SRR9703427)/2, 0)
matrix$SRR9703428_29 <- round((matrix$SRR9703428 + matrix$SRR9703429)/2, 0)
matrix$SRR9703430_31 <- round((matrix$SRR9703430 + matrix$SRR9703431)/2, 0)
matrix$SRR9703432_33 <- round((matrix$SRR9703432 + matrix$SRR9703433)/2, 0)
matrix <- matrix[ ,-c(14:33)]

Counts <- as.matrix(matrix[,-1]) 
colnames(Counts) <- names(matrix)[-1]
rownames(Counts) <- matrix$X
Sums <-rowSums(Counts)
Counts_info <- Counts[Sums >=2,]

coldata <- data.frame(genotype = c(rep('control',12), rep('MSA',10)))
rownames(coldata) <- colnames(Counts_info)
dds <- DESeqDataSetFromMatrix(countData = Counts_info, colData = coldata, design = ~ genotype)
for_pca <- rlog(dds)
plotPCA(object = for_pca, intgroup = 'genotype') + geom_text(aes(label=name),vjust=1,check_overlap = TRUE,size = 2)
DataSetAnalysis <- DESeq(dds, test = 'Wald', sfType = 'ratio')
res <- results(DataSetAnalysis, alpha = 0.05, lfcThreshold = 1)
resLFC2 <- lfcShrink(DataSetAnalysis, coef = "genotype_MSA_vs_control", type ="apeglm", res = res, lfcThreshold = 1)
plotMA(resLFC2)
EnhancedVolcano(resLFC2, lab = row.names(resLFC2), x = 'log2FoldChange', y = 'svalue', pCutoff = 0.05, FCcutoff = 1.0, titleLabSize = 8, subtitleLabSize = 8, captionLabSize = 8, axisLabSize = 8, legendLabSize = 8, legendIconSize = 1.0, legendPosition = 'none')

#Enhanded volcano with NO genes altering only
lab_italics <- paste0("italic('", rownames(resLFC2), "')")
selectLab_italics = paste0("italic('",
  +     c('CD36','AQP1','SPR','CCL2','TLR2','IFNG'),
  +     "')")
EnhancedVolcano(resLFC2, lab = lab_italics, x = 'log2FoldChange', y = 'svalue', selectLab = selectLab_italics, title = 'Multiple system atrophy, putamen', subtitle = ' ', caption = ' ', pCutoff = 0.005, FCcutoff = 1.0, titleLabSize = 25, axisLabSize = 20, pointSize = 3.0,labSize = 8.0,  drawConnectors = TRUE, widthConnectors = 1.0, colAlpha = 2/5, legendPosition = 'none', colConnectors = 'black', boxedLabels = TRUE, parseLabels = TRUE)

#Save results of differential gene expression
sorted = res[with(res, order(padj, -log2FoldChange)), ]

sorted.df = data.frame("id"=rownames(sorted),sorted)
write.table(sorted.df, file="result213215_replicates_as_sum.txt", sep="\t", col.names=NA, quote=FALSE)
#Save corrected results
sorted = resLFC2[with(resLFC2, order(svalue, -log2FoldChange)), ]
sorted.df = data.frame("id"=rownames(sorted),sorted)
genes = subset(sorted.df, svalue<0.005)
#save normalized counts 
nc = counts(DataSetAnalysis, normalized = TRUE)
dt = data.frame("id"=rownames(nc), nc)
write.table(dt, file="counts_norm_213215_replicates_as_sum.txt", sep="\t", col.names=NA, quote=FALSE)

#For 058181 dataset
matrix <- read.table(file = "C:/Users/User/Desktop/BI/Practicum/NO_project/058181_counts/counts_058181.csv", header = TRUE, sep = ",")
Counts <- as.matrix(matrix[,-1]) 
colnames(Counts) <- names(matrix)[-1]
rownames(Counts) <- matrix$X
Sums <-rowSums(Counts)
Counts_info <- Counts[Sums >=2,]
coldata <- data.frame(genotype = c(rep('control',42), rep('PD',29)))
rownames(coldata) <- colnames(Counts_info)
dds <- DESeqDataSetFromMatrix(countData = Counts_info, colData = coldata, design = ~ genotype)
DataSetAnalysis <- DESeq(dds, test = 'Wald', sfType = 'ratio')
res <- results(DataSetAnalysis, alpha = 0.05, lfcThreshold = 1)
resLFC2 <- lfcShrink(DataSetAnalysis, coef = "genotype_PD_vs_control", type ="apeglm", res = res, lfcThreshold = 1)

EnhancedVolcano(resLFC2, lab = lab_italics_058181, x = 'log2FoldChange', y = 'svalue', selectLab = selectLab_italics_058181, title = 'Parkinson disease, prefrontal cortex', subtitle = ' ', caption = ' ', pCutoff = 0.005, FCcutoff = 1.0, titleLabSize = 25, axisLabSize = 20, pointSize = 3.0,labSize = 8.0,  drawConnectors = TRUE, widthConnectors = 1.0, colAlpha = 2/5, legendPosition = 'none', colConnectors = 'black', boxedLabels = TRUE, parseLabels = TRUE)
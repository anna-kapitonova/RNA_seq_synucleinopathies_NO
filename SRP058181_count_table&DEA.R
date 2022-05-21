# BiocManager::install("DESeq2")
# BiocManager::install("apeglm")
# BiocManager::install("EnhancedVolcano")

library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)

#Constructs table with counts

#indicates directory with files containing counts from STAR alignment
dir <- "C:/Users/User/Desktop/BI/Practicum/NO_project/058181_Counts"

#collects all counts into one table, skipping four first lines (they contain technical info) and selecting second column of each file (because the library is non-stranded)
ff <- list.files(path = dir, pattern = "*ReadsPerGene.out.tab$", full.names = TRUE ) 
counts.files <- lapply(ff, read.table, skip = 4)
counts <- as.data.frame(sapply(counts.files, function(x) x[ , 2]))

#renames columns according to sample name
ff <- gsub("C:/Users/User/Desktop/BI/Practicum/NO_project/058181_Counts/", "", ff)
ff <- gsub("ReadsPerGene.out.tab", "", ff)
colnames(counts) <- ff 

#renames rows according to genes symbols written in first column of each file
row.names(counts) <- counts.files[[1]]$V1

#writes obtained table into file
write.csv(counts, "C:/Users/User/Desktop/BI/Practicum/NO_project/058181_Counts/counts_058181.csv")

#Differential expression analysis with DESeq2

#Creates table with counts from file
matrix <- read.table(file = "C:/Users/User/Desktop/BI/Practicum/NO_project/058181_Counts/counts_058181.csv", header = TRUE, sep = ",")


#creates matrix from table with counts
Counts <- as.matrix(matrix[,-1]) 
colnames(Counts) <- names(matrix)[-1]
rownames(Counts) <- matrix$X

#Deletes rows where the sum of counts is less than 2
Sums <-rowSums(Counts)
Counts_info <- Counts[Sums >=2,]

#creates table with experimental design
coldata <- data.frame(genotype = c(rep('control',42), rep('PD',29)))
rownames(coldata) <- colnames(Counts_info)

#creates DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = Counts_info, colData = coldata, design = ~ genotype)

#PCA
for_pca <- rlog(dds)
plotPCA(object = for_pca, intgroup = 'genotype') + geom_text(aes(label=name),vjust=1,check_overlap = TRUE,size = 2)

#Analysis of differential expression
DataSetAnalysis <- DESeq(dds, test = 'Wald', sfType = 'ratio')
res <- results(DataSetAnalysis, alpha = 0.05, lfcThreshold = 1)

#Log2FoldChange apeglm correction, lfc_threshold is > 1 or <-1
resLFC2 <- lfcShrink(DataSetAnalysis, coef = "genotype_PD_vs_control", type ="apeglm", res = res, lfcThreshold = 1)

#MA_plot
plotMA(resLFC2)

#Volcano-plot
EnhancedVolcano(resLFC2, 
                lab = row.names(resLFC2), 
                x = 'log2FoldChange', 
                y = 'svalue',
                title = ' ', 
                subtitle = ' ', 
                ylab = bquote(~Log[10]~ 'S'),
                caption = ' ',
                pCutoff = 0.05, 
                FCcutoff = 1.0, 
                titleLabSize = 8, 
                subtitleLabSize = 8, 
                captionLabSize = 8, 
                axisLabSize = 8, 
                legendLabSize = 8, 
                legendIconSize = 1.0, 
                legendPosition = 'none')


#Saves genes with significant results(lfc corrected <0.005)
sorted = resLFC2[with(resLFC2, order(svalue, -log2FoldChange)), ]
sorted.df = data.frame("id"=rownames(sorted),sorted)
genes = subset(sorted.df, svalue<0.005)
write.table(genes, file="DE_genes_total_085181", sep="\t", col.names=NA, quote=FALSE)

#Saves upregulated and downregulated genes
genes_up <- subset(genes, log2FoldChange > 0)
write.table(genes_up, file="genes_up_085181.txt", sep="\t", col.names=NA, quote=FALSE)
genes_down <- subset(genes, log2FoldChange < 0)
write.table(genes_down, file="genes_down_085181.txt", sep="\t", col.names=NA, quote=FALSE)

#saves normalized counts 
nc = counts(DataSetAnalysis, normalized = TRUE)
dt = data.frame("id"=rownames(nc), nc)
write.table(dt, file="counts_norm_058181.txt", sep="\t", col.names=NA, quote=FALSE)

#Enhanced volcano with NO genes altering only
lab_italics <- paste0("italic('", rownames(resLFC2), "')")
selectLab_italics = paste0("italic('", c('CD36','AQP1','SPR','CCL2','TLR2','IFNG'),"')")
EnhancedVolcano(resLFC, 
                lab = lab_italics, 
                x = 'log2FoldChange', 
                y = 'svalue',
                ylab = bquote(~Log[10]~ 'S'),
                selectLab = selectLab_italics, 
                title = 'Multiple system atrophy, putamen', 
                subtitle = ' ', 
                caption = ' ', 
                pCutoff = 0.005, 
                FCcutoff = 1.0, 
                titleLabSize = 25, 
                axisLabSize = 20, 
                pointSize = 3.0,
                labSize = 8.0,  
                drawConnectors = TRUE, 
                widthConnectors = 1.0, 
                colAlpha = 2/5, 
                legendPosition = 'none', 
                colConnectors = 'black', 
                boxedLabels = TRUE, 
                parseLabels = TRUE)


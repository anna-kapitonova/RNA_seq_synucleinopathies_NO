# BiocManager::install("DESeq2")
# BiocManager::install("apeglm")
# BiocManager::install("EnhancedVolcano")

#### loading of the packages (if needed, install first)
library(DESeq2)
library(dplyr)
library(tidyr)
library(apeglm)
library(VennDiagram)
library(EnhancedVolcano)
library(RColorBrewer)

###############################################################
#### Create a table with counts from ReadsPerGene.out.tab files
###############################################################

dir <- "/Users/anna/Desktop/BI/NO/results/SRP324001/alignments" # save the path to project directory
ff <- list.files(path = dir, pattern = "*ReadsPerGene.out.tab$", full.names = TRUE ) # find all files with gene counts
counts.files <- lapply(ff, read.table, skip = 4) # skip first 4 rows with general info
counts <- as.data.frame(sapply(counts.files, function(x) x[ , 2])) # we use the second column with counts for unstranded reads
ff <- gsub("/Users/anna/Desktop/BI/NO/results/SRP324001/alignments/", "", ff) # create column names from file names
ff <- gsub("ReadsPerGene.out.tab", "", ff)
colnames(counts) <- ff # columns are sample names
row.names(counts) <- counts.files[[1]]$V1 # rows are gene names
write.csv(counts, "/Users/anna/Desktop/BI/NO/results/SRP324001/alignments/counts.csv") # write the table into csv file


#################################################
#### Differential expression analysis with DESeq2
#################################################

# set up the conditions based on the experimental setup
cond <- c("сontrol", "сontrol", "DLB", "DLB", "DLB", "DLB", "DLB",
         "PDD", "PDD", "PD", "PDD", "PD", "PDD", "PD", "PDD", "PDD",
         "DLB", "PD", "PDD", "PD", "DLB", "PD", "PD", "сontrol",
         "сontrol", "сontrol", "сontrol", "сontrol")

# build the dataframe from the conditions
samples <- names(counts)
condition <- factor(c(cond))
colData <- data.frame(samples=samples, condition=condition)
levels(colData$condition) <- c('control','PD', 'PDD', 'DLB')

# remove rows with almost absent counts
Sums <- rowSums(counts)
table(Sums >= 2)
counts <- counts[Sums >= 2, ]

# create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = colData,
                              design = ~condition)

# set the reference to be compared
dds$condition <- relevel(dds$condition,"control")

# principal component analysis of samples
for_pca <- rlog(dds)
plotPCA(object = for_pca, intgroup = 'condition') +
  theme(text = element_text(size = 17))

# run analysis
dds <- DESeq(dds, test = 'Wald', sfType = 'ratio')

#### Normalized data matrix

# get normalized counts and write this to a file
nc = counts(dds, normalized=TRUE)

# turn it into a dataframe to have proper column names
dt = data.frame("id" = rownames(nc), nc)

# save the normalized data matrix
write.table(dt, file = "norm-matrix.txt",
            sep = "\t", row.names = FALSE,
            col.names = TRUE, quote = FALSE)


# format the results
res_PD <- results(dds, alpha = 0.05, lfcThreshold = 1,
                  contrast = c("condition", "PD", "control"))
res_PDD <- results(dds, alpha = 0.05, lfcThreshold = 1,
                   contrast = c("condition", "PDD", "control"))
res_DLB <- results(dds, alpha = 0.05, lfcThreshold = 1,
                   contrast = c("condition", "DLB", "control"))
summary(res_PD)
summary(res_PDD)
summary(res_DLB)

# logfoldchange correction and MA plots
resLFC2_PD <- lfcShrink(dds, coef = "condition_PD_vs_control",
                        type="apeglm", res = res_PD, lfcThreshold = 1)
resLFC2_PDD <- lfcShrink(dds, coef = "condition_PDD_vs_control",
                         type="apeglm", res = res_PDD, lfcThreshold = 1)
resLFC2_DLB <- lfcShrink(dds, coef = "condition_DLB_vs_control",
                         type="apeglm", res = res_DLB, lfcThreshold = 1)
summary(resLFC2_PD)
summary(resLFC2_PDD)
summary(resLFC2_DLB)

plotMA(resLFC2_PD)
plotMA(resLFC2_PDD)
plotMA(resLFC2_DLB)

# volcano plot
EnhancedVolcano(resLFC2_PD, lab = row.names(resLFC2_PD),
                x = 'log2FoldChange', y = 'svalue',
                pCutoff = 0.005, FCcutoff = 1.0,
                title = "",
                subtitle = "", caption = "",
                ylab = expression("-Log"[10]*"S"),
                pointSize = 3,
                labSize = 5,
                axisLabSize = 20,
                legendPosition = 'none')

EnhancedVolcano(resLFC2_PDD, lab = row.names(resLFC2_PDD),
                x = 'log2FoldChange', y = 'svalue',
                pCutoff = 0.005, FCcutoff = 1.0,
                title = "",
                subtitle = "", caption = "",
                ylab = expression("-Log"[10]*"S"),
                pointSize = 3,
                labSize = 5,
                axisLabSize = 20,
                legendPosition = 'none')

EnhancedVolcano(resLFC2_DLB, lab = row.names(resLFC2_DLB),
                x = 'log2FoldChange', y = 'svalue',
                pCutoff = 0.005, FCcutoff = 1.0,
                title = "",
                subtitle = "", caption = "",
                ylab = expression("-Log"[10]*"S"),
                pointSize = 3,
                labSize = 5,
                axisLabSize = 20,
                legendPosition = 'none')

#### Save the results

# sort the results dataframe by the padj(or svalue) and LFC columns
sorted_PD <- resLFC2_PD[with(resLFC2_PD, order(svalue, -log2FoldChange)), ]
sorted_PDD <- resLFC2_PDD[with(resLFC2_PDD, order(svalue, -log2FoldChange)), ]
sorted_DLB <- resLFC2_DLB[with(resLFC2_DLB, order(svalue, -log2FoldChange)), ]

# turn it into a dataframe to have proper column names
sorted_PD.df <- data.frame("id" = rownames(sorted_PD), sorted_PD)
sorted_PDD.df <<- data.frame("id" = rownames(sorted_PDD), sorted_PDD)
sorted_DLB.df = data.frame("id" = rownames(sorted_DLB), sorted_DLB)

# write the table out
write.table(sorted_PD.df,
            file = "/Users/anna/Desktop/BI/NO/results/SRP324001/result_PD_corrected.txt",
            sep = "\t", col.names = NA, quote = FALSE)
write.table(sorted_PDD.df,
            file = "/Users/anna/Desktop/BI/NO/results/SRP324001/result_PDD_corrected.txt",
            sep = "\t", col.names = NA, quote = FALSE)
write.table(sorted_DLB.df,
            file = "/Users/anna/Desktop/BI/NO/results/SRP324001/result_DLB_corrected.txt",
            sep = "\t", col.names = NA, quote = FALSE)

#### Venn diagrams with up- and downregulated genes
res_PD_noNA <- resLFC2_PD[complete.cases(resLFC2_PD),]
down_PD <- rownames(res_PD_noNA[res_PD_noNA$svalue < 0.005 & res_PD_noNA$log2FoldChange < -1, ])
up_PD <- rownames(res_PD_noNA[res_PD_noNA$svalue < 0.005 & res_PD_noNA$log2FoldChange > 1, ])

venn.diagram(x = list(rownames(resLFC2_PD), down_PD, up_PD),
             category.names = c("All", "Down", "Up"),
             filename = "/Users/anna/Desktop/BI/NO/results/SRP324001/plots/venn_PD_corrected.png",
             imagetype = "png",
             cex = 2, cat.cex = 3,
             fill = c("#B3E2CD", "#FDCDAC", "#96AAD1"),
             fontfamily = "sans",
             cat.fontfamily = "sans")

res_PDD_noNA <- resLFC2_PDD[complete.cases(resLFC2_PDD),]
down_PDD <- rownames(res_PDD_noNA[res_PDD_noNA$svalue < 0.005 & res_PDD_noNA$log2FoldChange < -1, ])
up_PDD <- rownames(res_PDD_noNA[res_PDD_noNA$svalue < 0.005 & res_PDD_noNA$log2FoldChange > 1, ])

venn.diagram(x = list(rownames(resLFC2_PDD), down_PDD, up_PDD),
             category.names = c("All", "Down", "Up"),
             filename = "/Users/anna/Desktop/BI/NO/results/SRP324001/plots/venn_PDD_corrected.png",
             imagetype = "png",
             cex = 2, cat.cex = 2,
             fill = c("#B3E2CD", "#FDCDAC", "#96AAD1"),
             fontfamily = "sans",
             cat.fontfamily = "sans")

res_DLB_noNA <- resLFC2_DLB[complete.cases(resLFC2_DLB),]
down_DLB <- rownames(res_DLB_noNA[res_DLB_noNA$svalue < 0.005 & res_DLB_noNA$log2FoldChange < -1, ])
up_DLB <- rownames(res_DLB_noNA[res_DLB_noNA$svalue < 0.005 & res_DLB_noNA$log2FoldChange > 1, ])

venn.diagram(x = list(rownames(resLFC2_DLB), down_DLB, up_DLB),
             category.names = c("All", "Down", "Up"),
             filename = "/Users/anna/Desktop/BI/NO/results/SRP324001/plots/venn_DLB_corrected.png",
             imagetype = "png",
             cex = 2, cat.cex = 3,
             fill = c("#B3E2CD", "#FDCDAC", "#96AAD1"),
             fontfamily = "sans",
             cat.fontfamily = "sans")

# extract lists of genes
write.table(c(up_PD, down_PD), file = "/Users/anna/Desktop/BI/NO/results/SRP324001/DEG_PD_corrected.txt",
            row.names = FALSE, quote = FALSE, col.names = FALSE)
write.table(down_PD, file = "/Users/anna/Desktop/BI/NO/results/SRP324001/genes_down_PD_corrected.txt",
            row.names = FALSE, quote = FALSE, col.names = FALSE)
write.table(up_PD, file = "/Users/anna/Desktop/BI/NO/results/SRP324001/genes_up_PD_corrected.txt",
            row.names = FALSE, quote = FALSE, col.names = FALSE)

write.table(c(up_PDD, down_PDD), file = "/Users/anna/Desktop/BI/NO/results/SRP324001/DEG_PDD_corrected.txt",
            row.names = FALSE, quote = FALSE, col.names = FALSE)
write.table(down_PDD, file = "/Users/anna/Desktop/BI/NO/results/SRP324001/genes_down_PDD_corrected.txt",
            row.names = FALSE, quote = FALSE, col.names = FALSE)
write.table(up_PDD, file = "/Users/anna/Desktop/BI/NO/results/SRP324001/genes_up_PDD_corrected.txt",
            row.names = FALSE, quote = FALSE, col.names = FALSE)

write.table(c(up_DLB, down_DLB), file = "/Users/anna/Desktop/BI/NO/results/SRP324001/DEG_DLB_corrected.txt",
            row.names = FALSE, quote = FALSE, col.names = FALSE)
write.table(down_DLB, file = "/Users/anna/Desktop/BI/NO/results/SRP324001/genes_down_DLB_corrected.txt",
            row.names = FALSE, quote = FALSE, col.names = FALSE)
write.table(up_DLB, file = "/Users/anna/Desktop/BI/NO/results/SRP324001/genes_up_DLB_corrected.txt",
            row.names = FALSE, quote = FALSE, col.names = FALSE)
# BiocManager::install("DESeq2")
# BiocManager::install("apeglm")
# BiocManager::install("EnhancedVolcano")
# BiocManager::install("clusterProfiler")


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

dir <- "/Users/anna/Desktop/BI/NO/results/SRP148970/alignments" # save the path to project directory
ff <- list.files(path = dir, pattern = "*ReadsPerGene.out.tab$", full.names = TRUE ) # find all files with gene counts
counts.files <- lapply(ff, read.table, skip = 4) # skip first 4 rows with general info
counts <- as.data.frame(sapply(counts.files, function(x) x[ , 2])) # we use the second column with counts for unstranded reads
ff <- gsub("/Users/anna/Desktop/BI/NO/results/SRP148970/alignments/", "", ff) # create column names from file names
ff <- gsub("ReadsPerGene.out.tab", "", ff)
colnames(counts) <- ff # columns are sample names
row.names(counts) <- counts.files[[1]]$V1 # rows are gene names

# averaging technical replicates
counts['SRR7218297_1'] <- (counts['SRR7218297_1'] + counts['SRR7218298_1']) %/% 2
counts['SRR7218301_1'] <- (counts['SRR7218301_1'] + counts['SRR7218302_1'] + counts['SRR7218303_1']) %/% 3
counts['SRR7218311_1'] <- (counts['SRR7218311_1'] +
                             counts['SRR7218312_1'] +
                             counts['SRR7218313_1'] +
                             counts['SRR7218314_1'] +
                             counts['SRR7218315_1']) %/% 5
counts['SRR7218316_1'] <- (counts['SRR7218316_1'] +
                             counts['SRR7218317_1'] +
                             counts['SRR7218318_1']) %/% 3
counts['SRR7218319_1'] <- (counts['SRR7218319_1'] + counts['SRR7218320_1']) %/% 2
counts['SRR7218324_1'] <- (counts['SRR7218324_1'] + counts['SRR7218325_1']) %/% 2
counts['SRR7218326_1'] <- (counts['SRR7218326_1'] + counts['SRR7218327_1']) %/% 2
counts['SRR7218331_1'] <- (counts['SRR7218331_1'] + counts['SRR7218332_1']) %/% 2

# or you can use summarization strategy

# counts['SRR7218297_1'] <- counts['SRR7218297_1'] + counts['SRR7218298_1']
# counts['SRR7218301_1'] <- counts['SRR7218301_1'] + counts['SRR7218302_1'] + counts['SRR7218303_1']
# counts['SRR7218311_1'] <- counts['SRR7218311_1'] +
  # counts['SRR7218312_1'] +
  # counts['SRR7218313_1'] +
  # counts['SRR7218314_1'] +
  # counts['SRR7218315_1']
# counts['SRR7218316_1'] <- counts['SRR7218316_1'] +
  # counts['SRR7218317_1'] +
  # counts['SRR7218318_1']
# counts['SRR7218319_1'] <- counts['SRR7218319_1'] + counts['SRR7218320_1']
# counts['SRR7218324_1'] <- counts['SRR7218324_1'] + counts['SRR7218325_1']
# counts['SRR7218326_1'] <- counts['SRR7218326_1'] + counts['SRR7218327_1']
# counts['SRR7218331_1'] <- counts['SRR7218331_1'] + counts['SRR7218332_1']

# then delete columns, that were averaged
counts <- counts[, !colnames(counts) %in% c('SRR7218298_1', 'SRR7218302_1', 'SRR7218303_1',
          'SRR7218312_1', 'SRR7218313_1', 'SRR7218314_1',
          'SRR7218315_1', 'SRR7218317_1', 'SRR7218318_1',
          'SRR7218320_1', 'SRR7218325_1', 'SRR7218327_1',
          'SRR7218332_1')]

# write table with counts into csv file
write.csv(counts, "/Users/anna/Desktop/BI/NO/results/SRP148970/counts.csv")


#################################################
#### Differential expression analysis with DESeq2
#################################################

# set up the conditions based on the experimental setup
cond_1 <- rep("PD_SN", 5)
cond_2 <- rep("PD_VTA", 4)
cond_3 <- c("control_SN", "control_SN", "control_VTA",
           "control_SN", "control_SN", "control_SN",
           "control_VTA", "control_SN", "control_VTA", "control_SN",
           "control_VTA", "control_SN", "control_VTA", "control_SN",
           "control_VTA", "control_SN", "control_VTA", "control_SN",
           "control_VTA", "control_SN", "control_VTA", "control_VTA",
           "control_SN", "control_VTA", "control_SN", "control_VTA",
           "control_SN", "control_VTA", "control_SN", "control_VTA")

# build the dataframe from the conditions
samples <- names(counts)
condition <- factor(c(cond_1, cond_2, cond_3))
colData <- data.frame(samples = samples, condition = condition)
levels(colData$condition) <- c('control_SN',
                              'control_VTA',
                              'PD_SN', 'PD_VTA')

# divide into two dataframes for SN and VTA neurons
colData_SN <- subset(colData, colData$condition == 'PD_SN' | colData$condition == 'control_SN')
colData_VTA <- subset(colData, colData$condition == 'PD_VTA' | colData$condition == 'control_VTA')
colData_SN$condition <- factor(colData_SN$condition)
levels(colData_SN$condition) <- c('control_SN','PD_SN')
colData_VTA$condition <- factor(colData_VTA$condition)
levels(colData_VTA$condition) <- c('control_VTA','PD_VTA')

# remove rows with almost absent counts
Sums <- rowSums(counts)
table(Sums >= 2)
counts <- counts[Sums >= 2, ]
counts_SN <- counts[colData_SN$samples]
counts_VTA <- counts[colData_VTA$samples]

# create DESeq2 datasets
dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~condition)
dds_SN <- DESeqDataSetFromMatrix(countData = counts_SN, colData = colData_SN, design = ~condition)
dds_VTA <- DESeqDataSetFromMatrix(countData = counts_VTA, colData = colData_VTA, design = ~condition)

# set the reference to be compared
dds_SN$condition <- relevel(dds_SN$condition,"control_SN")
dds_VTA$condition <- relevel(dds_VTA$condition,"control_VTA")

# principal component analysis of samples
for_pca <- rlog(dds)
pca_plot <- plotPCA(object = for_pca, intgroup = 'condition')
pca_plot + theme(text = element_text(size = 16))
# add sample names
pca_plot +
  geom_text(aes(label = name),
            position = position_nudge(x = 3, y = 2),
            size = 2) +
  theme(aspect.ratio = 0.7)

for_pca_SN <- rlog(dds_SN)
for_pca_VTA <- rlog(dds_VTA)
plotPCA(object = for_pca_SN, intgroup = 'condition')
plotPCA(object = for_pca_VTA, intgroup = 'condition')


# run deseq
dds_SN = DESeq(dds_SN)
dds_VTA = DESeq(dds_VTA)


#### Normalized data matrix

# get normalized counts and write this to a file
nc_SN <- counts(dds_SN,normalized = TRUE)
nc_VTA <- counts(dds_VTA,normalized = TRUE)
nc <- counts(dds, normalized = TRUE)

# turn it into a dataframe to have proper column names
dt_SN <- data.frame("id" = rownames(nc_SN), nc_SN)
dt_VTA <- data.frame("id" = rownames(nc_VTA), nc_VTA)
dt <- data.frame("id" = rownames(nc), nc)

# save the normalized data matrix
write.table(dt_SN,
            file = "/Users/anna/Desktop/BI/NO/results/SRP148970/norm-matrix-SN.txt",
            sep = "\t",  row.names = FALSE,
            col.names = TRUE, quote = FALSE)
write.table(dt_VTA,
            file = "/Users/anna/Desktop/BI/NO/results/SRP148970/norm-matrix-VTA.txt",
            sep = "\t",  row.names = FALSE,
            col.names = TRUE, quote = FALSE)
write.table(dt, file = "/Users/anna/Desktop/BI/NO/results/SRP148970/norm-matrix.txt",
            sep = "\t",  row.names = FALSE,
            col.names = TRUE, quote = FALSE)


# format the results
res_SN <- results(dds_SN, alpha = 0.05, lfcThreshold = 1)
res_VTA <- results(dds_VTA, alpha = 0.05, lfcThreshold = 1)
summary(res_SN)
summary(res_VTA)

# logfold change correction and MA plots
resLFC2_SN <- lfcShrink(dds_SN, coef = "condition_PD_SN_vs_control_SN",
                        type="apeglm", res = res_SN, lfcThreshold = 1)
resLFC2_VTA <- lfcShrink(dds_VTA, coef = "condition_PD_VTA_vs_control_VTA",
                         type="apeglm", res = res_VTA, lfcThreshold = 1)
summary(resLFC2_SN)
summary(resLFC2_VTA)

plotMA(resLFC2_SN)
plotMA(resLFC2_VTA)

# volcano plot
# plain ones
EnhancedVolcano(resLFC2_SN, lab = row.names(resLFC2_SN),
                x = 'log2FoldChange', y = 'svalue',
                pCutoff = 0.005, FCcutoff = 1.0,
                title = "",
                subtitle = "", caption = "",
                ylab = expression("-Log"[10]*"S"),
                pointSize = 3,
                labSize = 5,
                axisLabSize = 20,
                legendPosition = 'none')

EnhancedVolcano(resLFC2_VTA, lab = row.names(resLFC2_VTA),
                x = 'log2FoldChange', y = 'svalue',
                pCutoff = 0.005, FCcutoff = 1.0,
                title = "",
                subtitle = "", caption = "",
                ylab = expression("-Log"[10]*"S"),
                pointSize = 3,
                labSize = 5,
                axisLabSize = 20,
                legendPosition = 'none')

# volcano plots with names of NO-associated genes showed
lab_italics <- paste0("italic('", rownames(resLFC2_SN), "')")
selectLab_italics <- paste0("italic('", c('DNM3', 'SNTG1', 'JAK2'), "')")

EnhancedVolcano(resLFC2_SN, lab = lab_italics,
                selectLab = selectLab_italics,
                x = 'log2FoldChange', y = 'svalue',
                pCutoff = 0.005, FCcutoff = 1.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colAlpha = 2/5,
                colConnectors = 'black',
                labCol = 'black',
                boxedLabels = TRUE,
                parseLabels = TRUE,
                title = "Parkinson's disease, substantia nigra",
                subtitle = "", caption = "",
                ylab = expression("-Log"[10]*"S"),
                titleLabSize = 30,
                pointSize = 3,
                labSize = 10,
                axisLabSize = 35,
                legendPosition = 'none')

lab_italics <- paste0("italic('", rownames(resLFC2_VTA), "')")
selectLab_italics <- paste0("italic('", c('ADH5', 'TNFSF12'), "')")

EnhancedVolcano(resLFC2_VTA, lab = lab_italics,
                selectLab = selectLab_italics,
                x = 'log2FoldChange', y = 'svalue',
                pCutoff = 0.005, FCcutoff = 1.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colAlpha = 2/5,
                colConnectors = 'black',
                labCol = 'black',
                boxedLabels = TRUE,
                parseLabels = TRUE,
                title = "Parkinson's disease, ventral tegmental area",
                subtitle = "", caption = "",
                ylab = expression("-Log"[10]*"S"),
                titleLabSize = 30,
                pointSize = 3,
                labSize = 10,
                axisLabSize = 35,
                legendPosition = 'none')

#### Save the results

# sort the results data frame by the padj (or svalue) and LFC columns
sorted_SN = res_SN[with(res_SN, order(padj, -log2FoldChange)), ]
sorted_VTA = res_VTA[with(res_VTA, order(padj, -log2FoldChange)), ]
sorted_SN_lfc2 = resLFC2_SN[with(resLFC2_SN, order(svalue, -log2FoldChange)), ]
sorted_VTA_lfc2 = resLFC2_VTA[with(resLFC2_VTA, order(svalue, -log2FoldChange)), ]

# turn it into a dataframe to have proper column names.
sorted_SN.df = data.frame("id"=rownames(sorted_SN),sorted_SN)
sorted_VTA.df = data.frame("id"=rownames(sorted_VTA),sorted_VTA)
sorted_SN_lfc2.df = data.frame("id"=rownames(sorted_SN_lfc2),sorted_SN_lfc2)
sorted_VTA_lfc2.df = data.frame("id"=rownames(sorted_VTA_lfc2),sorted_VTA_lfc2)

# write the table out
write.table(sorted_SN.df,
            file="/Users/anna/Desktop/BI/NO/results/SRP148970/result_SN.txt",
            sep="\t", col.names=NA, quote=FALSE)
write.table(sorted_VTA.df,
            file="/Users/anna/Desktop/BI/NO/results/SRP148970/result_VTA.txt",
            sep="\t", col.names=NA, quote=FALSE)
write.table(sorted_SN_lfc2.df,
            file="/Users/anna/Desktop/BI/NO/results/SRP148970/result_SN_corrected.txt",
            sep="\t", col.names=NA, quote=FALSE)
write.table(sorted_VTA_lfc2.df,
            file="/Users/anna/Desktop/BI/NO/results/SRP148970/result_VTA_corrected.txt",
            sep="\t", col.names=NA, quote=FALSE)

#### Venn diagrams with up- and downregulated genes
res_SN_noNA <- resLFC2_SN[complete.cases(resLFC2_SN),]
down_SN <- rownames(res_SN_noNA[res_SN_noNA$svalue < 0.005 & res_SN_noNA$log2FoldChange < -1, ])
up_SN <- rownames(res_SN_noNA[res_SN_noNA$svalue < 0.005 & res_SN_noNA$log2FoldChange > 1, ])

venn.diagram(x = list(rownames(resLFC2_SN), down_SN, up_SN),
            category.names = c("All", "Down", "Up"),
            filename = "/Users/anna/Desktop/BI/NO/results/SRP148970/plots/venn_SN_corrected.png",
            imagetype = "png",
            cex = 2, cat.cex = 3,
            fill = c("#B3E2CD", "#FDCDAC", "#96AAD1"),
            fontfamily = "sans",
            cat.fontfamily = "sans")

res_VTA_noNA <- resLFC2_VTA[complete.cases(resLFC2_VTA),]
down_VTA <- rownames(res_VTA_noNA[res_VTA_noNA$svalue < 0.005 & res_VTA_noNA$log2FoldChange < -1, ])
up_VTA <- rownames(res_VTA_noNA[res_VTA_noNA$svalue < 0.005 & res_VTA_noNA$log2FoldChange > 1, ])

venn.diagram(x = list(rownames(resLFC2_VTA), down_VTA, up_VTA),
             category.names = c("All", "Down", "Up"),
             filename = "/Users/anna/Desktop/BI/NO/results/SRP148970/plots/venn_VTA.png",
             imagetype = "png",
             cex = 2, cat.cex = 3,
             fill = c("#B3E2CD", "#FDCDAC", "#96AAD1"),
             fontfamily = "sans",
             cat.fontfamily = "sans")

myCol <- brewer.pal(4, "Pastel2")
venn.diagram(x = list(down_SN, up_SN,
                      down_VTA, up_VTA),
             category.names = c("Down_SN", "Up_SN",
                                "Down_VTA", "Up_VTA"),
             filename = "/Users/anna/Desktop/BI/NO/results/SRP148970/plots/venn_all.png",
             imagetype = "png",
             cex = 1, cat.cex = 1,
             fill = myCol,
             fontfamily = "sans",
             cat.fontfamily = "sans")

# extract lists of genes
write.table(c(up_SN, down_SN),
            file = "/Users/anna/Desktop/BI/NO/results/SRP148970/DEG_SN_corrected.txt",
            row.names = FALSE, quote = FALSE, col.names = FALSE)
write.table(down_SN,
            file = "/Users/anna/Desktop/BI/NO/results/SRP148970/genes_down_SN_corrected.txt",
            row.names = FALSE, quote = FALSE, col.names = FALSE)
write.table(up_SN,
            file = "/Users/anna/Desktop/BI/NO/results/SRP148970/genes_up_SN_corrected.txt",
            row.names = FALSE, quote = FALSE, col.names = FALSE)

write.table(c(up_VTA, down_VTA),
            file = "/Users/anna/Desktop/BI/NO/results/SRP148970/DEG_VTA_corrected.txt",
            row.names = FALSE, quote = FALSE, col.names = FALSE)
write.table(down_VTA,
            file = "/Users/anna/Desktop/BI/NO/results/SRP148970/genes_down_VTA_corrected.txt",
            row.names = FALSE, quote = FALSE, col.names = FALSE)
write.table(up_VTA,
            file = "/Users/anna/Desktop/BI/NO/results/SRP148970/genes_up_VTA_corrected.txt",
            row.names = FALSE, quote = FALSE, col.names = FALSE)
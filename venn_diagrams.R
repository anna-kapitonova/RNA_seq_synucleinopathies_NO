library(VennDiagram)
library(RColorBrewer)

### generate lists with DEG names for all datasets

# SRP213215
genes213215 <- read.table("/Users/anna/Desktop/BI/NO/results/Sasha/genes213215corr_replicates_as_sum.txt")
res_213215_noNA <- genes213215[complete.cases(genes213215),]
down_213215_MSA <- rownames(res_213215_noNA[res_213215_noNA$svalue < 0.005 & res_213215_noNA$log2FoldChange < -1, ])
up_213215_MSA <- rownames(res_213215_noNA[res_213215_noNA$svalue < 0.005 & res_213215_noNA$log2FoldChange > 1, ])
all_213215_MSA <- c(down_213215_MSA, up_213215_MSA)

# SRP058181
all_058181 <- read.table("/Users/anna/Desktop/BI/NO/results/Sasha/gene_names_058181all_corr.txt")[, 1]

# SRP148970
all_148970_SN <- read.table("/Users/anna/Desktop/BI/NO/results/SRP148970/average_replicates/DEG_SN_corrected.txt")[, 1]
all_148970_VTA <- read.table("/Users/anna/Desktop/BI/NO/results/SRP148970/average_replicates/DEG_VTA_corrected.txt")[, 1]

# SRP324001
all_324001_PD <- read.table("/Users/anna/Desktop/BI/NO/results/SRP324001/DEG_PD_corrected.txt")[, 1]
all_324001_PDD <- read.table("/Users/anna/Desktop/BI/NO/results/SRP324001/DEG_PDD_corrected.txt")[, 1]
all_324001_DLB <- read.table("/Users/anna/Desktop/BI/NO/results/SRP324001/DEG_DLB_corrected.txt")[, 1]

### plot diagrams

my_col <- brewer.pal(4, "Pastel2")
venn.diagram(x = list(all_213215_MSA, all_058181,
                      c(all_148970_SN, all_148970_VTA),
                      c(all_324001_PD, all_324001_PDD, all_324001_DLB)),
             category.names = c("SRP213215_putamen_MSA",
                                "SRP058181_BA9_PD",
                                "SRP148970_SN+VTA_PD",
                                "SRP324001_ACC_PD+PDD+DLB"),
             filename = "/Users/anna/Desktop/BI/NO/results/venn_all_datasets.png",
             imagetype = "png",
             cex = 1, cat.cex = 0.5, cat.pos = c(-10, 10, -10, -5),
             fill = my_col,
             fontfamily = "sans",
             cat.fontfamily = "sans")

venn.diagram(x = list(all_058181, all_148970_SN,
                      all_148970_VTA, all_324001_PD),
             category.names = c("Brodmann area 9",
                                "Substantia nigra",
                                "Ventral tegmental area",
                                "Anterior cingulate cortex"),
             filename = "/Users/anna/Desktop/BI/NO/results/venn_all_PD.png",
             imagetype = "png",
             title = "DEG in Parkinson's disease",
             cex = 1, cat.cex = 0.7, cat.pos = c(-30, 1, -140, -200, 1),
             fill = my_col,
             fontfamily = "sans",
             cat.fontfamily = "sans")

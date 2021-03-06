## Analysis of differential expression of genes involved in NO-signaling in synucleinopathies
### 2022 Spring project, Bioinformatics Institute

Students: 

- Anna Kapitonova, Moscow State University
- Alexandra Livanova, Saint-Petersburg State University 

Supervisor:
- Stanislav Bondarev, Saint-Petersburg State University

### Introduction
Synucleinopathies are neurodegenerative diseases that include Parkinson's disease, multiple system atrophy, and dementia with Lewy bodies (Coon & Singer, 2020). During pathogenesis, protein aggregates, in particular, alpha-synuclein, are formed in brain neurons, which leads to loss of control of movements in patients. According to bioinformatic predictions, NOS1AP is also capable of forming protein aggregates in neurons and directly interacts with alpha-synuclein. Based on this, the hypothesis emerged, that NO signaling could be involved in the pathogenesis of synucleinopathies.


### Aim, tasks and data
The **aim** of this project was to evaluate changes in expression level of NO-signaling genes in brain samples from patients with synucleinopathies.

The following **tasks** were set in order to achieve the goal:

1. To assess differential expression of NOS1AP and other genes of NO-signaling in brain tissues of patients with synucleinopathies
2. To compare sets of differentially expressed genes in different brain regions of patients with synucleinopathies

The **available data** at the start of the project were four open RNA-seq datasets of raw reads from different brain tissues of patients with synucleinopathies: 
- SRP058181 - Brodmann area, prefrontal cortex, 42 control, 29 Parkinson's disease (Dumitriu *et al.*, 2016)
- SPR148970 - *substantia nigra* and ventral tegmental area midbrain dopamine neurons, 18 control, 5 Parkinson’s disease (Aguila *et al.*, 2021)
- SRP215213 - putamen, 12 control, 10 multiple system atrophy (data not published)
- SRP324001 - anterior cingulate cortex, 7 control, 7 dementia with Lewy bodies, 7 Parkinson’s disease, 7 Parkinson’s disease with dementia (Feleke *et al.*, 2021)


### Workflow
The workflow of the project is discussed below.

#### Preparing the raw reads
Raw reads in fastq format were downloaded to the server with SRA toolkit. The quality of raw reads was assessed in [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) v0.11.5 (Andrews, 2010) and summary reports were created in [MultiQC](https://multiqc.info/) 1.12 (Ewels *et al.*, 2016):

[1_download_data_and_QC.sh](./1_download_data_and_QC.sh)

#### Alignment
[STAR](https://github.com/alexdobin/STAR) 2.7.10a (Dobin *et al.*, 2013) alignment against GRCh38 human genome was performed with resulting bam files and GeneCounts option:

[2_alignment_to_reference.sh](./2_alignment_to_reference.sh)

#### Variant calling & filtering
[Samtools](https://github.com/samtools/samtools) htslib 1.10.2-3 (Danecek *et al.*, 2021) was used to sort .bam files with alignment, then mpileup was computed for sorted .bam files against GRCh38 human genome, followed by bcftools calling to generate .bcf files. Resulting .bcf files were filtered to obtain SNPs in SNCA, LRRK2, GBA, PRKN genes only:

[3_variant_calling.sh](./3_variant_calling.sh)

#### Constructing table with counts
For each dataset, a table with counts was constructed with R script from ReadsPerGene.out.tab files obtained after STAR alignment. Second columns corresponding to counts for non-stranded libraries were used. Two protocols were provided for technical replicates (if any), where they were summarized or averaged.

#### PCA & DE analysis
PCA using rlog transformation was performed to check the clustarization of groups and quality of replicates. DESeq2 (Love *et al.*, 2014) was used to perform differential expression analysis. Log2FoldChange was corrected with apeglm (Zhu *et al.*, 2019), thresholds for significant DE: s-value < 0.005, |lfc| > 1.

- [SRP058181_count_table&DEA.R](./SRP058181_count_table&DEA.R)
- [SRP215213_count_table&DEA.R](./SRP215213_count_table&DEA.R)
- [SRP148970_count_table&DEA.R](./SRP148970_count_table&DEA.R)
- [SRP324001_count_table&DEA.R](./SRP324001_count_table&DEA.R)

#### Gene enrichment
Lists of genes with significantly changed expression were obtained and uploaded to Gene Ontology, gsea and kobas databases to find out main signaling pathways upregulated and downregulated in synucleinopathies.

#### NOS1AP gene expression
The table with normalized counts was obtained after correction for library size for each of four datasets. Normalized counts for NOS1AP were compared in controls and patients with Parkinson's disease (Mann-Whitney test) in [GraphPad Prism 8 software](http://www.graphpad.com/faq/viewfaq.cfm?faq=1362).

#### Finding genes of NO signaling differentially expressed in synucleinopathies
The list of GO terms with 'nitric oxide' keyword was used to obtain the list of associated genes. Common positions between differentially expressed genes and NO-associated genes were found.

#### Finding intersections in DE in different brain tissues
Venn's diagram (R script) was used to estimate the number of common genes with differential expression between differenent brain tissues. 

[venn_diagrams.R](./venn_diagrams.R)

### Results
1. Expression of NOS1AP does not differ significantly in brain tissues of patients with synucleinopathies, although there is a decreasing tendency in prefrontal cortex (p = 0.08) and substantia nigra (p = 0.09) of patient with Parkinson's disease.

![](./plots/NOS1AP_expression.png)

2. Some NO-signaling genes are differentially expressed in patients with synucleinopathies.

![](./plots/volcano_NO_genes.png)

3. Patterns of DEGs differ tissue- and disease-specifically, only a few common genes were found.

![](./plots/venn_PD.png)

### Literature
1. Aguila, J., Cheng, S., Kee, N., Cao, M., Wang, M., Deng, Q., & Hedlund, E. (2021). Spatial RNA Sequencing Identifies Robust Markers of Vulnerable and Resistant Human Midbrain Dopamine Neurons and Their Expression in Parkinson's Disease. Frontiers in molecular neuroscience, 14, 699562. https://doi.org/10.3389/fnmol.2021.699562
2. Andrews, S. (2010). FastQC: A Quality Control Tool for High Throughput Sequence Data [Online].
3. Coon, E. A., & Singer, W. (2020). Synucleinopathies. Continuum (Minneapolis, Minn.), 26(1), 72–92. https://doi.org/10.1212/CON.0000000000000819
4. Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., Whitwham, A., Keane, T., McCarthy, S. A., Davies, R. M., & Li, H. (2021). Twelve years of SAMtools and BCFtools. GigaScience, 10(2), giab008. https://doi.org/10.1093/gigascience/giab008
5. Dobin, A., Davis, C. A., Schlesinger, F., Drenkow, J., Zaleski, C., Jha, S., Batut, P., Chaisson, M., & Gingeras, T. R. (2013). STAR: ultrafast universal RNA-seq aligner. Bioinformatics (Oxford, England), 29(1), 15–21. https://doi.org/10.1093/bioinformatics/bts635
6. Dumitriu, A., Golji, J., Labadorf, A. T., Gao, B., Beach, T. G., Myers, R. H., Longo, K. A., & Latourelle, J. C. (2016). Integrative analyses of proteomics and RNA transcriptomics implicate mitochondrial processes, protein folding pathways and GWAS loci in Parkinson disease. BMC medical genomics, 9, 5. https://doi.org/10.1186/s12920-016-0164-y
7. Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics (Oxford, England), 32(19), 3047–3048. https://doi.org/10.1093/bioinformatics/btw354
8. Feleke, R., Reynolds, R. H., Smith, A. M., Tilley, B., Taliun, S., Hardy, J., Matthews, P. M., Gentleman, S., Owen, D. R., Johnson, M. R., Srivastava, P. K., & Ryten, M. (2021). Cross-platform transcriptional profiling identifies common and distinct molecular pathologies in Lewy body diseases. Acta neuropathologica, 142(3), 449–474. https://doi.org/10.1007/s00401-021-02343-x
9. Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome biology, 15(12), 550. https://doi.org/10.1186/s13059-014-0550-8
10. R Core Team (2021). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.
11. Zhu, A., Ibrahim, J. G., & Love, M. I. (2019). Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences. Bioinformatics (Oxford, England), 35(12), 2084–2092. https://doi.org/10.1093/bioinformatics/bty895

## Analysis of differential expression of genes involved in NO-signaling in synucleinopathies
Authors: 

- Anna Kapitonova
- Alexandra Livanova 
- Stanislav Bondarev

### Introduction
Synucleinopathies are neurodegenerative diseases that include Parkinson's disease, multiple system atrophy, and dementia with Lewy bodies (Coon & Singer, 2020). 
During pathogenesis, protein aggregates, in particular, alpha-synuclein, are formed in brain neurons, which leads to loss of control of movements in patients. 
According to bioinformatic predictions, NOS1AP is also capable of forming protein aggregates in neurons and directly interacting with alpha-synuclein. 
Based on this, we proposed a hypothesis that NO signaling can probably be involved in the pathogenesis of synucleinopathies.

### Aim, tasks and data
The **aim** of this project was to to evaluate changes in expression level of NO-signaling genes in brain samples from patients with synucleinopathies

The following **tasks** were set in order to achive the goal:

0. To assess differential expression of NOS1AP and other genes of NO-signaling in brain tissues of patients with synucleinopathies
1. To compare sets of differentially expressed genes in different brain tissues of patient in synucleinopathies


The **available data** at the start of the project were four open RNA-seq datasets of raw reads from different brain tissues of patients with synucleinopathies: SRP058181 - Brodmann area, prefrontal cortex, 42 controls, 29 Parkinson's disease
SPR158970 - substantia nigra and ventral tegmental area midbrain dopamine neurons, 43  control, 9 Parkinson’s disease
SRP215213 - putamen, 12 controls, 10 multiple system atrophy
SRP324001 - anterior cingulate cortex, 7 control, 7 dementia with Lewy bodies, 7 Parkinson’s disease, 7 Parkinson’s disease with dementia


### Workflow

The workflow of the project is discussed below.

#### Fetching and dumping of raw reads
.fastq reads were fetched and dumped to the server with SRA toolkit

#### Preparing the raw reads
The quality of raw reads was assessed in [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and summary reports were created in [MultiQC](https://multiqc.info/).

#### Alignment
STAR alignment against GRCh38 human genome was performed with resulting bam files and genecounts option.

#### Variant calling & filtering
mpileup was computed for sorted .bam files against GRCh38 human genome, followed by bcftools calling to generate .bcf files. Resulting .bcf files were filtered to obtain SNPs at SNCA, LRRK2, GBA, PRKN regions only.

#### Constructing table with counts
For each of four datasets table with counts was constructed  with r script from ReadsPerGene.out.tab files obtained after STAR alignment. Second columns corresponding to counts obtained for non-stranded libraries were used. 
Two protocols were provided for technical replicates (if any), where they were summarized or averaged.

#### PCA & DE analysis
DEseq2 was used to perform differential expression analysis. Lof2FoldChange was corrected with apeglm, thresholds for significant DE: s-value<0.005, lfc > 1 or lfc < -1.
PCA using rlog transformation was performed to check the clustarization of groups and quality of replicates. 

#### Gene enrichment
Lists og genes with significantly changed expression were obtained and downloaded to Gene Ontology, gsea, kobas databases to find out main signaling pathways upregulated and downregulated in synucleinopathies.

#### NOS1AP gene expression
Table of normalized counts was obtained after correction for library size for each of four datasets. Normalized counts for NOS1AP were compared in controls and patients with parkinson's disease (Mann-Whitney test) in GraphPrism 8 software (http://www.graphpad.com/faq/viewfaq.cfm?faq=1362).

#### Finding genes of NO signaling differentially expressed in synucleinopathies
The list of GO terms with 'nitric oxide' keyword was obtained. The list of genes corresponding to these GO terms was obtained. Common genes between the one and the list of differentially expressed genes resulted from each of four datasets, were found.

#### Finding intersections in DE in different brain tissues
Venn's diagram (r script) was constructed to reveal common genes with differential expression between differenent brain tissues. 

### Results
0. Some NO-signaling genes are differentially expressed in patients with synucleinopathies
1. Patterns of DEGs differ tissue- and disease-specifically,  only a few common DEGs were found

### Literature
0. Coon, E. A., & Singer, W. (2020). Synucleinopathies. Continuum (Minneapolis, Minn.), 26(1), 72–92. https://doi.org/10.1212/CON.0000000000000819

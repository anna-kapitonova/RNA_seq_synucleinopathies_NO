## Analysis of differential expression of genes involved in NO-signaling in synucleinopathies
Authors: 

- Anna Kapitonova
- Alexandra Livanova 
- Stanislav Bondarev

### Introduction


### Aim, tasks and data
The **aim** of this project was to ... The following **tasks** were set in order to achive the goal:

0. 
1. 

The **available data** at the start of the project were: 4 RNA-seq datasets of raw reads ... 

### Workflow

The workflow of the project presented at the following scheme. Each part of scheme will be discussed below.

#### Preparing the raw reads
The quality of raw paired reads was assessed in [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and summary report were created in [MultiQC](https://multiqc.info/).

#### Alignment
STAR alignment with paired/single reads against reference genome resulted in the mean ~?% coverage of reference by reads.

#### Variant calling & filtering
Resulted .bam files were ...

SNP filtered by means of bcftools with following commands:

...Firstly we filter SNPs with Fisher Strand (FS), Strand Odds Ratio (SOR), Mapping Quality Rank Sum Test (MQRankSum), Read Position Rank Sum Test (ReadPosRankSum), Quality by depth (QD), RMS Mapping Quality (MQ). We run bcftools in mode -e - so it will exclude filtered SNPs...

### Results


### Literature

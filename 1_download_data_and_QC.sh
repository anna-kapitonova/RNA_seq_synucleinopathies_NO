#!/bin/bash

# Download samples (ids according to SRA database are listed in data/accession_list for each dataset)

cat accession_list.txt | while read id
do
/media/array/parkinsons_proj/tools/sratoolkit/sratoolkit.3.0.0-ubuntu64/bin/prefetch $id
echo "for $id prefetch done" 
/media/array/parkinsons_proj/tools/sratoolkit/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump --split-files $id/$id.sra
echo "for $id fastq-dump done"
done


# Quality control with FastQC for all reads in the folder

/media/array/parkinsons_proj/tools/fastqc/FastQC/fastqc *.fastq -o ./fastqc_reports


# Create one common report with MultiQC for sequences from one experiment

cd ./fastqc_reports
multiqc ./ 

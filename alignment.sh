#!/bin/bash

# Download and unzip latest human genome and annotation (GRCh38)

wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gff.gz
gunzip GRCh38_latest_genomic.fna.gz
gunzip GRCh38_latest_genomic.gff.gz


# Convert annotation file from GFF to GTF format




# Index genome with STAR

/media/array/parkinsons_proj/tools/star/STAR-2.7.10a/bin/Linux_x86_64/STAR --runThreadN 16 --runMode genomeGenerate --genomeDir /media/array/parkinsons_proj/data/reference/STAR_index_2 --genomeFastaFiles /media/array/parkinsons_proj/data/reference/GRCh38_latest_genomic.fna --sjdbGTFfile /media/array/parkinsons_proj/data/reference/GRCh38_latest_genomic.gtf  --sjdbOverhang 100


# ALign reads to genome with STAR (for paired reads, dataset SRP215213 example)

cat controls_list.txt | while read id
do
/media/array/parkinsons_proj/tools/star/STAR-2.7.10a/bin/Linux_x86_64/STAR --genomeLoad LoadAndKeep --genomeDir /media/array/parkinsons_proj/data/reference/STAR_index_2 --runThreadN 16 --readFilesIn $id"_1.fastq" $id"_2.fastq" --outFileNamePrefix /media/array/parkinsons_proj/data/Sasha/SRP215213/$id/$id --quantMode GeneCounts --outSAMtype BAM Unsorted
done


# ALign reads to genome with STAR (for single reads, dataset SRP148970 example)

for file in /media/array/parkinsons_proj/data/Anya/SRP148970/*.fastq
do
filename=$(basename -- "$file")
filename="${filename%.*}"
/media/array/parkinsons_proj/tools/star/STAR-2.7.10a/bin/Linux_x86_64/STAR --genomeLoad LoadAndKeep --genomeDir /media/array/parkinsons_proj/data/reference/STAR_index_2 --runThreadN 16 --readFilesIn $file --outFileNamePrefix /media/array/parkinsons_proj/data/Anya/SRP148970/alignments/$filename --quantMode GeneCounts --outSAMtype BAM Unsorted
echo "for $filename alignment done"
done

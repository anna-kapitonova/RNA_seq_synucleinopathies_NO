#!/bin/bash

#Download samples (ids according to SRA database are listed in acq_list)
cat acq_list.txt | while read id
do
/media/array/parkinsons_proj/tools/sratoolkit/sratoolkit.3.0.0-ubuntu64/bin/prefetch $id
echo "for $id prefetch done" 
/media/array/parkinsons_proj/tools/sratoolkit/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump --split-files $id/$id.sra
echo "for $id fastq-dump done"
done

#Quality control with FastQC for all reads in the folder
/media/array/parkinsons_proj/tools/fastqc/FastQC/fastqc *.fastq -o ./fastqc_reports

#Index genome with STAR
/media/array/parkinsons_proj/tools/star/STAR-2.7.10a/bin/Linux_x86_64/STAR --runThreadN 16  --runMode genomeGenerate --genomeDir /media/array/parkinsons_proj/data/reference/STAR_index_2 --genomeFastaFiles /media/array/parkinsons_proj/data/reference/GRCh38_latest_genomic.fna --sjdbGTFfile /media/array/parkinsons_proj/data/reference/GRCh38_latest_genomic.gtf  --sjdbOverhang 100

#ALign reads to genome with STAR 
cat controls_list.txt | while read id
do

/media/array/parkinsons_proj/tools/star/STAR-2.7.10a/bin/Linux_x86_64/STAR --genomeLoad LoadAndKeep --genomeDir /media/array/parkinsons_proj/data/reference/STAR_index_2 --runThreadN 16 --readFilesIn $id"_1.fastq" $id"_2.fastq" --outFileNamePrefix /media/array/parkinsons_proj/data/Sasha/SRP215213/$id/$id --quantMode GeneCounts --outSAMtype BAM Unsorted 

done



#Variant calling
for file in /media/array/parkinsons_proj/data/Sasha/SRP215213/*Aligned.out.bam
do

filename=$(basename -- "$file")
filename="${filename%Aligned.*}"

# sort bam file prior to SNP calling
samtools sort $file -o /media/array/parkinsons_proj/data/Sasha/SRP215213/$filename.out.sorted.bam

# SNP calling
bcftools mpileup -Ou -f /media/array/parkinsons_proj/data/reference/GRCh38_latest_genomic.fna /media/array/parkinsons_proj/data/Sasha/SRP215213/$filename.out.sorted.bam | bcftools call -mv -Ob -o $filename.bcf

# build index for bcf file
bcftools index $filename.bcf

# look at SNCA region
bcftools view -i '%QUAL>=20' $filename.bcf -r NC_000004.12:89724099-89838304 > $filename.txt

# look at LRRK2 region
bcftools view -i '%QUAL>=20' $filename.bcf -r NC_000012.12:40224997-40369285 >> $filename.txt

# look at GBA region
bcftools view -i '%QUAL>=20' $filename.bcf -r NC_000001.11:155234452-155244627 >> $filename.txt

# look at PRKN region
bcftools view -i '%QUAL>=20' $filename.bcf -r NC_000006.12:161347417-162727766 >> $filename.txt

echo "for $file SNP calling done"

done


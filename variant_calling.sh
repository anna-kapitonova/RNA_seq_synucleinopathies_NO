#!/bin/bash

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

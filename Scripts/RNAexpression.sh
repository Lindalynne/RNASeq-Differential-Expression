#!/bin/bash

# Change working directory
cd data


# STEP 1: Run fastqc
mkdir Fastqc
fastqc data/*.fastq -o Fastqc/

# Multi QC
multiqc Fastqc/ --outdir multiqcreports --filename report
echo "Multiqc finished running!"


# Step 2: Run HISAT2

#Index building the reference
hisat2-build Bos_taurus.fa Bostaurus

 
#Creating a list of the sample names
ls *.fastq | cut -d _ -f 1 | sort |uniq > samples.txt

#Mapping Using HISAT2
mkdir mapped
for f in $(<samples.txt);
do
hisat2 -p 8 -x Bostaurus -1 ${f}_1.fastq -2 ${f}_2.fastq -S mapped/${f}.sam
done

echo "HISAT2 finished running!"

#Alignment sorting
for f in $(<samples.txt);
do
samtools sort -@ 8 -o mapped/${f}.bam mapped/${f}.sam
done
 echo "Alignment sorting finished running!"

#removing sam files to create more space.
rm mapped/*.sam


#Get gtf file
#wget https://ftp.ensembl.org/pub/release-108/gtf/bos_taurus/Bos_taurus.ARS-UCD1.2.108.chr.gtf.gz
#gunzip Bos_taurus.ARS-UCD1.2.108.chr.gtf.gz


#Step 3: Trancsript Assembly using Stringtie
for f in $(<samples.txt);
do
stringtie -p 4 -G Bos_taurus.ARS-UCD1.2.108.chr.gtf -o ${f}.gtf -l $f mapped/${f}.bam
done

#Merging all data into a single gtf file
ls *.gtf > mergelist.txt
stringtie --merge -p 4 -G Bos_taurus.ARS-UCD1.2.108.chr.gtf -o stringtie_merged.gtf mergelist.txt
 

echo "Transcript assembly finished!"


# Step 4: Assembly statistics and check merged stats file for precision, senstivity of the transcripts.

gffcompare -r Bos_taurus.ARS-UCD1.2.108.chr.gtf -G -o merged_stats stringtie_merged.gtf

echo "Assembly statistics finished!"

# Step 5: Estimating Abundance

mkdir ballgown
for f in $(<samples.txt)
do
stringtie -e -B -p 4 -G Bos_taurus.ARS-UCD1.2.108.chr.gtf -o ballgown/${f}/${f}.gtf mapped/${f}.bam
done

echo " estimating abundance finished!"



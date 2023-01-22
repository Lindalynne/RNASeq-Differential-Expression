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

for f in $(<samples.txt);
do
#Mapping using HISAT2
hisat2 -p 8 -x Bostaurus -1 ${f}_1.fastq -2 ${f}_2.fastq -S ${f}.sam
done

#sort Bam
#Alignment sorting
for f in $(<samples.txt);
do
samtools sort -@ 8 -o ${f}.bam ${f}.sam
#Index 
samtools index ${f}_sorted.bam
#Mark duplicates
gatk MarkDuplicates --INPUT ${f}_sorted.bam --OUTPUT ${f}_sorted_dedup.bam --METRICS_FILE  ${f}_dup_metrics.txt --REMOVE_DUPLICATES true  --CREATE_INDEX true

done


#wget https://ftp.ensembl.org/pub/release-108/variation/vcf/bos_taurus/bos_taurus.vcf.gz
#gunzip bos_taurus.vcf.gz
tabix -f -p vcf bos_taurus.vcf
gatk CreateSequenceDictionary -R Bos_taurus.fa
samtools faidx Bos_taurus.fa

for f in $(<samples.txt);
do
 #split reads that have N cigar elements.
 gatk SplitNCigarReads -R Bos_taurus.fa -I ${f}_sorted_dedup.bam -O ${f}_sorted_dedupCIGAR.bam
 #perform BaseRecalibration 
 gatk BaseRecalibrator --input ${f}_sorted_dedupCIGAR.bam --output ${f}_recal_data.table --reference Bos_taurus.fa --known-sites bos_taurus.vcf
 #apply the BQSR recalibration to the input BAM file
 gatk ApplyBQSR --bqsr-recal-file ${f}_recal_data.table --input ${f}_sorted_dedupCIGAR.bam --output ${f}_sorted_dedupCIGAR_BQSR_recal.bam --reference Bos_taurus.fa
 
 gatk BaseRecalibrator --input ${f}_sorted_dedupCIGAR_BQSR_recal.bam --output ${f}_post_recal_data.table --reference Bos_taurus.fa --known-sites bos_taurus.vcf
 #call variants using GATK's HaplotypeCaller.
 gatk HaplotypeCaller --reference Bos_taurus.fa --input ${f}_sorted_dedupCIGAR_BQSR_recal.bam --output ${f}.g.vcf.gz --ERC GVCF
 #select SNPs from the input VCF file.
 gatk SelectVariants -R Bos_taurus.fa -V ${f}_raw_variants.vcf --select-type-to-include SNP -O ${f}_raw_snps.vcf
 #select Indels from the input VCF file.
 gatk SelectVariants -R Bos_taurus.fa -V ${f}_raw_variants.vcf --select-type-to-include INDEL -O ${f}_raw_indels.vcf
 #SNP Variant filtering
 gatk VariantFiltration -V ${f}_raw_snps.vcf -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "SOR > 3.0" --filter-name "SOR3" -filter"FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum <-12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum <-8.0" --filter-name "ReadPosRankSum-8" -O ${f}_snps_filtered.vcf
 #Indel variant filtering.
 gatk VariantFiltration -V ${f}_raw_indels.vcf -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "FS > 200.0" --filter-name "FS200" -filter "ReadPosRankSum <-20.0" --filter-name "ReadPosRankSum-20" -O ${f}_indels_filtered.vcf
done

echo "Variant calling finished!"



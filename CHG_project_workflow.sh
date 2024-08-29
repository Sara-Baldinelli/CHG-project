#!/bin/bash

cd Documents/HumanGenomics/

## Sorting and indexing
samtools sort Tumor.bam > Tumor.sorted.bam
samtools index Tumor.sorted.bam
samtools sort Control.bam > Control.sorted.bam
samtools index Control.sorted.bam

## Count reads
samtools view -c Tumor.sorted.bam #### 15039503
samtools view -c Control.sorted.bam #### 19720171
## Reads mapping to reverse strand (f)
samtools view -c -f 16 Tumor.sorted.bam #### 7518303
samtools view -c -f 16 Control.sorted.bam #### 9855853
## Reads mapping to forward strand (F)
samtools view -c -F 16 Tumor.sorted.bam #### 7521200
samtools view -c -F 16 Control.sorted.bam #### 9864318 
## General statistics 
samtools stats Tumor.sorted.bam > Stats_Tumor.txt ##number of properly paired reads 14979936
samtools stats Control.sorted.bam > Stats_Control.txt ##number of properly paired reads 19576046


## REALIGNMENT
java -jar ../Tools/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ../Annotations/human_g1k_v37.fasta -I Tumor.sorted.bam -o realigner.tumor.intervals -L Captured_Regions.bed

java -jar ../Tools/GenomeAnalysisTK.jar -T IndelRealigner -R ../Annotations/human_g1k_v37.fasta -I Tumor.sorted.bam -targetIntervals realigner.tumor.intervals -o Tumor.sorted.realigned.bam -L Captured_Regions.bed

samtools view Tumor.sorted.realigned.bam | grep OC | wc -l > count_aligned_reads_tumor.txt

java -jar ../Tools/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ../Annotations/human_g1k_v37.fasta -I Control.sorted.bam -o realigner.control.intervals -L Captured_Regions.bed

java -jar ../Tools/GenomeAnalysisTK.jar -T IndelRealigner -R ../Annotations/human_g1k_v37.fasta -I Control.sorted.bam -targetIntervals realigner.control.intervals -o Control.sorted.realigned.bam -L Captured_Regions.bed

samtools view Control.sorted.realigned.bam | grep OC | wc -l > count_aligned_reads_control.txt
##count the reads with specific keyword OC flag, which is added after realignment


## RECALIBRATION
## Whole-exome seq: tell GATK to process only the regions supposedly captured by the experiment, excluding off target reads --> -L CancerGenesSel.bed 
java -jar ../Tools/GenomeAnalysisTK.jar -T BaseRecalibrator -R ../Annotations/human_g1k_v37.fasta -I Tumor.sorted.realigned.bam -knownSites ../Annotations/hapmap_3.3.b37.vcf -o recal.tumor.table -L Captured_Regions.bed

java -jar ../Tools/GenomeAnalysisTK.jar -T PrintReads -R ../Annotations/human_g1k_v37.fasta -I Tumor.sorted.realigned.bam -BQSR recal.tumor.table -o Tumor.sorted.realigned.recalibrated.bam -L Captured_Regions.bed --emit_original_quals

java -jar ../Tools/GenomeAnalysisTK.jar -T BaseRecalibrator -R ../Annotations/human_g1k_v37.fasta -I Tumor.sorted.realigned.bam -knownSites ../Annotations/hapmap_3.3.b37.vcf -BQSR recal.tumor.table -o after_recal.tumor.table -L Captured_Regions.bed

java -jar ../Tools/GenomeAnalysisTK.jar -T AnalyzeCovariates -R ../Annotations/human_g1k_v37.fasta -before recal.tumor.table -after after_recal.tumor.table -csv recal_tumor.csv -plots recal_tumor.pdf


java -jar ../Tools/GenomeAnalysisTK.jar -T BaseRecalibrator -R ../Annotations/human_g1k_v37.fasta -I Control.sorted.realigned.bam -knownSites ../Annotations/hapmap_3.3.b37.vcf -o recal.control.table -L Captured_Regions.bed

java -jar ../Tools/GenomeAnalysisTK.jar -T PrintReads -R ../Annotations/human_g1k_v37.fasta -I Control.sorted.realigned.bam -BQSR recal.control.table -o Control.sorted.realigned.recalibrated.bam -L Captured_Regions.bed --emit_original_quals

java -jar ../Tools/GenomeAnalysisTK.jar -T BaseRecalibrator -R ../Annotations/human_g1k_v37.fasta -I Control.sorted.realigned.bam -knownSites ../Annotations/hapmap_3.3.b37.vcf -BQSR recal.control.table -o after_recal.control.table -L Captured_Regions.bed

java -jar ../Tools/GenomeAnalysisTK.jar -T AnalyzeCovariates -R ../Annotations/human_g1k_v37.fasta -before recal.control.table -after after_recal.control.table -csv recal_control.csv -plots recal_control.pdf


## VARIANT CALLING
java -jar ../Tools/GenomeAnalysisTK.jar -T UnifiedGenotyper -R ../Annotations/human_g1k_v37.fasta -I Tumor.sorted.realigned.bam -o Tumor.vcf 
java -jar ../Tools/GenomeAnalysisTK.jar -T UnifiedGenotyper -R ../Annotations/human_g1k_v37.fasta -I Control.sorted.realigned.bam -o Control.vcf 

vcftools --vcf Tumor.vcf --out Tumor.hetero --maf 0.2 --max-maf 0.8 --recode --recode-INFO-all
vcftools --vcf Control.vcf --out Control.hetero --maf 0.2 --max-maf 0.8 --recode --recode-INFO-all


## ANCESTRY
echo Control.sorted.realigned.recalibrated.bam > BAMs_List.txt
Rscript RunEthSEQ.R 

## Somatic copy number
samtools mpileup -q 1 -f ../Annotations/human_g1k_v37.fasta Control.sorted.realigned.recalibrated.bam Tumor.sorted.realigned.recalibrated.bam | java -jar ../Tools/VarScan.v2.3.9.jar copynumber --output-file SCNA --mpileup 1

java -jar ../Tools/VarScan.v2.3.9.jar copyCaller SCNA.copynumber --output-file SCNA.copynumber.called

Rscript CBS.R

## Somatic point mutations

samtools mpileup -q 1 -f ../Annotations/human_g1k_v37.fasta Control.sorted.realigned.recalibrated.bam > Control.sorted.realigned.recalibrated.pileup

samtools mpileup -q 1 -f ../Annotations/human_g1k_v37.fasta Tumor.sorted.realigned.recalibrated.bam > Tumor.sorted.realigned.recalibrated.pileup

java -jar ../Tools/VarScan.v2.3.9.jar somatic Control.sorted.realigned.recalibrated.pileup Tumor.sorted.realigned.recalibrated.pileup --output-snp somatic.pm --output-indel somatic.indel



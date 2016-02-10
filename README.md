# Finding whole genome (WG) SNPs and INDELs with NGS reads
Contributed by Chun-Yuan Huang, 2/10/2016

## Aims:
This workflow is to to find WG SNPs and short INDELs of genomic samples as compared to reference genome using NGS reads. It takes NGS reads (trimmed and quality filtered with TRegGA) and reference genome sequence as inputs, then generates VCF (Variant Call Format) files for SNPs and INDELs of the sample. The original workflow is described in the 3,000 rice genome project (3kRGP) [1] and the SNP-Seek paper [2] in order to discover allelic variants of 3kRGP. The original step-by-step workflow, published online in [3,4], is modified here to accomodate to IU large memory computer cluster Mason. The results can be a valuable extension to our TRegGA workflow for targeted rice variants studies.

## Workflow description [1]:
1. The clean reads were mapped to the temperate japonica Nipponbare reference genome IRGSP-1.0 using the BWA software.
2. The alignment results were then merged and indexed as BAM files. 
3. SNP calling was based on alignment using the Genome Analysis Toolkit 2.0-35 (GATK) and Picard package V1.71. 
4. To minimize the number of mismatched bases for SNP and InDel calling, all reads from each accession were further cleaned by: 
  1. deleting the reads that are unmapped to the reference in the alignment result; 
  2. deleting duplicate reads; 
  3. conducting alignment by the IndelRealigner package in GATK; and 
  4. recalibrating realignments using the BaseRecalibrator package in GATK. 
5. SNP and InDel calling for each sample were performed independently using the UnifiedGenotyper package in GATK with a minimum phred-scaled confidence threshold of 50, and a minimum phred-scaled confidence threshold for emitting variants at 10. 
6. To ensure the quality of variant calling, the conditions for every site in a genome were set at >20 for mapping quality, >50 for variant quality and >2 for the number of supporting reads for every base. 

## Execution description:
1. Edit and setup the parameters as described in 0SOURCE, then `source 0SOURCE`
2. Prepare for the prerequisite files and softwares as described in PREREQ.md
3. Generate qsub script to be run on Mason:
```
PPN=16
VMEM=40
WALLTIME=40
EMAIL=youremail@indiana.edu
JAVA_XMX=8g

##----------------
for i in ${SAMPLE} 
do
echo "
#!/bin/bash
#PBS -m abe
#PBS -l nodes=1:ppn=${PPN},vmem=${VMEM}gb,walltime=${WALLTIME}:00:00
#PBS -M ${EMAIL}
#PBS -N gatk_${i}-on-${REF_SEQNAME}
#PBS -j oe

module add bowtie2/2.2.3
module add bwa/0.7.6a
module add python
module add java
module add samtools/0.1.19
module add picard/1.52
module add gatk/3.4-0

cd ${WGvarSNP_DIR}

#1. Align the paired reads to reference genome using bwa mem.
bwa mem -M -t 8 \
${REF_SEQ} ${i}_1.fq ${i}_2.fq > ${i}.sam

#2. Sort SAM file and output as BAM file
java -Xmx${JAVA_XMX} -jar ${PICARD_DIR}/SortSam.jar \
INPUT=${i}.sam \
OUTPUT=${i}.sorted.bam \
SORT_ORDER=coordinate \
VALIDATION_STRINGENCY=LENIENT \
CREATE_INDEX=TRUE

#3. Fix mate using samtools fixmate, followed with sorting and indexing
# Fix mate with Picard/FixMateInformation.jar is causing errors for downstream steps, so use fixmate from samtools instead.
${SAMTOOLS_DIR}/samtools fixmate ${i}.sorted.bam ${i}.fxmt.bam
java -Xmx${JAVA_XMX} -jar ${PICARD_DIR}/SortSam.jar \
INPUT=${i}.fxmt.bam \
OUTPUT=${i}.fxmt.sorted.bam \
SORT_ORDER=coordinate \
VALIDATION_STRINGENCY=LENIENT \
CREATE_INDEX=TRUE

#4. Mark duplicate reads
java -Xmx${JAVA_XMX} -jar ${PICARD_DIR}/MarkDuplicates.jar \
INPUT=${i}.fxmt.sorted.bam \
OUTPUT=${i}.mkdup.bam \
METRICS_FILE=${i}.metrics \
VALIDATION_STRINGENCY=LENIENT \
CREATE_INDEX=TRUE \
MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000

#5. Add or replace read groups
java -Xmx${JAVA_XMX} -jar ${PICARD_DIR}/AddOrReplaceReadGroups.jar \
INPUT=${i}.mkdup.bam \
OUTPUT=${i}.addrep.bam \
RGID=${i} \
RGPU=${i} \
RGLB=${i} \
RGPL=Illumina \
RGSM=${i} \
RGCN=BGI \
SORT_ORDER=coordinate \
VALIDATION_STRINGENCY=LENIENT \
CREATE_INDEX=TRUE

#6. Realign Target using GenomeAnalysisTK (GATK)
java -Xmx${JAVA_XMX} -jar ${GATK_DIR}/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-I ${i}.addrep.bam \
-R ${REF_SEQ} \
-o ${i}.intervals \
-fixMisencodedQuals \
-nt 8

#7. Indel Realigner
java -Xmx${JAVA_XMX} -jar ${GATK_DIR}/GenomeAnalysisTK.jar \
-T IndelRealigner \
-I ${i}.addrep.bam \
-R ${REF_SEQ} \
-targetIntervals ${i}.intervals \
-fixMisencodedQuals \
-o ${i}.realn.bam

#8. Merge individual BAM files if there are multiple read pairs per sample
# ${SAMTOOLS_DIR}/samtools merge ${i}.merged.bam ${i}.*.realn.bam
# \cp ${i}.merged.bam ${i}.realn.bam

#9. Call variants using Unified Genotyper
java -Xmx${JAVA_XMX} -jar ${GATK_DIR}/GenomeAnalysisTK.jar \
-T UnifiedGenotyper \
-R ${REF_SEQ} \
-I ${i}.realn.bam \
-o ${i}.vcf \
-glm BOTH \
--genotyping_mode DISCOVERY \
-out_mode EMIT_ALL_SITES \
--sample_ploidy 2 \
--min_base_quality_score 20 \
--standard_min_confidence_threshold_for_emitting 10 \
--standard_min_confidence_threshold_for_calling 50 \
--min_indel_count_for_genotyping 2

#10. Select Variants by SNPS and Indels using SelectVariants
java -Xmx${JAVA_XMX} -jar ${GATK_DIR}/GenomeAnalysisTK.jar \
-T SelectVariants \
-R ${REF_SEQ} \
--variant ${i}.vcf \
--excludeFiltered \
--excludeNonVariants \
-o ${i}_env.vcf

java -Xmx${JAVA_XMX} -jar ${GATK_DIR}/GenomeAnalysisTK.jar \
-T SelectVariants \
-R ${REF_SEQ} \
--selectTypeToInclude SNP \
--variant ${i}_env.vcf \
-o ${i}_env_snp.vcf

java -Xmx${JAVA_XMX} -jar ${GATK_DIR}/GenomeAnalysisTK.jar \
-T SelectVariants \
-R ${REF_SEQ} \
--selectTypeToInclude INDEL \
--minIndelSize 5 \
--variant ${i}_env.vcf \
-o ${i}_env_indel.vcf

" > gatk_${i}-on-${REF_SEQNAME}.qsub

done
```
4. Submit jobs on Mason
```
for i in ${SAMPLE}
do
qsub gatk_${i}-on-${REF_SEQNAME}.qsub
done
```

## Reference:
1. The 3,000 rice genomes project. Gigascience. 2014 May 28;3:7.
2. Alexandrov N, et al. SNP-Seek database of SNPs derived from 3000 rice genomes. Nucleic Acids Res. 2015 Jan;43(Database issue):D1023-7.
3. http://aws.amazon.com/public-data-sets/3000-rice-genome/
4. http://s3.amazonaws.com/3kricegenome/README-snp_pipeline.txt


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
2. Prepare for the prerequisite files and softwares by `sh PREREQ.sh`
3. Generate qsub script to be run on Mason: `sh x1-WGvarSNP`
4. Submit jobs on Mason
```
for i in ${SAMPLE}
do
qsub gatk_${i}-on-${REF_SEQNAME}.qsub
done
```
5. cleanup with `sh xcleanup`

## Reference:
1. The 3,000 rice genomes project. Gigascience. 2014 May 28;3:7.
2. Alexandrov N, et al. SNP-Seek database of SNPs derived from 3000 rice genomes. Nucleic Acids Res. 2015 Jan;43(Database issue):D1023-7.
3. http://aws.amazon.com/public-data-sets/3000-rice-genome/
4. http://s3.amazonaws.com/3kricegenome/README-snp_pipeline.txt


# Finding whole genome SNPs and short INDELs with NGS reads
Contributed by Chun-Yuan Huang, 2/15/2016

## Aims:
This workflow is to to find whole genome (WG) SNPs and short INDELs of genomic samples as compared to reference genome using NGS reads. It takes NGS reads (trimmed and quality filtered with TRegGA) and reference genome sequence as inputs, then generates VCF (Variant Call Format) files for SNPs and INDELs of the sample. The original workflow is described in the 3,000 rice genome project (3kRGP) [1] and the SNP-Seek paper [2] in order to discover allelic variants of 3kRGP. The original step-by-step workflow, published online in [3,4], is modified here to accomodate to IU large memory computer cluster Mason. The results can be a valuable extension to our TRegGA workflow for targeted rice variants studies.

## Workflow description [1]:
1. The TRegGA-cleaned reads are mapped to the temperate japonica Nipponbare reference genome IRGSP-1.0 using the BWA software.
2. All reads are further cleaned by the following to minimize the number of mismatched bases for variant calling: 
  * remove unmapped and duplicate reads in the bwa alignment result
  * realign the raw gapped alignment to reduce the number of miscalls of INDELs.
3. Variant (SNP and InDel) calling is based on alignment using the UnifiedGenotyper package of Genome Analysis Toolkit (GATK), with the following constrains (phred-scaled scores):
  * minimum confidence threshold for calling variant at 50, 
  * minimum confidence threshold for emitting variants at 10, 
  * minimam base quality at 20,
  * minimum supporting read counts at 2.

## Workflow execution:
1. Edit and setup the parameters as described in 0SOURCE, then `source 0SOURCE`
2. Edit and prepare for the prerequisite files and softwares as described in PREREQ.sh, then `sh PREREQ.sh`
3. Generate qsub script to be run on Mason: `sh x1-WGvarSNP`
4. Submit jobs on Mason: `sh x2-qsub`. 
5. Cleanup files with `sh xcleanup`
6. Find final outputs in *data/*.

## Sub-directories for workflow implementation:
1. *prereq/*: prerequisite inputs such as retrieval and storage of TRegGA processed reads; retrieval and storage of reference genomes, preparation of BLAST+ database for reference genome.
2. *doc/*: reference and tutorial documents.
3. *bin/*: ancillary codes and scripts.
4. *run/*: main scripts and execution results.
5. *data/*: final outputs and reports.

## Notes: 
1. The workflow default to run a test case using 10% reads from rice cultivar Zhengshan97 against reference rice Japponica Chr10. 
2. PREREQ.sh submits a job to retrieve and index the reference genome. Make sure the job is done before proceeding to the next step.
3. If x2-qsub encounter an error message of "SAM/BAM/CRAM file xxxxx appears to be using the wrong encoding for quality scores", you need to turn on option "-fixMisencodedQuals" in step #6 Realign Target and step #7 Indel Realigner of the qsub script located in run/. Turning on this option universally, however, might incur the risk of having another error "Bad input: while fixing mis-encoded base qualities we encountered a read that was correctly encoded".

## Reference:
1. The 3,000 rice genomes project. Gigascience. 2014 May 28;3:7.
2. Alexandrov N, et al. SNP-Seek database of SNPs derived from 3000 rice genomes. Nucleic Acids Res. 2015 Jan;43(Database issue):D1023-7.
3. http://aws.amazon.com/public-data-sets/3000-rice-genome/
4. http://s3.amazonaws.com/3kricegenome/README-snp_pipeline.txt


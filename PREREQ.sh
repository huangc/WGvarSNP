#!/bin/bash
source ./0SOURCE

# Setup sub-directory for workflow
cd ${WORK_DIR}
mkdir -p ${WORK_DIR}/prereq
mkdir -p ${WORK_DIR}/doc
mkdir -p ${WORK_DIR}/bin
# mkdir -p ${WORK_DIR}/src
mkdir -p ${WORK_DIR}/run
mkdir -p ${WORK_DIR}/data
mkdir -p ${WORK_DIR}/scratch

# Preparation for the reference sequence
cd ${prereq_DIR}
# ln -s ${REF_DIR}/${REFSEQ} .
sh ${bin_DIR}/xgetseq

# Preparation for the sample reads with TRegGA
# Note that reads have been trimmed and quality filtered with TRegGA previously, and are just linked here.
for i in ${SAMPLE}
do
ln -s ${READ_DIR}/${i}/${i}_1.fq .
ln -s ${READ_DIR}/${i}/${i}_2.fq .
done

# Index the reference genome
echo "
#!/bin/bash
#PBS -m abe
#PBS -l nodes=1:ppn=8,vmem=20gb,walltime=00:20:00
#PBS -N prereq-on-${REFSEQNAME}
#PBS -j oe

module add bwa/0.7.6a
module add java
module add samtools/0.1.19
module add picard/1.52

cd ${prereq_DIR}
# Index the reference genome for bwa
bwa index ${REFSEQ}

# Index the reference genome for Picard and GATK
${SAMTOOLS_DIR}/samtools faidx ${REFSEQ}
java -Xmx8g -jar ${PICARD_DIR}/CreateSequenceDictionary.jar REFERENCE=${REFSEQ} OUTPUT=${REFSEQNAME}.dict
" > prereq-on-${REFSEQNAME}.qsub
qsub prereq-on-${REFSEQNAME}.qsub


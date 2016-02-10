### Setup sub-directory for workflow
```
source 0SOURCE
cd ${WORK_DIR}
mkdir -p ${WORK_DIR}/prereq
mkdir -p ${WORK_DIR}/doc
mkdir -p ${WORK_DIR}/bin
mkdir -p ${WORK_DIR}/src
mkdir -p ${WORK_DIR}/run
mkdir -p ${WORK_DIR}/data
mkdir -p ${WORK_DIR}/scratch
```
### Preparation for the reference sequence
```
ln -s ${REF_DIR}/${REFSEQ} .
```
### Index the reference genome for bwa
```
bwa index ${REF_SEQ}
```
### Index the reference genome for Picard and GATK
```
${SAMTOOLS_DIR}/samtools faidx ${REF_SEQ}
java -Xmx8g -jar ${PICARD_DIR}/CreateSequenceDictionary.jar REFERENCE=${REF_SEQ} OUTPUT=${REF_SEQNAME}.dict
```
### Preparation for the sample reads with TRegGA
Note that reads have been trimmed and quality filtered with TRegGA previously, and are just linked here.
```
for i in ${SAMPLE}
do
ln -s ${READ_DIR}/${i}/${i}_1.fq .
ln -s ${READ_DIR}/${i}/${i}_2.fq .
done
```

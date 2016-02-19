## *src/* directory contains prerequisite softwares:
- bwa
- samtools
- picard
- gatk

WGvarSNP is implemented as a collection of shell scripts and ancillary Python codes, so no compilation is required. However, the workflow depends on several third-party programs, and many of which do require compiling and/or additional configuration for your particular system. Please see the cited URLs below for details on the software installation. *src/* is assumed for the installation path, but should be replaced with the actual path.

In the case on the IU Mason cluster, the prerequisite softwares can be loaded from the system:
- module add bwa/0.7.6a
- module add samtools/0.1.19
- module add picard/1.52
- module add gatk/3.4-0

### BWA
See http://bio-bwa.sourceforge.net.
Last update: December 3, 2015.
```bash
cd ${src_DIR}
mkdir bwa cd bwa
wget https://github.com/lh3/bwa/archive/0.7.12.tar.gz
tar -xzf 0.7.12.tar.gz
cd bwa-0.7.12/
make
```

### Samtools
See http://www.htslib.org.
Last update: December 3, 2015.
```bash
cd ${src_DIR}
mkdir samtools
cd samtools
wget https://github.com/samtools/samtools/releases/download/1.2/samtools-1.2.tar.bz2
tar -xjf samtools-1.2.tar.bz2
cd samtools-1.2/
make
```

### Picard
See https://github.com/broadinstitute/picard for the source code and installation instruction.
```bash
cd ${src_DIR}
mkdir picard
cd picard
wget https://github.com/broadinstitute/picard/releases/download/2.1.0/picard-tools-2.1.0.zip
unzip picard-tools-2.1.0.zip

```

### GATK
* See https://github.com/broadinstitute/gatk for the source code and installation instruction.
* See https://www.broadinstitute.org/gatk/download/index for the pre-compiled GATK download instruction
```bash
cd ${src_DIR}
mkdir gatk
cd gatk
# Register, then download the pre-compiled GenomeAnalysisTK-3.5.tar.bz2
tar -jxf GenomeAnalysisTK-3.5.tar.bz2
chmod u+x GenomeAnalysisTK.jar
```

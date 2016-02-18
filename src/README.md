## *src/* directory contains prerequisite softwares:
- bwa/0.7.6a
- samtools/0.1.19
- picard/1.52
- gatk/3.4-0

WGvarSNP is implemented as a collection of shell scripts and ancillary Python codes, so no compilatin is required. However, the workflow depends on several third-party programs, and many of which do require compiling and/or additional configuration for your particular system. Source of the softwares are listed. Please see the cited URLs for details on the software and installation. *src/* is assumed below, but should be replaced with the actual installation path.

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
wget https://github.com/lh3/bwa/archive/0.7.12.tar.gz
tar -xzf 0.7.12.tar.gz
cd bwa-0.7.12/
make
cp bwa ${bin_DIR}
```

### Samtools
See http://www.htslib.org.
Last update: December 3, 2015.
```bash
cd ${src_DIR}
wget https://github.com/samtools/samtools/releases/download/1.2/samtools-1.2.tar.bz2
tar -xjf samtools-1.2.tar.bz2
cd samtools-1.2/
make prefix=${bin_DIR}
make prefix=${bin_DIR} install
```



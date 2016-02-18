## *src/* directory contains prerequisite softwares:
- bwa/0.7.6a
- samtools/0.1.19
- picard/1.52
- gatk/3.4-0

Source of the softwares are listed. Please see the cited URLs for details on the software and installation.

BWA

See http://bio-bwa.sourceforge.net. Last update: December 3, 2015.

cd $TRegGA_DIR/local/src/
wget https://github.com/lh3/bwa/archive/0.7.12.tar.gz
tar -xzf 0.7.12.tar.gz
cd bwa-0.7.12/
make
cp bwa $TRegGA_DIR/local/bin/

Samtools

See http://www.htslib.org. Last update: December 3, 2015.

cd $TRegGA_DIR/local/src/
wget https://github.com/samtools/samtools/releases/download/1.2/samtools-1.2.tar.bz2
tar -xjf samtools-1.2.tar.bz2
cd samtools-1.2/
make prefix=$TRegGA_DIR/local
make prefix=$TRegGA_DIR/local install



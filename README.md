# Reference
https://doi.org/10.1038/s41467-020-15571-8

# RIP-seq pipeline

## 1. Fastq

### 1-1. Fastqc

### 1-2. Add UMI to qnames of fastq 1/2
`umiToQname2 Read1.fastq.gz UMI.fastq.gz Read2.fastq.gz SAMPLE`
umiToQname2 generates two files named SAMPLE.R1.fastq.gz and SAMPLE.R2.fastq.gz.

## 2. Alignment

### 2-1. STAR

(1) Run STAR

`STAR --genomeDir STAR_INDEX --runThreadN N --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --outSAMattributes NH HI AS nM NM MD --outFileNamePrefix SAMPLE SAMPLE.R1.fastq.gz SAMPLE.R2.fastq.gz`

(2) BAM index
`samtools index SAMPLE.Aligned.sortedByCoord.out.bam`

### 2-2. Filter

2-3. Remove PCR duplicates

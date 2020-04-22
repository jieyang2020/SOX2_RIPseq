# Reference
https://doi.org/10.1038/s41467-020-15571-8

# RIP-seq pipeline

## 1. Fastq

### 1-1. Fastqc
Check quality by running fastqc.

### 1-2. Add UMI to read 1 and 2
`umiToQname2 Read1.fastq.gz UMI.fastq.gz Read2.fastq.gz SAMPLE`\
This generates two files named SAMPLE.R1.fastq.gz and SAMPLE.R2.fastq.gz.

## 2. Alignment

### 2-1. STAR

(1) Run STAR

`STAR --genomeDir STAR_INDEX --runThreadN N --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --outSAMattributes NH HI AS nM NM MD --outFileNamePrefix SAMPLE SAMPLE.R1.fastq.gz SAMPLE.R2.fastq.gz`

(2) BAM index
`samtools index SAMPLE.Aligned.sortedByCoord.out.bam`

### 2-2. Remove PCR duplicates by using UMI
This procedue removes PCR duplicates among uniquely-mapped reads.

(1) Sort BAM by qname
`samtools sort -@5 -n SAMPLE.Aligned.sortedByCoord.out.bam -o SAMPLE.nameSorted.bam`\
Optionally, -@5 can be omitted if multiple threads (here the number of threads is 5) are not available.

(2) Remove PCR duplicates 
`samtools view -h -q 255 SAMPLE.nameSorted.bam | python ~/Source/Desmond/dedupped.py | samtools view -Sb - | samtools sort -@5 - > SAMPLE.uniq.dedupped.bam`\
Optionally, -@5 can be omitted if multiple threads (here the number of threads is 5) are not available.


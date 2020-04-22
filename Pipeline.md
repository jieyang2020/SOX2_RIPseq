# RIP-seq pipeline

## 0. Prerequsite

### Softwares
[samtools](http://www.htslib.org/)

[fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

[star](https://github.com/alexdobin/STAR)

[featureCounts](http://subread.sourceforge.net/)

### Genome & annotation
[Gencode vM16](https://www.gencodegenes.org/mouse/release_M16.html)\
Select "comprehensive gene annotation as gtf".

## 1. Fastq

### 1-1. Fastqc
Check quality by running fastqc.

### 1-2. Add UMI to read 1 and 2
`umiToQname2 Read1.fastq.gz UMI.fastq.gz Read2.fastq.gz SAMPLE`

This generates two files named SAMPLE.R1.fastq.gz and SAMPLE.R2.fastq.gz.

## 2. Alignment

### 2-1. STAR

(1) Run STAR

`STAR --genomeDir STAR_INDEX --runThreadN N --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --outSAMattributes NH HI AS nM NM MD --outFileNamePrefix SAMPLE SAMPLE.R1.fastq.gz SAMPLE.R2.fastq.gz`

(2) BAM index

`samtools index SAMPLE.Aligned.sortedByCoord.out.bam`

### 2-2. Remove PCR duplicates by using UMI
This step removes PCR duplicates among uniquely-mapped reads.

(1) Sort BAM by qname

`samtools sort -@5 -n SAMPLE.Aligned.sortedByCoord.out.bam -o SAMPLE.nameSorted.bam`\
Optionally, -@5 can be omitted if multiple threads (here the number of threads is 5) are not available.

(2) Remove PCR duplicates 

`samtools view -h -q 255 SAMPLE.nameSorted.bam | python ~/Source/Desmond/dedupped.py | samtools view -Sb - | samtools sort -@5 - > SAMPLE.uniq.dedupped.bam`\
Optionally, -@5 can be omitted if multiple threads (here the number of threads is 5) are not available.

## 3. Count
This step counts sequencing reads in gene body across genomes.

`featureCounts SAMPLE.uniq.dedupped.bam -a ANNOTATION.gtf -p -s 1 -t gene -o SAMPLE.fCounts`


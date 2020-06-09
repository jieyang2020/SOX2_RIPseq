# RIP-seq pipeline

## 0. Prerequsite

### Softwares
[Samtools](http://www.htslib.org/)

[FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

[STAR](https://github.com/alexdobin/STAR)

[FeatureCounts](http://subread.sourceforge.net/)

Additionally, download all the manual programs in the folder "Programs" in this github.

### Genome & annotation
[Gencode vM16](https://www.gencodegenes.org/mouse/release_M16.html)\
Select "comprehensive gene annotation as gtf".

## 1. Fastq

### 1-1. Fastqc
Check quality by running fastqc.

### 1-2. Add UMI information to read 1 file and read 2 file.
Add UMI to header in fastq. The resulting header will have format of [PASS/FAIL]:[UMI]::[orginal qname]. If any nucletodie in UMI has quality lower than 10, I put FAIL, otherwise PASS.

`umiToQname2 Read1.fastq.gz Read2.fastq.gz UMI.fastq.gz SAMPLE`\
Read1.fastq.gz: Read 1 file of paired-end sequencing.\
UMI.fastq.gz: UMI read file.\
Read2.fastq.gz: Read 2 file of paired-end sequencing.\
SAMPLE: output file prefix.

This generates two files named SAMPLE.R1.fastq.gz and SAMPLE.R2.fastq.gz.

## 2. Alignment

### 2-1. STAR

(1) Run STAR

`STAR --genomeDir STAR_INDEX --runThreadN N --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix SAMPLE SAMPLE.R1.fastq.gz SAMPLE.R2.fastq.gz`\
--genomeDir STAR_INDEX: see STAR manual.\
--runThreadN N: see STAR manual.

(2) BAM index

`samtools index SAMPLE.Aligned.sortedByCoord.out.bam`

### 2-2. Remove PCR duplicates by using UMI
This step removes PCR duplicates among uniquely-mapped reads.

(1) Sort BAM by qname

`samtools sort -@5 -n SAMPLE.Aligned.sortedByCoord.out.bam -o SAMPLE.nameSorted.bam`\
"-@5" should be omitted if multiple threads (here the number of threads is 5) are not available.

(2) Remove PCR duplicates 

`samtools view -h -q 255 SAMPLE.nameSorted.bam | python dedupped2.py | samtools view -Sb - | samtools sort -@5 - > SAMPLE.uniq.dedupped.bam`\
"-@5" should be omitted if multiple threads (here the number of threads is 5) are not available.

## 3. Count
This step counts sequencing reads in gene bodies across genome.

`featureCounts SAMPLE.uniq.dedupped.bam -a ANNOTATION.gtf -p -s 1 -t gene -o SAMPLE.fCounts`


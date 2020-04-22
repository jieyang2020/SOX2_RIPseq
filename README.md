# Reference
https://doi.org/10.1038/s41467-020-15571-8

# RIP-seq pipeline

### 1. Fastq

**1-1. Fastqc**

**1-2. Add UMI to qnames of fastq 1/2**
`umiToQname R1.fastq.gz R3.fastq.gz R2.fastq.gz OUTPUT_PREFIX`

1-2. FASTQC

2. Alignment

2-1. STAR

2-2. Filter alignment

2-3. Remove PCR duplicates

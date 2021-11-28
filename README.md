# Bioinformatics
Sequencing was conducted by the Beijing Genomics Institute (BGI) using the DNBseq platform to produce 100bp paired-end reads.

## Quality filtering
We received forward and reverse files for each individual that had been cleaned by BGI to remove: barcode sequences, DNBseq adaptors, low quality reads (50% of quality scores <10), and reads containing >10% unidentified bases. 

We ran basic quality control reports for each file using ```FastQC```

```
FastQC code here
```

A report for all individuals was compiled using ```MultiQC``` (which metrics were checked?)

```
MultiQC code here
```


## Alignment

The *Senecio lautus* reference genome was indexed using ```BWA```.

```
bwa index reference.fasta
```

For each individual, we aligned reads to the reference genome and added read groups using the ```BWA-MEM v0.7.13``` algorithm [(Li and Durbin, 2009)] (https://pubmed.ncbi.nlm.nih.gov/19451168/) under default parameters. The resulting BAM files were sorted using ```Samtools v1.3```.


```
bwa mem \
-M -R "@RG\tSM:ind1\tID:ind1\tLB:ind1\tPL:ILLUMINA\tPU:ind1" \
        reference.fasta \
        ind1_1.fq.gz \
        ind2_2.fq.gz |
samtools sort -@ 12 -T ind1 -o ind1_sorted.bam 
```


## Cleaning BAMs

## Calculating allele frequency


# Bioinformatics
Sequencing was conducted by the Beijing Genomics Institute (BGI) using the DNBseq platform to produce 100bp paired-end reads at approximately ~5X coverage.
Samples were plated randomly with despect to Dune/Headland morphology (for  natural population samples) or gravitropic/agravitropic status (for advanced recombinant
samples). X lanes were used.

##### A note on differences between populations in this pipeline
The bioformatic processing of the datasets (natural populations and advancted recombinant populations) was conducted by different researchers at separate times. Extremely similar pipelines were used overall, with common programs used for all major steps, as is evidenced below. The few minor cleaning/processing steps where the piplines diverge or use different program reflect personal tool preferences rather than a difference in data structure or needs.

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

The *Senecio lautus* reference genome was indexed using ```BWA v0.7.13```.

```
bwa index reference.fasta
```

For each individual, we aligned reads to the reference genome and added read groups using the ```BWA-MEM v0.7.13``` algorithm [(Li and Durbin, 2009)](https://pubmed.ncbi.nlm.nih.gov/19451168/) under default parameters. The resulting BAM files were sorted using ```Samtools v1.3```[(Li et al. 2009)](https://pubmed.ncbi.nlm.nih.gov/19505943/).


```
bwa mem \
-M -R "@RG\tSM:ind1\tID:ind1\tLB:ind1\tPL:ILLUMINA\tPU:ind1" \
        reference.fasta \
        ind1_1.fq.gz \
        ind2_2.fq.gz |
samtools sort -T ind1 -o ind1_sorted.bam 
```


## Cleaning BAMs

For the advanced recombinant population, sorted BAM files were cleaned using ```Picard v2.22``` CleanSam to softclip reads extending beyond the reference genome, and set unmapped quality scores to 0

java -Xmx2g -jar /home/uqralls1/programs/picard.jar CleanSam \
	INPUT=${SAMPLE}_sorted.bam \
	OUTPUT=${SAMPLE}_cln.sorted.bam

### Marking PCR duplicates

For natural population samples, ```samblaster v.0.1.24``` was used on default parameter (using the option to work from a Samtools-sorted file) to mark PCR duplicates for removal.

```
samblaster -M
```

For natural population samples, ```Picard``` was used on default parameter (using the option to work from a sorted file) to mark PCR duplicates for removal.
java -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1 -Xmx4g -jar picard.jar MarkDuplicates \
	INPUT=${SAMPLE}_cln.sorted.bam \
	OUTPUT=${SAMPLE}_mdup.cln.sorted.bam \
	REMOVE_DUPLICATES=false \
	ASSUME_SORTED=true \
	VALIDATION_STRINGENCY=SILENT \
	READ_NAME_REGEX=null \
	MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=900 \
	METRICS_FILE=${SAMPLE}_mdup.cln.sorted.metrics
        

Regions around indels were realigned using GenomeAnalysisToolKit v3.8 (Broad Institute 2018). Using the low-coverage variant caller ANGSD v0.930 (Korneliussen et al. 2014) variable sites were called in regions from the auxin gene set across all study populations, then allele frequency at these sites calculated jointly within each study population.

## Calculating allele frequency


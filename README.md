# Data availability

Sequencing data files generated at each step of this process are available on The University of Queensland Research Data Manager repository ().

# Bioinformatics
Sequencing was conducted by the Beijing Genomics Institute (BGI) using the DNBseq platform to produce 100bp paired-end reads at approximately ~5X coverage.
Samples were plated randomly with respect to Dune/Headland ecotype (for natural population samples) or gravitropic/agravitropic status (for advanced recombinant
samples). X lanes were used.

##### A note on differences in processing of natural and recombinant population data
The bioinformatic processing of these datasets (natural populations and advanced recombinant populations) was conducted by different researchers at separate times. Extremely similar pipelines were used overall, with common programs used for all major steps, as is evidenced below. The few minor file cleaning/processing steps where the pipelines diverge or use different programs reflects only personal tool preferences.

## Quality filtering
We received two FASTQ files for each individual (one containing forward reads, the other containing reverse reads). These FASTQ files  had been cleaned by BGI to remove: barcode sequences, DNBseq adaptors, low quality reads (50% of quality scores <10), and reads containing >10% unidentified bases. 

Basic quality control reports were run for each file using ```FastQC``` [(Andrews, 2010)](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) on default parameters.

```
fastqc ind1_1.fq.gz
```

A report for all individuals was compiled from the ```FastQC``` outputs using ```MultiQC v1.8``` [(Ewels et al. 2016)](https://academic.oup.com/bioinformatics/article/32/19/3047/2196507) to screen for any quality-control anomalies (e.g. per-base sequence quality, per-base N content, per-sequence GC content, sequence duplication levels, overrepresented sequences etc.). More information on these measures can be found in these resources from [Michigan State University Research Technology Support Facility](https://rtsf.natsci.msu.edu/genomics/tech-notes/fastqc-tutorial-and-faq/) and [Harvard Chan Bioinformatics Core](https://hbctraining.github.io/Intro-to-rnaseq-hpc-salmon/lessons/qc_fastqc_assessment.html).

```
multiqc directory/with/fastqc/outputs
```
All files were found to be of sound quality, so no further pre-alignment cleaning was undertaken.

## Alignment

The copy of the *Senecio lautus* reference genome used in this study, produced by [James et al. (2021)](https://academic.oup.com/mbe/article/38/11/4805/6319724#309303427), was indexed using ```BWA v0.7.13``` [(Li and Durbin, 2009)](https://pubmed.ncbi.nlm.nih.gov/19451168/). 

```
bwa index reference.fasta
```

For each individual, reads were aligned to the reference genome and read groups added using the ```BWA-MEM v0.7.13``` algorithm [(Li and Durbin, 2009)](https://pubmed.ncbi.nlm.nih.gov/19451168/) under default parameters. The resulting BAM files were sorted using ```Samtools v1.3``` [(Li et al. 2009)](https://pubmed.ncbi.nlm.nih.gov/19505943/) sort function.


```
bwa mem \
-M -R "@RG\tSM:ind1\tID:ind1\tLB:ind1\tPL:ILLUMINA\tPU:ind1" \
        reference.fasta \
        ind1_1.fq.gz \
        ind2_2.fq.gz |
samtools sort -T ind1 -o ind1_sorted.bam 
```

### Alignment statistics
Alignment statistics were calculated for each individual using ```Samtools v1.3``` flagstat function to check alignment had proceeded successfully (checking total number of reads aligned, total number of reads properly paired etc.). More information on this can be found in the [Samtools flagstat documentation](http://www.htslib.org/doc/samtools-flagstat.html).

```
samtools flagstat ind1_sorted.bam &> ind1_stats.txt
```


## Cleaning BAMs

### Basic cleaning
For the advanced recombinant population, sorted BAM files for each individual were cleaned using ```Picard v2.22``` [(Broad Institute, 2018)](http://broadinstitute.github.io/picard/) CleanSam to softclip reads extending beyond the reference genome, and set unmapped quality scores to 0.

```
java -Xmx2g -jar picard.jar CleanSam \
	INPUT=ind1_sorted.bam \
	OUTPUT=ind1_cln.sorted.bam
```

### Marking PCR duplicates

For natural population samples, ```samblaster v.0.1.24``` [(Faust and Hall, 2014)](https://academic.oup.com/bioinformatics/article/30/17/2503/2748175) was used with default parameters (using the option to work from a Samtools-sorted file) to mark PCR duplicates for removal.

```
samblaster -M ind1_sorted.bam > ind1_mdup.cln.sorted.bam
```

For recombinant population samples, ```Picard v2.22``` MarkDuplicates was used on default parameters (using the option to work from a sorted file) to mark PCR duplicates for removal.

```
java -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1 -Xmx4g -jar picard.jar MarkDuplicates \
	INPUT=ind1_cln.sorted.bam \
	OUTPUT=ind1_mdup.cln.sorted.bam \
	REMOVE_DUPLICATES=false \
	ASSUME_SORTED=true \
	VALIDATION_STRINGENCY=SILENT \
	READ_NAME_REGEX=null \
	MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=900 \
	METRICS_FILE=ind1_mdup.cln.sorted.metrics
```

### Indexing BAMs
Cleaned and sorted BAM files with PCR duplicates marked were then indexed using ```Samtools v1.3``` index function.

```
samtools index ind1_mdup.cln.sorted.bam 
```

### Realigning around indels
Regions around indels were realigned using ```GenomeAnalysisToolKit (GATK) v3.8``` [(Van der Auwera and O'Connor, 2020)](https://www.oreilly.com/library/view/genomics-in-the/9781491975183/) CreateSequenceDictionary.

First, a 'dictionary file' of the *S. lautus* reference genome was created using ```Picard v2.22``` CreateSequenceDictionary.

```
java -jar picard.jar CreateSequenceDictionary \
  -REFERENCE reference.fasta
```

Then, the *S. lautus* reference genome was indexed with ```Samtools v1.3``` faidx.

```
samtools faidx reference.fasta
```

Targets for realignment were identified using ```GATK v3.8``` RealignerTargetCreator.

```
java -jar GenomeAnalysisTK.jar \
        -T RealignerTargetCreator \
        -R reference.fasta \
        -I ind1_mdup.cln.sorted.bam \
        -nt 2 \
        -o ind1.intervals
```


Finally, the realignment was performed with ```GATK v3.8``` IndelRealigner, using the targets for realignment, dictionary, and reference files generated in the previous steps.

```
java -jar /home/uqralls1/programs/GATK/GenomeAnalysisTK.jar \
        -T IndelRealigner \
        -R reference.fasta \
        -I ind1_mdup.cln.sorted.bam \
        -targetIntervals ind1.intervals \
        -o ind1_rln.mdup.cln.bam
```

### Final index and validation of BAM files

The final BAM files were indexed  using ```Samtools v1.3``` index function.

```
samtools index ind1_rln.mdup.cln.sorted.bam 
```

Files were validated using ```Picard v2.22``` ValidateSamFile.

```
java -jar /opt/biotools/picard/picard.jar ValidateSamFile \
        INPUT= ind1_rln.mdup.cln.sorted.bam \
        OUTPUT=ind1.out \
        MODE=VERBOSE \
        MAX_OPEN_TEMP_FILES=1000	
```

## Calculating allele frequencies

### Calling variable sites in target gene-regions
Using the low-coverage variant caller ```ANGSD v0.930``` [(Korneliussen et al. 2014)](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-014-0356-4), variable sites were called in regions from the auxin and shoot gravitropism gene set. 

Variable sites were called across all individuals within the population type (i.e. all six natural populations called together; both recombinant populations called together).

```
angsd -bam bam-file-paths.txt \
        -GL 1 \
        -rf gene-regions.txt \
        -out snp-sites \
        -doMaf 2 \
        -SNP_pval 1e-6 \
        -doMajorMinor 1 \
        -minMaf 0.05 \
        -nThreads 10
```

*Parameter notes (more information available in [ANGSD documentation](http://www.popgen.dk/angsd/index.php/ANGSD)):*
* *GL 1: calculates genotype likelihood with the Samtools method*
* *doMaf 2: assumes fixed major allele inferred from genotype likelihoods (GLs), unknown minor (sums GLs of alleles to determine)*
* *SNP_pval 1e-6: keeps only sites with a p-value less than 1e-6*
* *doMajorMinor 1: uses a maximum likelihood approach to choose major and minor alleles*
* *minMaf 0.05: filters for sites with minimum minor allele freq >0.05*



A text file (.txt) version of the variable sites files generated in the previous step (.mafs.gz) was created by extracting contig, snp position, major and minor allele (in columns 1 - 4, respectively)

```
zcat snp-sites.mafs.gz | awk 'NR!=1{print $1"\t"$2"\t"$3"\t"$4}' > snp-sites.txt
```

The .txt file from the previous step was then sorted to ensure contigs and snp positions were in order.

```
sort -k1 snp-sites.txt > snp-sites-sorted.txt
```

The sorted variable sites file was indexed using ```ANGSD v0.930``` sites index.
```
angsd sites index snp-sites-sorted.txt
```

### Joint allele frequency calling
Allele frequency at these sites calculated jointly within each population, again using ```ANGSD v0.930```, at the variable sites found in the auxin and gravitropism gene regions, as identified in the previous step.

```
${ANGSD}/angsd -bam ${BAMS} \
        -GL 1 \
        -doMaf 1 \
        -doMajorMinor 3 \
        -rf regions.txt \
        -sites snp-sites.txt \
        -out outputfile \
        -nThreads 10
```	
	
*Parameter notes (more information available in [ANGSD documentation](http://www.popgen.dk/angsd/index.php/ANGSD)):*
+ *GL 1: calculates genotype likelihood with the Samtools method*
+ *doMaf 1: uses fixed major and minor alleles (specified by -sites argument using the snp-sites-sorted.txt file produced earlier)*
+ *doMajorMinor 3: uses major and minor alleles (provided in the -sites argument as the .txt file produced earlier)*

This gave one allele frequency file [(.mafs.gz)](http://www.popgen.dk/angsd/index.php/Allele_Frequencies#.mafs.gz) per population.

# File processing

Allele frequency files for all populations were combined in R (scripts available in files) to ensure that major and minor alleles were consistent across populations, and only variable sites that were successfully called in all populations were kept. This was done separately for the natural populations and the recombinant populations, such that there was one final combined allele frequency file for the natural populations, and one final combined allele frequency file for the recombinant population.

Auxin and gravitropism gene names and functions from the original gene set were applied to each site in the combined files using a perl script (available in files).

Final filtering and analysis was done using [JMP](https://www.jmp.com/en_au/home.html) statistical software. 

For the natural populations, sites were retained if they were sequenced in a minimum of three individuals per population. For the recombinant population, sites were retained if they were sequenced in a minimum of 5 individuals per tail.
The following files were used for data analysis:

* The allele frequencies for the natural populations (76,716 SNPs in 560 genes), used for calculating the outlier SNPs and genes: [alleleFreqsNaturalPops.xlsx](Data Files/ alleleFreqsNaturalPops.xlsx)

* The allele frequencies for the recombinant population (79,556 SNPs in 568 genes), used for calculating the outlier SNPs and genes: [alleleFreqsF11s.xlsx](Data Files/ alleleFreqsF11s.xlsx)

* The combined allele frequencies for the natural populations and recombinant populations (merging of the above two files), used for plotting the change in allele frequency graphs: [alleleFreqsCombined.xlsx](Data Files/alleleFreqsCombined.xlsx)

* Summary table of the outlier genes across localities and the recombinant population: [outlierGenes.xlsx](Data Files/outlierGenes.xlsx)

* Summary table of the functions of the outlier genes: [functions.xlsx](Data Files/functions.xlsx)






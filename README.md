# Bioinformatics
Sequencing was conducted by the Beijing Genomics Institute (BGI) using the DNBseq platform to produce 100bp paired-end reads at approximately ~5X coverage.
Samples were plated randomly with respect to Dune/Headland morphology (for  natural population samples) or gravitropic/agravitropic status (for advanced recombinant
samples). X lanes were used.

##### A note on differences in processing of natural and recombinant population data
The bioinformatic processing of these datasets (natural populations and advanced recombinant populations) was conducted by different researchers at separate times. Extremely similar pipelines were used overall, with common programs used for all major steps, as is evidenced below. The few minor file cleaning/processing steps where the pipelines diverge or use different programs reflects only personal tool preferences.

## Quality filtering
We received forward and reverse files for each individual that had been cleaned by BGI to remove: barcode sequences, DNBseq adaptors, low quality reads (50% of quality scores <10), and reads containing >10% unidentified bases. 

Basic quality control reports were run for each file using ```FastQC``` [(Andrews, 2010)](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

```
fastqc ind1_1.fq.gz
```

A report for all individuals was compiled using ```MultiQC v1.8``` [(Ewels et al. 2016)](https://academic.oup.com/bioinformatics/article/32/19/3047/2196507) to screen for any quality-control anomalies.

```
multiqc directory/with/fastqc/outputs
```


## Alignment

The *Senecio lautus* reference genome was indexed using ```BWA v0.7.13``` on default parameters [(Li and Durbin, 2009)](https://pubmed.ncbi.nlm.nih.gov/19451168/).

```
bwa index reference.fasta
```

For each individual, reads were aligned to the reference genome and read groups added using the ```BWA-MEM v0.7.13``` algorithm [(Li and Durbin, 2009)](https://pubmed.ncbi.nlm.nih.gov/19451168/) under default parameters. The resulting BAM files were sorted using ```Samtools v1.3```[(Li et al. 2009)](https://pubmed.ncbi.nlm.nih.gov/19505943/) sort function.


```
bwa mem \
-M -R "@RG\tSM:ind1\tID:ind1\tLB:ind1\tPL:ILLUMINA\tPU:ind1" \
        reference.fasta \
        ind1_1.fq.gz \
        ind2_2.fq.gz |
samtools sort -T ind1 -o ind1_sorted.bam 
```

### Alignment statistics
Alignment statistics were calculated for each individual using ```Samtools v1.3```flagstat function on default parameters.

```
samtools flagstat ind1_sorted.bam &> ind1_stats.txt
```

## Cleaning BAMs

### Basic cleaning
For the advanced recombinant population, sorted BAM files were cleaned using ```Picard v2.22``` [(Broad Institute, 2018)](http://broadinstitute.github.io/picard/) CleanSam to softclip reads extending beyond the reference genome, and set unmapped quality scores to 0.

```
java -Xmx2g -jar picard.jar CleanSam \
	INPUT=ind1_sorted.bam \
	OUTPUT=ind1_cln.sorted.bam
```

### Marking PCR duplicates

For natural population samples, ```samblaster v.0.1.24``` [(Faust and Hall, 2014)](https://academic.oup.com/bioinformatics/article/30/17/2503/2748175) was used on default parameters (using the option to work from a Samtools-sorted file) to mark PCR duplicates for removal.

```
samblaster -M ind1_sorted.bam > ind1_mdup.cln.sorted.bam
```

For natural population samples, ```Picard v2.22``` MarkDuplicates was used on default parameters (using the option to work from a sorted file) to mark PCR duplicates for removal.

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
samtools faidx /90days/uqralls1/Reference/Senecio.contigs.fasta
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

The realignment was performed using ```GATK v3.8``` IndelRealigner.

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

And validated using ```Picard v2.22``` ValidateSamFile.

```
java -jar /opt/biotools/picard/picard.jar ValidateSamFile \
        INPUT= ind1_rln.mdup.cln.sorted.bam \
        OUTPUT=ind1.out \
        MODE=VERBOSE \
        MAX_OPEN_TEMP_FILES=1000	
```

## Calculating allele frequency

### Calling variable sites in target gene-regions
Using the low-coverage variant caller ```ANGSD v0.930``` [(Korneliussen et al. 2014)](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-014-0356-4), variable sites were called in regions from the auxin and shoot gravitropism gene set (regions available in files). This was done for all populations (separately for natural and recombinants). 

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


Parameter notes (more information available in [ANGSD documentation](http://www.popgen.dk/angsd/index.php/ANGSD)):
* GL 1: calculates genotype likelihood with the Samtools method
* doMaf 2: assumes fixed major allele inferred from genotype likelihoods (GLs), unknown minor (sums GLs of alleles to determine)
* SNP_pval 1e-6: keeps only sites with a p-value less than 1e-6
* doMajorMinor 1: uses a maximum likelihood approach to choose major and minor alleles
* minMaf 0.05: filters for sites with minimum minor allele freq >0.05

A text file (.txt) version of the variable sites file generated in the previous step (.mafs.gz) was created by extracting contig, snp position, major and minor allele (in columns 1 - 4, respectively)

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
Allele frequency at these sites calculated jointly within each population, again using ```ANGSD v0.930```.

```
${ANGSD}/angsd -bam ${BAMS} \
        -GL 1 \
        -doMaf 1 \
        -doMajorMinor 3 \
        -rf ${REGIONS} \
        -sites ${SITES} \
        -out ${OUTFILE} \
        -nThreads 10
```	
	
*Parameter notes (more information available in [ANGSD documentation](http://www.popgen.dk/angsd/index.php/ANGSD)):
* GL 1: calculates genotype likelihood with the Samtools method
* doMaf 1: uses fixed major and minor alleles (specified by -sites argument using the snp-sites-sorted.txt file produced earlier)
* doMajorMinor 3: uses major and minor alleles (provided in the -sites argument as the .txt file produced earlier)*


# File processing

Allele frequency files for all populations, both natural and recombinant, were combined in R (script available in files). Only variable sites that were successfully called in all populations were kept. 

Auxin and gravitropism gene names and functions from the original gene set were applied to each site in the combined file using a perl script (available in files).



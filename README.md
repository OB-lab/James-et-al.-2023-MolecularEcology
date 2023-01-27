# Data availability

Sequencing data files generated at each step of this process will be available on The University of Queensland Research Data Manager repository at the time of publication.

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

### Calling variant and invariant allele frequencies in the auxin gene-regions
Using the low-coverage variant caller ```ANGSD v0.930``` [(Korneliussen et al. 2014)](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-014-0356-4), sites were called in regions from the auxin and shoot gravitropism gene set. 

Generation of target sites:

```
angsd -bam bam-file-paths.txt \
      -GL 2 \
      -rf gene-regions.txt \
      -out snp-sites \
      -doMaf 2 \
      -doMajorMinor 1 \
      -remove_bads 1 \
      -minMapQ 30 \
      -minQ 20 \
      -nThreads 24
```

*Parameter notes (more information available in [ANGSD documentation](http://www.popgen.dk/angsd/index.php/ANGSD)):*
* *bam bam-file-paths.txt: the paths to all bam files (natural and MAGIC populations)*
* *rf gene-regions.txt: the regions of the auxin genes(contig:start-finish)*
* *GL 2: calculates genotype likelihood with the GATK method*
* *doMaf 2: assumes fixed major allele inferred from genotype likelihoods (GLs), unknown minor (sums GLs of alleles to determine)*
* *doMajorMinor 1: uses a maximum likelihood approach to choose major and minor alleles*
* *remove_bads 1: remove bad reads*
* *minMapQ 30: minimum mapping quality of 30*
* *minQ 20: minimum per-base quality of 20*

This generates snp-sites.txt, which contains the major and minor alleles across all populations for all sites in the auxin gene-regions. Allele frequencies were then called for each population separately (i.e., each of the six natural populations and each tail of the MAGIC population).

Calling allele frequencies (done separately per population):

```
angsd -bam [pop]-bam-paths.txt \
      -GL 2 \
      -doMaf 1 \
      -doMajorMinor 3 \
      -rf gene-regions.txt \
      -sites snp-sites.txt \
      -out [pop]-freq \
      -remove_bads 1 \
      -minMapQ 30 \
      -minQ 20 \
      -nThreads 24
```

* *doMajorMinor 3: uses pre-defined major and minor alleles from the snp-sites.txt file*

This generated a [pop]-freq.mafs file per population. Custom R scripts (available upon request) were used to combine files from each population together. We retained sites if they were sampled in at least three individuals in each natural population and at least five individuals in each tail of the MAGIC population, in which 565 auxin genes were retained. To filter for a MAF > 0.05 while retaining invariant sites, we used custom R scripts to replace all allele frequencies that are less than 0.05 with 0, which was done per population or tail. File: [alleleFreqsAll .tar.gz](Data%20files/alleleFreqsAll.tar.gz)

# Genetic clustering of populations

## Phylogeny

We generated a maximum likelihood phylogeny in ```IQ-TREE v1.6.12``` [(Nguyen et al., 2015)](https://pubmed.ncbi.nlm.nih.gov/25371430/) using the polymorphisms-aware phylogenetic model. We used sites that were variable across all populations and the MAGIC population with MAF > 0.05. We used custom R scripts to thin SNPs by retaining one unlinked SNP per auxin gene – to examine robustness of results, this was repeated five times to obtain five different sets of unlinked SNPs. Results remained consistent across the different sets. Input files: [Phylogeny input](Data%20files/Phylogeny input)

We first used ```ModelFinder``` [(Kalyaanamoorthy et al., 2017)]( https://www.nature.com/articles/nmeth.4285) to determine the best-fit substitution model for the data:

```
iqtree-linux -s countsAllUnlinked-1.txt -m MF -nt 10 -pre unlinked1
```

The best-fit substitution model for the data was found to be: TVMe+FQ+P+N9.

Tree construction:

```
iqtree-linux -s countsAllUnlinked-1.txt -m TVMe+FQ+P+N9 -nt 10 -bb 10000 -alrt 10000 -pre unlinked1-run1
```

* *bb 10000: 10,000 replicates of UFboot*
* *alrt 10000: 10,000 replicates of SH-aLRT*

To assess convergence, we undertook ten separate runs of ```IQ-TREE``` for each unlinked dataset and examined tree topology (which remained unchanged across the ten independent runs). We also ensured that the log-likelihood values were stable at the end of each run.

## PCA

We explored the genetic clustering of populations using ```PCAngsd v1.10``` [(Meisner and Albrechtsen, 2018)]( https://pubmed.ncbi.nlm.nih.gov/30131346/). We used the same unlinked dataset as the phylogeny above. 

Estimate genotype likelihoods:

```
angsd	-bam bam-file-paths.txt \
	-GL 2 \
	-out allGenoLike \
	-doGlf 2 \
	-doMajorMinor 1 \
	-doMaf 1 \
	-rf unlinked-sites.txt \
	-remove_bads 1 \
	-minMapQ 30 \
	-minQ 20 \
	-nThreads 10
```

This generated the beagle file to use for the PCA.

PCA:

```
pcangsd -b allGenoLike.beagle.gz -o PCAunlinked
```
This generated the covariance matrix for the PCA (PCAunlinked.cov), which was plotted in R. R code:
```
# Read the covariance file as a matrix
cov <- as.matrix(read.table("path-to-file/PCAunlinked.cov"))
# Compute the eigenvalues and eigenvectors
eig <- eigen(C)

# Calculate the % variance the first PC explains
eig$values[1]/sum(eig$values)*100

# Read a text file to specify the colours (first column) and shape size (second column) for all individuals in the covariance file
pop <- read.table("path-to-file-/IndPopCols.txt")

# Plot the PCA
plot(eig$vectors[,1:2],col=pop[,1], pch=pop[,2], xlab="PC1",ylab="PC2")
```
# Detection of outlier SNPs

For the natural populations, we used three approaches to detect outlier SNPs separately for each locality: ```PCAngsd```, ```BayeScan```, and the top 1% of SNPs with the highest change in allele frequencies between the Dune and Headland. For the MAGIC population, we only used the top 1% of SNPs with the highest change in allele frequencies between the two tails of the gravitropic distribution. 

## PCAngsd

We undertook ```PCAngsd``` for the natural populations per locality. 

Estimate genotype likelihoods for variable sites (done per locality):

```
angsd 	-bam [locality]-bam-paths.txt \
	-GL 2 \
	-out [locality]GenoLike \
	-doGlf 2 \
	-doMajorMinor 1 \
	-SNP_pval 1e-6 \
	-doMaf 1 \
	-rf gene-regions.txt \
	-remove_bads 1 \
	-minMapQ 30 \
	-minQ 20 \
	-nThreads 10 \

* *SNP_pval 1e-6: P-value for defining a SNP*
```

Calculate selection statistics for each SNP among individuals without defining disjoint groups:

```
pcangsd -b [locality]GenoLike.beagle.gz -o [locality] --selection --sites_save
```

We used R to plot the PCAs per locality. R code (shown for LH population pair):

```
# Read the covariance file as a matrix
LH-cov <- as.matrix(read.table("path-to-file/LH.cov"))
# Compute the eigenvalues and eigenvectors
LH-eig <- eigen(C)

# Calculate the % variance the first PC explains
LH-eig$values[1]/sum(eig$values)*100

# Read a text file to specify the colours (first column) and shape size (second column) for all individuals in the covariance file
pop <- read.table("path-to-file-/LH-IndPopCols.txt")

# Plot the first 4 PCs
plot(eig$vectors[,1:2],col=pop[,1], pch=pop[,2], xlab="PC1",ylab="PC2")
plot(eig$vectors[,3:4],col=pop[,1], pch=pop[,2], xlab="PC3",ylab="PC4")
```
For each locality we visually examined which PC best separated the two populations (orange = Dune ecotype, green = Headland ecotype). 

![alt text](https://github.com/OB-lab/James-et-al.-submission-to-Current-Biology/blob/main/Images/PerLocalityPCAs.pdf?raw=true)

For LH we chose PC2, PC2 for CB, and PC1 for CH.  

We used R to calculate the outliers for each locality. R code (shown for LH population pair):

```
# Load Numpy library
library(RcppCNPy)

# Define the function for QQplot
qqchi<-function(x,...){
  lambda<-round(median(x)/qchisq(0.5,1),2)
  qqplot(qchisq((1:length(x)-0.5)/(length(x)),1),x,ylab="Observed",xlab="Expected",...);abline(0,1,col=2,lwd=2)
  legend("topleft",paste("lambda=",lambda))
}

#Set working directory
setwd("-path-to-working-directory")

# Read in selection statistics (chi2 distributed)
LH-s <- npyLoad("LH.selection.npy")

# Make QQ plot to QC the test statistics
qqchi(LH-s)

# Convert test statistics to P-values
LH-pval <- 1-pchisq(LH-s,1)

# Make a vector for the sites
LH-sites <- seq(1, dim(LH-s)[1], by=1)

# Make Manhattan plot, using PC2 (as it distinguished the ecotypes in the above PCs)
plot(-log10(LH-pval[,2])~LH-sites, xlab = "Sites",ylab = "-log10 P-value", main = "Manhattan plot")

# Put in an outlier line, that's 2.7 (corresponding to a P-value of 0.002). We considers SNPs as highly differentiated if they fell above this line
abline(2.7, 0, col="red"). 
```

Manhattan plots for each locality:

![alt text](https://github.com/OB-lab/James-et-al.-submission-to-Current-Biology/blob/main/Images/PerLocalityManhattan.pdf?raw=true)

## BayeScan

We undertook ```BayeScan v2.1``` for the natural populations per locality on sites with a MAF > 0.05. We used custom R scripts and excel to generate the BayeScan input files: [BayeScan input](Data%20files/BayeScan input)
```
bayescan_2.1 [locality]-BayeScan.txt -threads 24 -pr_odds 10 -o [locality] 
```
* *pr_odds 10: prior odds of 10, meaning the neutral model is ten times more likely than the model with selection (as recommended for datasets of our size)*

We plotted the results per population using the [(plot_R.r)]( https://github.com/mfoll/BayeScan/blob/master/R_functions/plot_R.r) R script. We categorized SNPs as highly differentiated if they contained a Q-value (the False Discovery Rate analogue of the P-value) < 0.05. R code:

```
source("path-to-file/plot_R.r")
setwd("path-to-BayeScan-output")
plot_bayescan("[locality]",FDR=0.05)
```

This produced a list of outlier SNPs per locality.

## Top 1%

We calculated the absolute allele frequency difference (|∆p|) for each allele between the Dune and Headland ecotypes at each locality and only considered sites with a MAF > 0.05 per locality. We classified SNPs as highly differentiated if they had an allele frequency difference in the top 1% quantile at each locality (corresponding to |∆p| ≥ 0.63 for Lennox Head, |∆p| ≥ 0.53 for Cabarita Beach and |∆p| ≥ 0.52 for Coffs Harbour). For the MAGIC population We calculated the absolute allele frequency difference (|∆p|) for each allele between the agravitropic and gravitropic tails, and we considered SNPs as highly differentiated if they had an allele frequency difference in the top 1% quantile (corresponding to |∆p| ≥ 0.18); we refer to these SNPs as ‘architecture SNPs’.

# Summary of outlier SNPs

For each locality of the natural populations, we considered a SNP as an outlier if it was detected as highly differentiated in at least one of the three above approaches. For the highly stringent dataset, SNPs were outliers if they were highly differentiated in at least two of the three approaches. For each locality we considered a gene as an outlier if it contained at least one outlier SNP. We compared the number of common outlier SNPs and genes across the three localities, and classified a SNP or gene as ‘parallel’ if it was an outlier in all three localities; we refer to these as ‘parallel SNPs’ and ‘parallel genes’, respectively. 

Summary file of the SNPs across the natural populations and the MAGIC population: [alleleFreqsAllSummary.tar.gz](Data%20files/alleleFreqsAllSummary.tar.gz)
Summary file of the genes across the natural populations and the MAGIC population: [genesAllSummary.txt](Data%20files/genesAllSummary.txt)
To calculate whether the number of shared outlier SNPs and genes between populations was greater than expected by chance, we used a hypergeometric distribution function, phyper, in R:

```
# Calculation of the shared outlier genes between LH and CB

# phyper (q, m, n, k...)
# q = size of overlap minus 1
# m = #outliers in group 1
# n = total # of genes minus m
# k = #outliers in group 2

phyper(85, 179, 386, 145, lower.tail=FALSE)
```

# Molecular signatures of selection

To identify molecular signatures of selection we used ```ANGSD``` to calculate Tajima’s D for each gene-region per population. The following code was run separately per population. 

Calculate the site allele frequency spectrum:

```
angsd -bam [pop]-bam-paths.txt \
	-doSaf 1 \
	-anc Senecio.contigs.fasta \
	-rf gene-regions.txt \
	-out [pop]  \
	-remove_bads 1 \
	-minMapQ 30 \
	-minQ 20 \
	-GL 2 \
	-P 24 \
	-doCounts 1 \
	-doMajorMinor 1 \
	-underFlowProtect 1 \
```

This produces a .saf.idx file. 

Obtain the folded estimate of the site allele frequency spectrum:

```
realSFS [pop].saf.idx' -P 24 -fold 1 > [pop].sfs
```

Calculate thetas for each site:

```
angsd -pest [pop].sfs \
	-bam [POP]-bam-paths.txt \
	-anc Senecio.contigs.fasta \
	-rf gene-regions.txt \
	-doThetas 1 \
	-P 24 \
	-dosaf 1 \
	-GL 2 -doMajorMinor 1 \
	-underFlowProtect 1 \
	-remove_bads 1 \
	-minMapQ 30 \
	-minQ 20 \
	-fold 1 \
	-out [pop]
```

Calculate Tajima’s D:

```
thetaStat do_stat [pop].TajD.idx 
```

Using custom R scripts we combined the output files per population. Tajima’s D values are found in the [genesAllSummary.txt](Data%20files/genesAllSummary.txt) file.

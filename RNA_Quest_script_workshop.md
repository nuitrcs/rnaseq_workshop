# RNA-seq Workshop Script
The first half of this workshop involves commands typed into the command line on Quest.  The second half is done in RStudio, on the Quest Analytics nodes.  This workshop requires that you have an account on Quest.  Information on applying for an account on Quest can be found [here](https://www.it.northwestern.edu/research/user-services/quest/allocation-guidelines.html).

## Running Commands on Quest
### Bash environment - Setup 
On your local computer, open a terminal on your local computer and [connect to Quest](https://kb.northwestern.edu///internal/page.php?id=70705) with your netID, for example from Terminal:
```
ssh YOUR_NETID@quest.it.northwestern.edu
```
Once you have logged in and are in your home directory, copy the workshop directory into your home directory: 
``` 
cd ~                                      
cp -R /projects/genomicsshare/RNAseq_workshop .   
cd RNAseq_workshop                        
```
#### Load the necessary modules
Software on Quest can either be installed locally or loaded from [system-wide modules](https://kb.northwestern.edu/page.php?id=70718).
```
module load fastqc/0.11.5;
module load hisat2/2.0.4;  
module load samtools/1.6; 
module load stringtie/1.3.4; 
```
### Step1. Analyze raw reads’ quality with FastQC  
[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is a quality control tool which provides a report on the quality of your samples.
From the command line on Quest:						
```
fastqc --outdir ./qualitycheck/ ./samples/*_chrX_*.fastq.gz 	 
```
### Step2. Filtering raw reads with Trimmomatic
[Trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic) is a flexible read trimming tool for Illumina NGS data.
```
java -jar trimmomatic-0.33.jar PE -threads 1 -phred33 ./samples/ERR188273_chrX_1.fastq.gz ./samples/ERR188273_chrX_2.fastq.gz ./ERR188273_chrX_1_paired_filtered.fastq.gz ./ERR188273_chrX_1_unpaired_filtered.fastq.gz ./ERR188273_chrX_2_paired_filtered.fastq.gz ./ERR188273_chrX_2_unpaired_filtered.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:70:20 MINLEN:30; 
```
### Step3. Re-analyze the quality of filtered reads with FastQC
Run fastqc again to confirm that this sample has enough quality to be worth studying further.
```
fastqc --outdir ./qualitycheck/ *_chrX_*_filtered.fastq.gz
```
At this point you can look at the .html file outputted by fastqc to inspect the quality control report.
### Step4. Alignment of RNA-seq reads to the genome with HISAT2
[HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml) is a fast and sensitive alignment program for mapping next-generation sequencing reads (both DNA and RNA).
```
hisat2 -p 1 --dta -x ./indexes/chrX_tran -1 ./ERR188044_chrX_1.fastq.gz -2 ./ERR188044_chrX_2.fastq.gz -S ERR188044_chrX.sam
```
### Step 5. Sort and convert the SAM file to BAM with samtools
[Samtools](https://github.com/samtools/samtools) convertx sam files into different formats, in this case a [BAM file](http://software.broadinstitute.org/software/igv/bam).
```
samtools sort -@ 1 -o ERR188044_chrX.bam ERR188044_chrX.sam
```
### Step 6. Assemble and quantify expressed genes and transcripts with StringTie
[StringTie](https://ccb.jhu.edu/software/stringtie/) is a fast and highly efficient assembler of RNA-Seq alignments into potential transcripts.
### 6-a. Stringtie assembles transcripts for each sample:
```
stringtie -p 1 -G ./genes/chrX.gtf -o ERR188044_chrX.gtf -l ERR188044 ERR188044_chrX.bam
```
### 6-b. Stringtie merges transcripts from all samples:
```
stringtie --merge -p 1 -G ./genes/chrX.gtf -o stringtie_merged.gtf mergelist.txt
```
### 6-c. Stringtie estimates transcript abundances and create table counts for Ballgown (a program downstream in the pipeline):
```
stringtie -e -B -p 1 -G stringtie_merged.gtf -o ./ballgown/ERR188044/ERR188044_chrX.gtf ERR188044_chrX.bam
```

## Submit a job to run the pipeline to this point on all the samples  
This part is done on Quest's compute nodes, where we can request multiple cores to run more threads.  Start by looking at the submission script to see how the commands we've run so far can be run on the compute nodes.  Use `msub` to submit the job to run the pipeline on all samples on the compute nodes:
```
more RNAseq_workshop_submit.sh
msub RNAseq_workshop_submit.sh
```
## Analysis and Visualization performed in RStudio
Go to: https://rstudio.questanalytics.northwestern.edu/auth-sign-in and login with your Quest credentials.  Once RStudio has started up, type:
```
setwd("/home/<YOUR_NETID>/RNAseq_workshop/") 
```
where "<YOUR_NETID>" is your Quest login netID.

### Step 7. Run  differential expression analysis with Ballgown
[Ballgown](https://bioconductor.org/packages/release/bioc/html/ballgown.html) provides flexible, isoform-level differential expression analysis.  Ballgown is part of the Bioconductor suite of tools available as R packages.

#### 7-a. RStudio environment setup (takes ~ 10 mins)
```
library("devtools") 
source("http://www.bioconductor.org/biocLite.R")
biocLite(c("alyssafrazee/RSkittleBrewer","ballgown", "genefilter","dplyr","devtools"))

library("ballgown")
library("RSkittleBrewer") # for color setup
library("genefilter") # faster calculation of mean/variance
library("dplyr") # to sort/arrange results
library("devtools")  # reproducibility/installing packages
```
#### 7-b. Load phenotype data into RStudio
```
pheno_data = read.csv("phenodata.csv")
head(pheno_data)
```
#### 7-c. Read in the expression data that were calculated by Stringtie in previous step 6-(c)
```
chrX <- ballgown(dataDir="ballgown", samplePattern="ERR", pData=pheno_data)
str(chrX)
```
#### 7-d. Filter to remove low-abundance genes 
```
chrX_filtered <- subset(chrX, "rowVars(texpr(chrX)) >1", genomesubset=TRUE)
str(chrX_filtered)
```
### Step 8. Identify transcripts/genes that show statistically significant differences with Ballgown in RStudio
#### 8-a. Identify transcripts that show statistically significant differences between groups:
```
results_transcripts <- stattest(chrX_filtered, feature="transcript", covariate="sex", adjustvars=c("condition"), getFC=TRUE, meas="FPKM")
head(results_transcripts)
```
Add gene names and gene IDs to the results:
```
results_transcripts <- data.frame(geneNames=ballgown::geneNames(chrX_filtered), geneIDs=ballgown::geneIDs(chrX_filtered), results_transcripts)
head(results_transcripts)
```
#### 8-b. Identify genes that show statistically significant differences between groups 
```
results_genes <- stattest(chrX_filtered, feature="gene", covariate="sex", adjustvars=c("condition"), getFC=TRUE, meas="FPKM")
head(results_genes)
```
### Step 9. Explore the results in RStudio
What is the most ‘differentially’ expressed transcript/genes between sexes?
Sort the results from the smallest P value to the largest:
```
results_transcripts <- results_transcripts[order(results_transcripts$pval),]
results_genes <- results_genes[order(results_genes$pval),]
```
What are the top transcript/gene expressed differently between sexes? 
```
head(results_transcripts)
head(results_genes)
```
Save the analysis results to csv files:
```
write.csv(results_transcripts, file="DifferentialExpressionAnalysis_transcript_results.csv", row.names=FALSE)
write.csv(results_genes, file="DifferentialExpressionAnalysis_gene_results.csv", row.names=FALSE)
save.image()			
```
Your workspace will be saved as '.RData' in current working directory.
### Step 10. Visualization in RStudio
#### 10-a. Plot for distribution of gene abundances across samples:
```
fpkm <- texpr(chrX, meas='FPKM')
fpkm <- log2(fpkm +1)
boxplot(fpkm, col=as.numeric(pheno_data$sex), las=2,ylab='log2(FPKM+1)')
```
#### 10-b. Plot for individual expression of a certain transcript between groups: 
Setup palette with your favorite colors
```
coloring <- c('darkgreen', 'skyblue', 'hotpink', 'orange', 'lightyellow')
palette(coloring)
```
#### In this example, by looking head(results_transcripts), I randomly choose to draw the 13th most differientially expressed transcript. (gene name "XIST") You can also decide the transcript/gene of your interest. What you need to know is its genename or transcript name! 
```
which(ballgown::geneNames(chrX)=="XIST")	
ballgown::transcriptNames(chrX)[1484]		
plot(fpkm[1484,] ~ pheno_data$sex, border=c(1,2), main=paste(ballgown::geneNames(chrX)[1484], ' : ',ballgown::transcriptNames(chrX)[1484]), pch=19, xlab="sex", ylab='log2(FPKM+1)')
points(fpkm[1484,] ~ jitter(as.numeric(pheno_data$sex)), col=as.numeric(pheno_data$sex))
```
#### 10-c/d. Plot the average expression levels for all transcripts of a gene within different groups:
```
geneIDs(chrX)[1484] 
plotMeans('MSTRG.495', chrX_filtered, groupvar="sex", legend=FALSE)
plotMeans(ballgown::geneIDs(chrX)[1484], chrX, groupvar="sex", legend=FALSE)
```
Play with it! 




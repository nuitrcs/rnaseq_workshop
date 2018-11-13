#RNA-seq Workshop Script

### Bash environment - Setup 

### * It is highly enouraged to install the required packages before attending the workshop
On your local computer, open a terminal on your local computer and connect to Quest with your netID (ssh YOUR_NETID@quest.it.northwestern.edu), then copy the workshop directory into your home directory: 
``` 
cd ~                                      
cp -R /projects/genomicsshare/RNAseq_workshop .   
cd RNAseq_workshop                        
```
### Load the necessary modules
```
module load fastqc/0.11.5;
module load hisat2/2.0.4;  
module load samtools/1.6; 
module load stringtie/1.3.4; 
```
# Step1. Analyze raw reads’ quality with FastQC  
### STEPS 1 - 3 DEMO ONLY due to time constraints
### Bash environment 						
```
fastqc --outdir ./qualitycheck/ ./samples/*_chrX_*.fastq.gz 	 
```
### <Step4. Alignment of RNA-seq reads to the genome with HISAT>
```
hisat2 -p 1 --dta -x ./indexes/chrX_tran -1 ./ERR188044_chrX_1.fastq.gz -2 ./ERR188044_chrX_2.fastq.gz -S ERR188044_chrX.sam
```

### <Step 5. Sort and convert the SAM file to BAM with samtools>
```
samtools sort -@ 1 -o ERR188044_chrX.bam ERR188044_chrX.sam
```

### <Step 6. Assemble and quantify expressed genes and transcripts with StringTie>
### 6-a. Stringtie assembles transcripts for each sample:
```
stringtie -p 1 -G ./genes/chrX.gtf -o ERR188044_chrX.gtf -l ERR188044 ERR188044_chrX.bam
```

### 6-b. Stringtie merges transcripts from all samples:
```
stringtie --merge -p 1 -G ./genes/chrX.gtf -o stringtie_merged.gtf mergelist.txt
```
### 6-c. Stringtie estimates transcript abundances and create table counts for Ballgown:
```
stringtie -e -B -p 1 -G stringtie_merged.gtf -o ./ballgown/ERR188044/ERR188044_chrX.gtf ERR188044_chrX.bam
```
## Run the pipeline to this point on all the samples  
This part is done on the compute nodes, where we can request multiple cores to run more threads.  Start by looking at the submission script to see how the commands we've run so far can be run on the compute nodes.  Use msub to submit the job to run the pipeline on all samples on the compute nodes:
```
more RNAseq_workshop_submit.sh
msub RNAseq_workshop_submit.sh
```
### VISUALIZATION 
Go to: https://rstudio.questanalytics.northwestern.edu/auth-sign-in and login with your Quest credentials.  Once RStudio has started up, type:
```
setwd("/home/<YOUR_NETID>/RNAseq_workshop/") 
```
where <YOUR_NETID> is your Quest login netID.

### <Step 7. Run  differential expression analysis with Ballgown>

### 7-a. RStudio environment setup (takes ~ 10 mins)
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
### 7-b. Load phenotype data into RStudio
```
pheno_data = read.csv("phenodata.csv")
head(pheno_data)
```
### 7-c. Read in the expression data that were calculated by Stringtie in previous step 6-(c)
```
chrX <- ballgown(dataDir="ballgown", samplePattern="ERR", pData=pheno_data)
str(chrX)
```
### 7-d. Filter to remove low-abundance genes 
```
chrX_filtered <- subset(chrX, "rowVars(texpr(chrX)) >1", genomesubset=TRUE)
str(chrX_filtered)
```
### Step 8. Identify transcripts/genes that show statistically significant differences with Ballgown in RStudio
### 8-a. Identify transcripts that show statistically significant differences between groups:
```
results_transcripts <- stattest(chrX_filtered, feature="transcript", covariate="sex", adjustvars=c("condition"), getFC=TRUE, meas="FPKM")
head(results_transcripts)
```
Add gene names and gene IDs to the results:
```
results_transcripts <- data.frame(geneNames=ballgown::geneNames(chrX_filtered), geneIDs=ballgown::geneIDs(chrX_filtered), results_transcripts)
head(results_transcripts)
```
### 8-b. Identify genes that show statistically significant differences between groups 
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
save.image()			# your workspace will be saved as '.RData' in current working directory
```
### Step 10. Visualization in RStudio
### 10-a. Plot for distribution of gene abundances across samples:
```
fpkm <- texpr(chrX, meas='FPKM')
fpkm <- log2(fpkm +1)
boxplot(fpkm, col=as.numeric(pheno_data$sex), las=2,ylab='log2(FPKM+1)')
```
### 10-b. Plot for individual expression of a certain transcript between groups: 
### Setup palette with your favorite colors
```
coloring <- c('darkgreen', 'skyblue', 'hotpink', 'orange', 'lightyellow')
palette(coloring)
```
### In this example, by looking head(results_transcripts), I randomly choose to draw the 13th most differientially expressed transcript. (gene name "XIST") You can also decide the transcript/gene of your interest. What you need to know is its genename or transcript name! 
```
which(ballgown::geneNames(chrX)=="XIST")	# 1484
ballgown::transcriptNames(chrX)[1484]		# NR_001564
plot(fpkm[1484,] ~ pheno_data$sex, border=c(1,2), main=paste(ballgown::geneNames(chrX)[1484], ' : ',ballgown::transcriptNames(chrX)[1484]), pch=19, xlab="sex", ylab='log2(FPKM+1)')
points(fpkm[1484,] ~ jitter(as.numeric(pheno_data$sex)), col=as.numeric(pheno_data$sex))
```
### 10-c/d. Plot the average expression levels for all transcripts of a gene within different groups:
```
geneIDs(chrX)[1484] # MSTRG.495
plotMeans('MSTRG.495', chrX_filtered, groupvar="sex", legend=FALSE)
plotMeans(ballgown::geneIDs(chrX)[1484], chrX, groupvar="sex", legend=FALSE)
```
Play with it! 



#!/bin/bash
#MSUB -A b1042
#MSUB -q genomics
#MSUB -l walltime=08:00:00
#MSUB -l nodes=1:ppn=6
#MSUB -N RNAseq_workshop
 
WORKDIR=/home/$USER
cp -r /projects/genomicsshare/RNAseq_workshop $WORKDIR/.
WORKDIR=/home/$USER/RNAseq_workshop
cd $WORKDIR

### Load the necessary modules 
module load fastqc/0.11.5
module load hisat2/2.0.4  
module load samtools/1.6 
module load stringtie/1.3.4 

### Step1. Analyze raw reads’ quality with FastQC 
echo "Step 1";
fastqc --outdir $WORKDIR/qualitycheck/ $WORKDIR/samples/*_chrX_*.fastq.gz

### Step2. Filtering raw reads with Trimmomatic
SAMPLES=('ERR188273' 'ERR188044' 'ERR188104' 'ERR188234' 'ERR188454' 'ERR204916');
echo "Step 2";
for sample in "${SAMPLES[@]}"; do
    java -jar trimmomatic-0.33.jar PE -threads 6 -phred33 $WORKDIR/samples/${sample}_chrX_1.fastq.gz $WORKDIR/samples/${sample}_chrX_2.fastq.gz $WORKDIR/${sample}_chrX_1_paired_filtered.fastq.gz $WORKDIR/${sample}_chrX_1_unpaired_filtered.fastq.gz $WORKDIR/${sample}_chrX_2_paired_filtered.fastq.gz $WORKDIR/${sample}_chrX_2_unpaired_filtered.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:70:20 MINLEN:30;
done

### Step3. Re-analyze the quality of filtered reads with FastQC
echo "Step 3";
fastqc --outdir $WORKDIR/qualitycheck/ $WORKDIR/*_chrX_*_filtered.fastq.gz

### Step4. Alignment of RNA-seq reads to the genome with HISAT
echo "Step 4";
for sample in "${SAMPLES[@]}"; do
    hisat2 -p 6 --dta -x $WORKDIR/indexes/chrX_tran -1 $WORKDIR/${sample}_chrX_1_paired_filtered.fastq.gz -2 $WORKDIR/${sample}_chrX_2_paired_filtered.fastq.gz -S $WORKDIR/${sample}_chrX.sam;
done

### Step 5. Sort and convert the SAM file to BAM with samtools
echo "Step 5";
for sample in "${SAMPLES[@]}"; do
    samtools sort -@ 6 -o $WORKDIR/${sample}_chrX.bam $WORKDIR/${sample}_chrX.sam;
done


### Step 6. Assemble and quantify expressed genes and transcripts with tringTie
#### 6-a. Stringtie assembles transcripts for each sample:
echo "Step 6a";
for sample in "${SAMPLES[@]}"; do
    stringtie -p 6 -G $WORKDIR/genes/chrX.gtf -o $WORKDIR/${sample}_chrX.gtf -l $WORKDIR/${sample} $WORKDIR/${sample}_chrX.bam;
done

#### 6-b. Stringtie merges transcripts from all samples:
echo "Step 6b";
stringtie --merge -p 6 -G $WORKDIR/genes/chrX.gtf -o $WORKDIR/stringtie_merged.gtf $WORKDIR/mergelist.txt

#### 6-c. Stringtie estimates transcript abundances and create table counts for Ballgown
echo "Step 6c";
for sample in "${SAMPLES[@]}"; do
    stringtie -e -B -p 6 -G $WORKDIR/stringtie_merged.gtf -o $WORKDIR/ballgown/${sample}/${sample}_chrX.gtf $WORKDIR/${sample}_chrX.bam;
done

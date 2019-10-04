#!/bin/bash

printf "\n***Description: This script takes fastq files of RNA-Seq data in the current directory, runs STAR alignment and FastQC, and outputs a counts matrix for all.***\n"

# Exit if an error occurs
set -e

# To debug, display the command before executing it by "set -x"; to turn this option off by "set +x"

# Decide on threads
printf "\n=========================================================================================\n"
read -p "***How many threads do you want to use?*** Answer: " nThread
echo "***Your input thread is $nThread***"
printf "=========================================================================================\n"

# Create directories
# -p option will make sure that mkdir will create the directory only if it does not exist, and it won’t throw an error if it does exist
mkdir -p starOutput countOutput fastqcOutput

printf "=========================================================================================\n"
echo "STAR"
printf "=========================================================================================\n"

# STAR alignment
STARindex=/data/Genomes/mouse/mm10/star/index101bp/
vds2Dir=/data/home/lingzili/VDS2/NGS_Runs/190627_A00316_0084_BHCJ7LDRXX/Data/Intensities/BaseCalls/Demultiplex/SMB_KR

ls ${vds2Dir} | cut -d"_" -f 1,2 | sort -u | while read id
do 
	echo "***The sample ID is $id***"

	# Set up output filenames and locations
	sortedByCoordBAM=./starOutput/${id}_Aligned.sortedByCoord.out.bam
	sortedByNameBAM=./starOutput/${id}_Aligned.out.bam
	indexedBAM=./starOutput/${id}_Aligned.sortedByCoord.out.bam.bai
	statBAM=./starOutput/${id}.stat.txt
	read1=/data/home/lingzili/VDS2/NGS_Runs/190627_A00316_0084_BHCJ7LDRXX/Data/Intensities/BaseCalls/Demultiplex/SMB_KR/${id}_R1_001.fastq.gz
	read2=/data/home/lingzili/VDS2/NGS_Runs/190627_A00316_0084_BHCJ7LDRXX/Data/Intensities/BaseCalls/Demultiplex/SMB_KR/${id}_R2_001.fastq.gz

	echo "***STAR alignment of ${read1} and ${read2}***"

	/data/home/lingzili/.conda/envs/lingzili/bin/STAR --genomeDir $STARindex \
	--runThreadN $nThread \
	--readFilesIn $read1 $read2 \
	--readFilesCommand zcat \
	--outFileNamePrefix ./starOutput/${id}_ \
	--outBAMsortingThreadN $nThread \
	--outSAMtype BAM Unsorted SortedByCoordinate
	echo "***${id} is aligned***"
	echo "***Note: The unsorted BAM file from STAR is actually sorted by name***"

	echo "***Start index BAM file of ${id}***"
	samtools index -@ $nThread $sortedByCoordBAM $indexedBAM

	echo "***Start extract statistics from BAM file of ${id}***"
	samtools flagstat -@ $nThread $sortedByCoordBAM >$statBAM
done

printf "=========================================================================================\n"
echo "featureCounts"
printf "=========================================================================================\n"

geneAnnotation=/data/Genomes/mouse/mm10/mm10.refGene.gtf

featureCounts -T $nThread --verbose -p -a $geneAnnotation -o ./countOutput/featureCount ./starOutput/*.sortedByCoord.out.bam

# Quality control with FastQC and MultiQC
printf "=========================================================================================\n"
echo "FastQC"
printf "=========================================================================================\n"

# Disable pop-up display
unset DISPLAY

# Run FastQC on fastq and BAM files
fq=/data/home/lingzili/VDS2/NGS_Runs/190627_A00316_0084_BHCJ7LDRXX/Data/Intensities/BaseCalls/Demultiplex/SMB_KR/*
BAM=./starOutput/*.sortedByCoord.out.bam

fastqc -t $nThread -o fastqcOutput ${fq} ${BAM}

# multiqc to summarize all stats
printf "=========================================================================================\n"
echo "MultiQC"
printf "=========================================================================================\n"

/usr/local/bin/multiqc -o fastqcOutput .

# Extract information from count table
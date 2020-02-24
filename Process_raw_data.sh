#!/bin/bash

#Shell script that performs quality checking of the reads before and after trimming, maps the reads onto reference genome generating indexed BAM files
#Requires: fastqc, trimmomatic, bowtie2, samtools

WD="~/data/DroAr_data/"
echo $WD
cd $WD

#Path to file with sequencing adapters
Adapters="TruSeq3-SE.fa"
RefSeqs="refseqs/Cba_and_phage.fa"
RefSeqs_index="refseqs/index/Cba_and_phage_index"

mkdir QC
mkdir QC/Trimmed

fastqc -t 4 -o QC Data/*

mkdir Data/Trimmed

#Removing adapter sequences and low-quality reads fragments
for f in Data/*fastq.gz;
do
	echo $f;
	java -jar ~/Trimmomatic-0.38/trimmomatic-0.38.jar SE -phred33 $f Data/Trimmed/`basename $f .fastq.gz`_trimmed.fastq.gz ILLUMINACLIP:$Adapters:2:30:10 LEADING:0 TRAILING:0 SLIDINGWINDOW:4:15 MINLEN:36;
done

fastqc -t 4 -o QC/Trimmed Data/Trimmed/*

#Mapping reads onto reference sequences
mkdir refseqs/index
mkdir alignments

#Build bowtie2 index
bowtie2-build $RefSeqs $RefSeqs_index

for f in Data/Trimmed/*fastq.gz; 
do 
	echo $f; 
	bowtie2 -x $RefSeqs_index -U $f -q -p 4 | samtools view -S -b -u - | samtools sort - -o alignments/`basename $f .fastq.gz`.sorted; 
done

#Count mapped reads for each sample
for f in alignments/*bam; 
do
	printf "%s\t%s\n" $f `samtools view -F 4 -c $f` >> reads_per_lib.tsv; 
done

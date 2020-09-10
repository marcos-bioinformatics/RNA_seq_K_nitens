## This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International$
## Author: Marcos Elizalde Horcada
## Date: February 2020
## Contact: Francisco José Romero Campero - fran@us.es

#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -j yes

SAMPLE_FOLDER=$1
ACCESSION=$2
SAMPLE_NUMBER=$3

#Example: qsub -N sampleN -o sampleN.bam sample_processing.sh DIRECTORY ACCESSION_NUMBER SAMPLE_NUMBER

# Downloading sample file
cd ${SAMPLE_FOLDER} 
fastq-dump --split-files ${ACCESSION} 

# Sample quality control and read mapping to reference genome
if [ -f ${ACCESSION}_2.fastq ] 
then
   fastqc ${ACCESSION}_1.fastq
   fastqc ${ACCESSION}_2.fastq

   hisat2 --dta -x ../../annotation/index -1 ${ACCESSION}_1.fastq -2 ${ACCESSION}_2.fastq -S sample_${SAMPLE_NUMBER}.sam
else
   fastqc ${ACCESSION}_1.fastq

   hisat2 --dta -x ../../annotation/index -U ${ACCESSION}_1.fastq -S sample_${SAMPLE_NUMBER}.sam
fi

# Generating sorted bam file
samtools sort -o sample_${SAMPLE_NUMBER}.bam sample_${SAMPLE_NUMBER}.sam
rm sample_${SAMPLE_NUMBER}.sam
rm *.fastq			
samtools index sample_${SAMPLE_NUMBER}.bam
bamCoverage -bs 10 --normalizeUsing CPM --bam sample_${SAMPLE_NUMBER}.bam -o sample_${SAMPLE_NUMBER}.bw

# Transcript assembly
stringtie -G ../../annotation/annotation.gtf -o sample_${SAMPLE_NUMBER}.gtf -l sample_${SAMPLE_NUMBER} sample_${SAMPLE_NUMBER}.bam

# Preparing merge list file for transcriptome merging
echo ${SAMPLE_FOLDER}/sample_${SAMPLE_NUMBER}.gtf >> ../../results/merge_list.txt

# Gene Expression Quantification
stringtie -e -B -G ../../annotation/annotation.gtf -o sample_${SAMPLE_NUMBER}.gtf sample_${SAMPLE_NUMBER}.bam


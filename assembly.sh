#!/bin/bash
path="/home/eppicenter/Downloads/kleb_pneumoniae/MARINA/Klebfastq/neonatal1/fastq/noassemblies/fastq/1"

reads="${path}/fastq"
cd ${reads}
mkdir -p ${path}/assembly
results="${path}/assembly"
mkdir -p ${results}/trim
mkdir -p ${results}/all

for sample in `cat ${reads}/list.txt`
do
    R1=`echo ${sample}_*R1*fastq.gz`
    R2=`echo ${sample}_*R2*fastq.gz`
    trimmomatic PE -threads 24 ${reads}/$R1 ${reads}/$R2\
		-trimlog trim.log -baseout ${results}/trim/${sample}.fastq.gz\
		ILLUMINACLIP:/home/eppicenter/Downloads/kleb_pneumoniae/MARINA/kleb_07_04/Fastq/Sequencing_adaptors.fasta:2:30:10\
		LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    
    mkdir -p ${results}/${sample}
    unicycler -1 ${results}/trim/${sample}_1P.fastq.gz\
	       -2 ${results}/trim/${sample}_2P.fastq.gz\
	       -o ${results}/${sample} --keep 0 --mode conservative -t 24
    
    cp ${results}/${sample}/assembly.fasta ${results}/all/${sample}.fasta
    

   # skesa --cores 4 --memory 10 --reads $R1,$R2 > ${sample}.fasta
done


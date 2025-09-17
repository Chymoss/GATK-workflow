#!/bin/bash 
#SBATCH -J BATCH
#SBATCH -p Acluster 
#SBATCH -n 4
#SBATCH --output=%j.out 
#SBATCH --error=%j.err

#set environment
source /Share/samples/Acluster.sh
source /Share/home/dengzg1993/bin/DZG_env.sh

#my_code
#Parameter(variable)
CURRdir=/Share/home/dengzg1993/Analyses/RNASeq/Arabidopsis/WJJ_RNA_seq/RNAseq_hormone_20240326/Batch-20240326
Reference=/Data/ChenhaodongLab/Data_deposited/Reference_genome/Arabidopsis/TAIR10/assembly/Ath
GTFfile=/Data/ChenhaodongLab/Data_deposited/Reference_genome/Arabidopsis/TAIR10/assembly/Arabidopsis_thaliana.TAIR10.31.gtf
Thread=4

#Main
cd $CURRdir
mkdir fastqc
mkdir export
for x in $(cat sample.list | cut -f 1)
do
	Sample=${x##*/}
	Filename=${x%/*}
	READ1=/Data/ChenhaodongLab/Data_deposited/RNASeq/RNASeq/Arabidopsis/WJJ_RNA_seq/X101SC23105095-Z01-J073/${Filename}_good_1.fq.gz
	READ2=/Data/ChenhaodongLab/Data_deposited/RNASeq/RNASeq/Arabidopsis/WJJ_RNA_seq/X101SC23105095-Z01-J073/${Filename}_good_2.fq.gz
	
	fastqc --noextract -o ./fastqc ${READ1} ${READ2} 2> ./fastqc/${Sample}_qcreport.txt &
	
	hisat2 -p ${Thread} -t -x ${Reference} \
	-1 ${READ1} -2 ${READ2} \
	-S ./${Sample}_mapping.sam
	
	samtools view -b -S ${Sample}_mapping.sam > ${Sample}_mapping.bam
	
	samtools sort ${Sample}_mapping.bam -o ${Sample}_sort.bam
	
	samtools index ${Sample}_sort.bam
	
	FILE_SIZE=`du ${Sample}_sort.bam | awk '{print $1}'`
	
	if [ $FILE_SIZE -ge 1000000 ]
	then
		rm ${Sample}_mapping.sam
		rm ${Sample}_mapping.bam
	fi
	
	htseq-count \
	--format=bam \
	--order=pos \
	--stranded=no \
	-t exon \
	-i gene_id \
	${Sample}_sort.bam ${GTFfile} > ${Sample}_htseqreport-gene_id_no.txt
	
	cp ${Sample}_htseqreport-gene_id_no.txt $export
	
done



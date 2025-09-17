#!/bin/bash 
#SBATCH -J 5uM-Pic-3h-1
#SBATCH -p Acluster 
#SBATCH -n 4
#SBATCH --output=%j.out 
#SBATCH --error=%j.err

#set environment
source /Share/samples/Acluster.sh
source /Share/home/dengzg1993/bin/DZG_env.sh

#my_code
#Parameter(variable)
CURRdir=/Share/home/dengzg1993/Analyses/RNASeq/Arabidopsis/WJJ_RNA_seq/RNAseq_hormone_20240326/DMSO/5uM-Pic-3h-1
Sample=5uM-Pic-3h-1

READ1=/Data/ChenhaodongLab/Data_deposited/RNASeq/RNASeq/Arabidopsis/WJJ_RNA_seq/X101SC23105095-Z01-J073/3-P-1_1.fq.gz
READ2=/Data/ChenhaodongLab/Data_deposited/RNASeq/RNASeq/Arabidopsis/WJJ_RNA_seq/X101SC23105095-Z01-J073/3-P-1_2.fq.gz

Reference=/Data/ChenhaodongLab/Data_deposited/Reference_genome/Arabidopsis/TAIR10/assembly/Ath
Splicesites=/Data/ChenhaodongLab/Data_deposited/Reference_genome/Arabidopsis/TAIR10/assembly/Splicesite.txt
GTFfile=/Data/ChenhaodongLab/Data_deposited/Reference_genome/Arabidopsis/TAIR10/assembly/Arabidopsis_thaliana.TAIR10.31.gtf

#Parameter(fixed)
fullname1=${READ1##*/}
filename1=${fullname1%.fq*}
READ_v1="${filename1}_val_1.fq.gz"

fullname2=${READ2##*/}
filename2=${fullname2%.fq*}
READ_v2="${filename2}_val_2.fq.gz"

Thread=4

#preprocess
cd ${CURRdir}
mkdir Trimmed
cd Trimmed

echo " trim_galore cut adapters started at $(date)"
trim_galore --length 100 --paired ${READ1} ${READ2} --gzip -o ${CURRdir}/TrimReads        
echo "trim_galore cut adapters finished at $(date)"

#QC
mkdir fastqc
fastqc --noextract -o ./fastqc ${READ_v1} ${READ_v2}  2> ./fastqc/${Sample}_qcreport_Trim.txt &

#align
hisat2 -p ${Thread} -t -x ${Reference} \
-1 ${READ_v1} -2 ${READ_v2} \
--known-splicesite-infile ${Splicesites} \
-S ./${Sample}_mapping.sam 2> runreport.txt

samtools view -b -S ${Sample}_mapping.sam > ${Sample}_mapping.bam
samtools sort ${Sample}_mapping.bam -o ${Sample}_sort.bam
samtools index ${Sample}_sort.bam

FILE_SIZE=`du ${Sample}_sort.bam | awk '{print $1}'`
if [ $FILE_SIZE -ge 1000000 ]
then
	rm ${Sample}_mapping.sam
	rm ${Sample}_mapping.bam
fi

#Reads-counting 如果确定方法，就不用跑3个程序(链特异性建库分为 --stranded = reverse，stranded，非特异建库对应 no，一般都是非特异建库)
htseq-count \
--format=bam \
--order=pos \
--stranded=no \
-t exon \
-i gene_id \
${Sample}_sort.bam ${GTFfile} > ${Sample}_htseqreport-gene_id_no.txt

#cat > JobDone.txt << END_TEXT
#Job is done.
#END_TEXT


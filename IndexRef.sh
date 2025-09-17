#!/bin/bash 
#SBATCH -J Idx
#SBATCH -p Acluster
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --output=%j.out
#SBATCH --error=%j.err

#set environment
source /Share/samples/Acluster.sh
source /Share/home/dengzg1993/bin/DZG_env.sh


###################################
#####  0.Parameters
#Dir
WORKDIR=/Data/ChenhaodongLab/Data_deposited/Reference_genome/Ceratodon/507/customer_backup/02.assemble
Species=PP

#Tools
RefFasta=/Data/ChenhaodongLab/Data_deposited/Reference_genome/Ceratodon/507/customer_backup/02.assemble/final.genome.fasta
RefGff3=/Data/ChenhaodongLab/Data_deposited/Reference_genome/Ceratodon/507/customer_backup/03.annotation/02.gene_prediction/EVM.final.all.transcripts.gene.gff3
my_picard=/Share/home/dengzg1993/tools/picard.jar
my_ANOV=/Share/home/dengzg1993/tools/annovar

#bwa index
cd $WORKDIR
mkdir Anovar_db
mkdir Blast_db

samtools faidx $RefFasta
java -jar $my_picard CreateSequenceDictionary -R $RefFasta
bwa index $RefFasta
bwa-mem2 index $RefFasta
#annovar index
cd $WORKDIR/Anovar_db

gff3ToGenePred $RefGff3 ${Species}_refGene.txt
$my_ANOV/retrieve_seq_from_fasta.pl --format refGene --seqfile $RefFasta ${Species}_refGene.txt --out ${Species}_refGeneMrna.fa 

#blast Index
cd $WORKDIR/Blast_db
makeblastdb -dbtype nucl -in $RefFasta -out ${Species}.blastdb

#hisat2 index
cd $WORKDIR
gffread $RefGff3 -T -o $Species.gtf
hisat2_extract_splice_sites.py $Species.gtf > Splicesite.txt
hisat2_extract_exons.py $Species.gtf > Exon.txt
hisat2-build -f --ss Splicesite.txt --exon Exon.txt $RefFasta $Species



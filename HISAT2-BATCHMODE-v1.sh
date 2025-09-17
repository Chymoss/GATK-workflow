#!/bin/bash
#SBATCH -J GATK_BATCH
#SBATCH -p Acluster
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --array=1-12
#SBATCH --output=%j.out
#SBATCH --error=%j.err

# set environment
source /Share/samples/Acluster.sh
source /Share/home/dengzg1993/bin/DZG_env.sh

# Parameters
WORKDIR=$(pwd)
Reference=/Data/ChenhaodongLab/Data_deposited/Reference_genome/Arabidopsis/TAIR10/assembly/Ath
GTFfile=/Data/ChenhaodongLab/Data_deposited/Reference_genome/Arabidopsis/TAIR10/assembly/Arabidopsis_thaliana.TAIR10.31.gtf

#############################################
# 函数：处理单个样本
#############################################
process_sample() {
    local sample=$1
    local read1=$2
    local read2=$3

    local sample_dir=$WORKDIR/$sample
    mkdir -p $sample_dir/report $sample_dir/1.mapping $sample_dir/2.output

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Start processing $sample"

    # FastQC
    fastqc --noextract -o $sample_dir/report $read1 $read2 > $sample_dir/report/${sample}_fastqc.log 2>&1

    # Mapping & sorting
    hisat2 -p 4 -t -x $Reference -1 $read1 -2 $read2 \
        | samtools view -Sb - \
        | samtools sort -o $sample_dir/1.mapping/${sample}.sorted.bam

    local bam=$sample_dir/1.mapping/${sample}.sorted.bam

    # HTSeq-count
    htseq-count --format=bam --order=pos --stranded=no -t exon -i gene_id $bam $GTFfile \
        > $sample_dir/2.output/${sample}_htseqreport-gene_id_no.txt

    if [ $? -eq 0 ]; then
        # 删除大文件以节省空间
        local file_size=$(du -k $bam | awk '{print $1}')
        if [ "$file_size" -ge 1000000 ]; then
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] File size >= 1GB, deleting $bam"
            rm $bam
        else
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] File size < 1GB, keeping $bam"
        fi
    else
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] HTSeq-count failed for $sample"
    fi

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Finished $sample"
}

#############################################
# 主逻辑：读取 array job 对应的行
#############################################
dos2unix sample_list.tsv
TASK_ID=${SLURM_ARRAY_TASK_ID:-1}
LINE=$(sed -n "${TASK_ID}p" sample_list.tsv)

SAMPLE=$(echo $LINE | awk '{print $1}')
READ1=$(echo $LINE | awk '{print $2}')
READ2=$(echo $LINE | awk '{print $3}')

process_sample $SAMPLE $READ1 $READ2

#!/bin/bash
#SBATCH -J GTRC
#SBATCH -p Acluster
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --output=%j.out
#SBATCH --error=%j.err


#set environment
source /Share/samples/Acluster.sh
source /Share/home/dengzg1993/bin/DZG_env.sh

# 参数设置
#############################################
WORKDIR=/Share/home/dengzg1993/Analyses/DNASeq/Physcomitrella/2.Pool_sequencing/GTRC_pool/2024.12.10
TMPDIR=/Share/home/dengzg1993/tmp

SAMPLE_MUT=DONG-pool
MUT_PLOIDY=48
READ1_MUT=/Data/ChenhaodongLab/Data_deposited/DNASeq/Physcomitrella/2.Pool-seq/GTRC-F2pool-mut/clean_data/DONG-pool_R1.fq.gz
READ2_MUT=/Data/ChenhaodongLab/Data_deposited/DNASeq/Physcomitrella/2.Pool-seq/GTRC-F2pool-mut/clean_data/DONG-pool_R2.fq.gz

SAMPLE_WT=WT-pool
WT_PLOIDY=48
READ1_WT=/Data/ChenhaodongLab/Data_deposited/DNASeq/Physcomitrella/2.Pool-seq/GTRC-F2pool-WT/clean_data/WT-pool_R1.fq.gz
READ2_WT=/Data/ChenhaodongLab/Data_deposited/DNASeq/Physcomitrella/2.Pool-seq/GTRC-F2pool-WT/clean_data/WT-pool_R2.fq.gz

DEPTH=50

REF_FASTA=/Data/ChenhaodongLab/Data_deposited/Reference_genome/Physcomitrella/Gd_v6.1/Ppatens_870_V6.fasta
SNP_SITES=/Share/home/dengzg1993/Analyses/DNASeq/Physcomitrella/1.WGS/05.Villersexel/2024.11.25_v6.1/4.annotate/Results/Vx.cohort.snp.final.vcf.gz
#############################################



# 初始化工作目录
#############################################
cd $WORKDIR || exit
mkdir -p report 1.mapping 2.variants 3.output

# 记录时间
log_step() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}
log_step "开始任务"
#############################################



# 快速质量检查 (FastQC)
#############################################
log_step "执行 FastQC 检查"
fastqc -o report --noextract $READ1_MUT $READ2_MUT &> report/fastqc_m.log &
fastqc -o report --noextract $READ1_WT $READ2_WT &> report/fastqc_w.log &
wait
#############################################


# 映射和去重 (Mapping and Deduplication)
#############################################
log_step "开始 Reads 映射和去重"
cd $WORKDIR/1.mapping || exit
process_sample() {
    local sample_name=$1
    local read1=$2
    local read2=$3

    bwa mem -t 4 -M \
        -R "@RG\tID:$sample_name\tLB:$sample_name\tPL:ILLUMINA\tPM:HISEQ\tSM:$sample_name" \
        $REF_FASTA $read1 $read2 \
        | samtools view -Sb - \
        | samtools sort -o $sample_name.sorted.bam

    picard MarkDuplicates \
        -I $sample_name.sorted.bam \
        -M sorted_dedup_metrics-${sample_name}.txt \
        --TMP_DIR $TMPDIR \
        -O $sample_name.sorted_dedup_reads.bam

    samtools index $sample_name.sorted_dedup_reads.bam

    # 清理空间
    FILE_SIZE=$(du -k $sample_name.sorted_dedup_reads.bam | awk '{print $1}')
    if [ $FILE_SIZE -ge 1000000 ]; then
        rm $sample_name.sorted.bam
    fi
}
process_sample $SAMPLE_MUT $READ1_MUT $READ2_MUT &
process_sample $SAMPLE_WT $READ1_WT $READ2_WT &
wait
#############################################


# 变异检测和注释 (Variant Calling and Annotation)
#############################################
log_step "执行变异检测"
cd $WORKDIR/2.variants || exit
variant_calling() {
    local sample_name=$1
    local ploidy=$2

    gatk HaplotypeCaller -R $REF_FASTA \
        -I $WORKDIR/1.mapping/$sample_name.sorted_dedup_reads.bam \
        --intervals $SNP_SITES --sample-ploidy $ploidy -ERC GVCF \
        -O $sample_name.g.vcf.gz
}
variant_calling $SAMPLE_MUT $MUT_PLOIDY &
variant_calling $SAMPLE_WT $WT_PLOIDY &
wait
log_step "合并变异文件"
gatk CombineGVCFs -R $REF_FASTA \
    --variant $SAMPLE_MUT.g.vcf.gz \
    --variant $SAMPLE_WT.g.vcf.gz \
    -O Combine.g.vcf.gz
gatk GenotypeGVCFs -R $REF_FASTA -V Combine.g.vcf.gz -O Combine.vcf.gz
mv Combine.vcf.gz $WORKDIR/3.output
mv Combine.vcf.gz.tbi $WORKDIR/3.output
#############################################

# 提取等位基因频率 (Extract Allele Frequency)
#############################################
cd $WORKDIR/3.output

gatk VariantFiltration -R $REF_FASTA -V Combine.vcf.gz \
    -O Combine.filter.vcf.gz \
    -filter-name "hard_filter" \
    -filter "QD < 2.0 || SOR > 3.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"

gatk SelectVariants --exclude-filtered --select-type-to-include SNP \
    --restrict-alleles-to BIALLELIC \
    -select "vc.getGenotype('$SAMPLE_MUT').getDP() > 0.3 * $DEPTH " \
    -select "vc.getGenotype('$SAMPLE_MUT').getDP() < 3 * $DEPTH " \
	-select "vc.getGenotype('$SAMPLE_WT').getDP() > 0.3 * $DEPTH " \
    -select "vc.getGenotype('$SAMPLE_WT').getDP() < 3 * $DEPTH " \
    -V Combine.filter.vcf.gz \
    -O Combine.snp.cleaned.vcf.gz

gatk VariantsToTable -V Combine.snp.cleaned.vcf.gz \
    -F CHROM -F POS -GF AD -O Combine.results_var_AD.table
	
log_step "任务结束"
#############################################
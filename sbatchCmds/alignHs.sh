#!/bin/sh
#SBATCH --account=gdkendalllab
#SBATCH --array=0-4
#SBATCH --error=slurmOut/alignHs-%j.txt
#SBATCH --output=slurmOut/alignHs-%j.txt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --job-name alignHs
#SBATCH --wait
#SBATCH --partition=general,himem
#SBATCH --mail-user=matthew.cannon@nationwidechildrens.org
#SBATCH --mail-type=FAIL,REQUEUE,TIME_LIMIT_80
#SBATCH --time=2-00:00:00

set -e ### stops bash script if line ends with error

echo ${HOSTNAME} ${SLURM_ARRAY_TASK_ID}

ml purge
ml load GCCcore/8.3.0 \
        Trim_Galore/0.6.5-Java-11.0.2-Python-3.7.4

nameArray=(EGAR00001508614_SARC061
           EGAR00001508618_SARC065
           EGAR00001508623_SARC070-Primary
           EGAR00001508624_SARC070-Relapse1
           EGAR00001508656_SARC102)

baseName=${nameArray[${SLURM_ARRAY_TASK_ID}]}

inputPath=/gpfs0/home1/gdkendalllab/lab/raw_data/fastq/2016_12_14

trim_galore \
    --length 30 \
    -j 8 \
    --paired \
    -o /gpfs0/scratch/mvc002/kendall \
    ${inputPath}/${baseName}.R1.fastq.gz \
    ${inputPath}/${baseName}.R2.fastq.gz

ml purge
ml load GCC/7.3.0-2.30 \
        OpenMPI/3.1.1 \
        SAMtools/1.9

STAR \
    --runMode alignReads \
    --outSAMtype BAM SortedByCoordinate \
    --runThreadN 10 \
    --outFilterMultimapNmax 1 \
    --readFilesCommand zcat \
    --genomeDir ref/starhg38.p4 \
    --readFilesIn /gpfs0/scratch/mvc002/kendall/${baseName}.R1_val_1.fq.gz \
                  /gpfs0/scratch/mvc002/kendall/${baseName}.R2_val_2.fq.gz \
    --outFileNamePrefix output/aligned/Hs/${baseName}

samtools index output/aligned/Hs/${baseName}Aligned.sortedByCoord.out.bam

rm /gpfs0/scratch/mvc002/kendall/${baseName}.R1_val_1.fq.gz \
   /gpfs0/scratch/mvc002/kendall/${baseName}.R2_val_2.fq.gz \
   /gpfs0/scratch/mvc002/kendall/${baseName}.R1.fastq.gz_trimming_report.txt \
   /gpfs0/scratch/mvc002/kendall/${baseName}.R2.fastq.gz_trimming_report.txt

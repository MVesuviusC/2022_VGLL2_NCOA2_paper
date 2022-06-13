#!/bin/sh
#SBATCH --account=gdkendalllab
#SBATCH --array=0-36
#SBATCH --error=slurmOut/alignFusion-%j.txt
#SBATCH --output=slurmOut/alignFusion-%j.txt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --job-name alignFusion
#SBATCH --wait
#SBATCH --mail-user=matthew.cannon@nationwidechildrens.org
#SBATCH --mail-type=FAIL,REQUEUE,TIME_LIMIT_80

set -e ### stops bash script if line ends with error

echo ${HOSTNAME} ${SLURM_ARRAY_TASK_ID}

module purge

ml load SAMtools/1.15

fileArray=(/home/gdkendalllab/lab/raw_data/fastq/2016_12_14/*.R1.fastq.gz
           /home/gdkendalllab/lab/raw_data/fastq/2015_12_22/*.R1.fastq.gz
           /gpfs0/home1/gdkendalllab/lab/raw_data/fastq/2016_01_01/*R1.fastq.gz)

fileName=${fileArray[${SLURM_ARRAY_TASK_ID}]}
inputPath=${fileName%/*}
baseName=${fileName##*/}
baseName=${baseName%.R1.fastq.gz}

STAR \
    --runMode alignReads \
    --runThreadN 3 \
    --readFilesCommand zcat \
    --genomeDir ref/fusionRef \
    --outTmpDir ${TMPDIR}/temp \
    --outSAMtype BAM Unsorted \
    --readFilesIn ${inputPath}/${baseName}.R1.fastq.gz \
                  ${inputPath}/${baseName}.R2.fastq.gz \
    --outFileNamePrefix output/fusionAligned/${baseName}

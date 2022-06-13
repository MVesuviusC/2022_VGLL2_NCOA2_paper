#!/bin/sh
#SBATCH --account=gdkendalllab
#SBATCH --error=slurmOut/trim-%j.txt
#SBATCH --output=slurmOut/trim-%j.txt
#SBATCH --mem=2G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --job-name trim
#SBATCH --wait
#SBATCH --array=0-7
#SBATCH --time=1-00:00:00
#SBATCH --mail-user=matthew.cannon@nationwidechildrens.org
#SBATCH --mail-type=FAIL,REQUEUE,TIME_LIMIT_80

set -e ### stops bash script if line ends with error

echo ${HOSTNAME} ${SLURM_ARRAY_TASK_ID}

module purge

ml load GCCcore/8.3.0 \
        Trim_Galore/0.6.5-Java-11.0.2-Python-3.7.4

inputPath=/home/gdkendalllab/lab/raw_data/fastq/2015_12_22

fileArray=(${inputPath}/*R1.fastq.gz)

baseName=${fileArray[${SLURM_ARRAY_TASK_ID}]%R1.fastq.gz}

trim_galore \
    --length 50 \
    -o output/trimmed \
    -j 10 \
    --paired \
    ${baseName}R1.fastq.gz \
    ${baseName}R2.fastq.gz

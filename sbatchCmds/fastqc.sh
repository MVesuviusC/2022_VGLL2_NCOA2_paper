#!/bin/sh
#SBATCH --account=gdkendalllab
#SBATCH --error=slurmOut/fastqc-%j.txt
#SBATCH --output=slurmOut/fastqc-%j.txt
#SBATCH --mem=2G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --job-name fastqc
#SBATCH --wait
#SBATCH --array=0-15
#SBATCH --time=1-00:00:00
#SBATCH --mail-user=matthew.cannon@nationwidechildrens.org
#SBATCH --mail-type=FAIL,REQUEUE,TIME_LIMIT_80

set -e ### stops bash script if line ends with error

echo ${HOSTNAME} ${SLURM_ARRAY_TASK_ID}

module purge

module load FastQC/0.11.9-Java-11.0.2

inputPath=/home/gdkendalllab/lab/raw_data/fastq/2015_12_22

fileArray=(${inputPath}/*fastq.gz)

fastqc \
  -o output/fastqc \
  -t 5 \
  --extract \
  ${fileArray[${SLURM_ARRAY_TASK_ID}]}


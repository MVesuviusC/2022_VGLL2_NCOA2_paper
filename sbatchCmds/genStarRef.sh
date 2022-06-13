#!/bin/sh
#SBATCH --account=gdkendalllab
#SBATCH --error=slurmOut/genRef-%j.txt
#SBATCH --output=slurmOut/genRef-%j.txt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --job-name genRef
#SBATCH --wait
#SBATCH --mail-user=matthew.cannon@nationwidechildrens.org
#SBATCH --mail-type=FAIL,REQUEUE,TIME_LIMIT_80

set -e ### stops bash script if line ends with error

echo ${HOSTNAME} ${SLURM_ARRAY_TASK_ID}

module purge

STAR \
    --runMode genomeGenerate \
    --runThreadN 20 \
    --genomeDir ref/starmm10 \
    --genomeFastaFiles ref/starmm10/mm10.fasta \
    --sjdbGTFfile /reference/mus_musculus/mm10/ucsc_assmebly/illumina_download/Annotation/Genes/genes.gtf

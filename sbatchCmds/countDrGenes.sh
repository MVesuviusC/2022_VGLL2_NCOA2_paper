#!/bin/sh
#SBATCH --account=gdkendalllab
#SBATCH --error=slurmOut/countGenesDr-%j.txt
#SBATCH --output=slurmOut/countGenesDr-%j.txt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --job-name countGenesDr
#SBATCH --wait
#SBATCH --array=0-42
#SBATCH --mail-user=matthew.cannon@nationwidechildrens.org
#SBATCH --mail-type=FAIL,REQUEUE,TIME_LIMIT_80
#SBATCH --time=1-24:00:00

set -e ### stops bash script if line ends with error

echo ${HOSTNAME} ${SLURM_ARRAY_TASK_ID}

fileArray=(output/aligned/drKras/*.bam
           output/aligned/drVGLL/*.bam)

inFile=${fileArray[${SLURM_ARRAY_TASK_ID}]}
baseName=${inFile##*/}
baseName=${baseName%Aligned.sortedByCoord.out.bam}

ml purge
ml load HTSeq/0.12.4

htseq-count \
    -n 4 \
    -t exon \
    -c output/geneCounts/geneCountsDr_${baseName}.txt \
    -r pos \
    -s no \
    ${inFile} \
    /home/gdkendalllab/lab/references/annotations/danRer11.ensGene.gtf.gz

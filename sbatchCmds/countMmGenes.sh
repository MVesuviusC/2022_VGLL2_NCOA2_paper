#!/bin/sh
#SBATCH --account=gdkendalllab
#SBATCH --error=slurmOut/countGenesMm-%j.txt
#SBATCH --output=slurmOut/countGenesMm-%j.txt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --job-name countGenesMm
#SBATCH --wait
#SBATCH --array=0-27
#SBATCH --mail-user=matthew.cannon@nationwidechildrens.org
#SBATCH --mail-type=FAIL,REQUEUE,TIME_LIMIT_80
#SBATCH --time=1-24:00:00

set -e ### stops bash script if line ends with error

echo ${HOSTNAME} ${SLURM_ARRAY_TASK_ID}

fileArray=(output/aligned/mmVGLL/*.bam
           output/aligned/mmSkMu/*Aligned.sortedByCoord.out.bam)

inFile=${fileArray[${SLURM_ARRAY_TASK_ID}]}
baseName=${inFile##*/}
baseName=${baseName%Aligned.sortedByCoord.out.bam}

ml purge
ml load HTSeq/0.12.4

htseq-count \
    -n 4 \
    -t exon \
    -c output/geneCounts/geneCountsMm_${baseName}.txt \
    -r pos \
    -s no \
    ${inFile} \
    /reference/mus_musculus/mm10/ucsc_assmebly/illumina_download/Annotation/Genes/genes.gtf

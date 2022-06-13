#!/bin/sh
#SBATCH --account=gdkendalllab
#SBATCH --array=0-7
#SBATCH --error=slurmOut/alignMmVGLL-%j.txt
#SBATCH --output=slurmOut/alignMmVGLL-%j.txt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --job-name alignMmVGLL
#SBATCH --wait
#SBATCH --mail-user=matthew.cannon@nationwidechildrens.org
#SBATCH --mail-type=FAIL,REQUEUE,TIME_LIMIT_80
#SBATCH --time=2-00:00:00

set -e ### stops bash script if line ends with error

echo ${HOSTNAME} ${SLURM_ARRAY_TASK_ID}

ml purge

ml load GCCcore/8.3.0 \
        Trim_Galore/0.6.5-Java-11.0.2-Python-3.7.4

inputPath=/home/gdkendalllab/lab/raw_data/fastq/2015_12_22
tempPath=/gpfs0/scratch/mvc002/kendall

fileArray=(${inputPath}/*R1.fastq.gz)

R1=${fileArray[${SLURM_ARRAY_TASK_ID}]}
R2=${R1/R1.fastq/R2.fastq}

baseName=${R1%.R1.fastq.gz}
baseName=${baseName##*/}

trim_galore \
    --length 30 \
    -o ${tempPath} \
    -j 10 \
    --paired \
    ${R1} \
    ${R2}

# Hard trim to 50bp to keep consistent across samples
trimLen=150
perl ~/oldCode/oldScripts/trimFastq.pl \
    --max ${trimLen} \
    --fastq ${tempPath}/${baseName}.R1_val_1.fq.gz |
    pigz \
        > ${tempPath}/${baseName}_1.fastq.gz

perl ~/oldCode/oldScripts/trimFastq.pl \
    --max ${trimLen} \
    --fastq ${tempPath}/${baseName}.R2_val_2.fq.gz |
    pigz \
        > ${tempPath}/${baseName}_2.fastq.gz

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
    --genomeDir ref/starmm10 \
    --readFilesIn ${tempPath}/${baseName}_1.fastq.gz \
                  ${tempPath}/${baseName}_2.fastq.gz \
    --outFileNamePrefix output/aligned/mmVGLL/${baseName##*/}

samtools index output/aligned/mmVGLL/${baseName##*/}Aligned.sortedByCoord.out.bam

rm ${tempPath}/${baseName}.R[12]_val_[12].fq.gz \
   ${tempPath}/${baseName}.R[12].fastq.gz_trimming_report.txt \
   ${tempPath}/${baseName}_[12].fastq.gz

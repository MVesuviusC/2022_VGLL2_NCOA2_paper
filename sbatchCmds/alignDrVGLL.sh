#!/bin/sh
#SBATCH --account=gdkendalllab
#SBATCH --array=0-36
#SBATCH --error=slurmOut/alignDrVGLL-%j.txt
#SBATCH --output=slurmOut/alignDrVGLL-%j.txt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --job-name alignDrVGLL
#SBATCH --wait
#SBATCH --mail-user=matthew.cannon@nationwidechildrens.org
#SBATCH --mail-type=FAIL,REQUEUE,TIME_LIMIT_80
#SBATCH --time=3-00:00:00

set -e ### stops bash script if line ends with error

echo ${HOSTNAME} ${SLURM_ARRAY_TASK_ID}

module purge
ml load GCCcore/8.3.0 \
        Trim_Galore/0.6.5-Java-11.0.2-Python-3.7.4

R1Array=(/home/gdkendalllab/lab/raw_data/fastq/2016_01_02/*R1*gz
         /home/gdkendalllab/lab/raw_data/fastq/2016_01_01/*R1*gz)

R1=${R1Array[${SLURM_ARRAY_TASK_ID}]}
R2=${R1/R1.fastq/R2.fastq}

baseName=${R1##*/}
baseName=${baseName%.R1.fastq.gz}

trim_galore \
    --length 30 \
    -j 8 \
    --paired \
    ${R1} \
    ${R2}

# Hard trim to 50bp to keep consistent across samples
trimLen=50
perl ~/oldCode/oldScripts/trimFastq.pl \
    --max ${trimLen} \
    --fastq ${baseName}.R1_val_1.fq.gz |
    pigz \
        > ${baseName}_1.fastq.gz

perl ~/oldCode/oldScripts/trimFastq.pl \
    --max ${trimLen} \
    --fastq ${baseName}.R2_val_2.fq.gz |
    pigz \
        > ${baseName}_2.fastq.gz

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
    --genomeDir /gpfs0/home1/gdkendalllab/lab/references/star/danRer11 \
    --readFilesIn ${baseName}_1.fastq.gz \
                  ${baseName}_2.fastq.gz \
    --outFileNamePrefix output/aligned/drVGLL/${baseName}

samtools index output/aligned/drVGLL/${baseName}Aligned.sortedByCoord.out.bam

rm ${baseName##*/}.R[12]_val_[12].fq.gz \
   ${baseName##*/}.R[12].fastq.gz_trimming_report.txt \
   ${baseName}_[12].fastq.gz

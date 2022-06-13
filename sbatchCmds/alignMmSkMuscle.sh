#!/bin/sh
#SBATCH --account=gdkendalllab
#SBATCH --array=0-19%10
#SBATCH --error=slurmOut/alignMmSkMu-%j.txt
#SBATCH --output=slurmOut/alignMmSkMu-%j.txt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --job-name alignMmSkMu
#SBATCH --wait
#SBATCH --mail-user=matthew.cannon@nationwidechildrens.org
#SBATCH --mail-type=FAIL,REQUEUE,TIME_LIMIT_80

set -e ### stops bash script if line ends with error

echo ${HOSTNAME} ${SLURM_ARRAY_TASK_ID}

module purge

ml load GCCcore/8.3.0 \
        Trim_Galore/0.6.5-Java-11.0.2-Python-3.7.4

curPath=$(pwd)

sraArray=($(cut -f 2 misc/mouse_skeletal_muscle_data.txt | grep -v SRA))

sraNum=${sraArray[${SLURM_ARRAY_TASK_ID}]}

cd /gpfs0/scratch/mvc002/kendall/

prefetch --max-size 100G ${sraNum}

fasterq-dump -S ${sraNum}

rm -r ${sraNum}/

trim_galore \
    --length 30 \
    -j 8 \
    --paired \
    ${sraNum}_1.fastq \
    ${sraNum}_2.fastq

trimLen=150
perl ~/oldCode/oldScripts/trimFastq.pl \
    --max ${trimLen} \
    --fastq ${sraNum}_1_val_1.fq |
    pigz \
        > ${sraNum}_1.fastq.gz

perl ~/oldCode/oldScripts/trimFastq.pl \
    --max ${trimLen} \
    --fastq ${sraNum}_2_val_2.fq |
    pigz \
        > ${sraNum}_2.fastq.gz

rm \
    ${sraNum}_1.fastq \
    ${sraNum}_2.fastq \
    ${sraNum}_1_val_1.fq \
    ${sraNum}_2_val_2.fq \
    ${sraNum}_1.fastq_trimming_report.txt \
    ${sraNum}_2.fastq_trimming_report.txt

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
    --genomeDir ${curPath}/ref/starmm10 \
    --readFilesIn ${sraNum}_1.fastq.gz \
                  ${sraNum}_2.fastq.gz \
    --outFileNamePrefix ${curPath}/output/aligned/mmSkMu/${sraNum}

samtools index ${curPath}/output/aligned/mmSkMu/${sraNum}Aligned.sortedByCoord.out.bam

rm ${sraNum}_1.fastq.gz ${sraNum}_2.fastq.gz

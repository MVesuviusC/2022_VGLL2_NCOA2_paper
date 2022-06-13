#!/bin/sh
#SBATCH --account=gdkendalllab
#SBATCH --array=0-5
#SBATCH --error=slurmOut/alignKras-%j.txt
#SBATCH --output=slurmOut/alignKras-%j.txt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --job-name alignKras
#SBATCH --wait
#SBATCH --mail-user=matthew.cannon@nationwidechildrens.org
#SBATCH --mail-type=FAIL,REQUEUE,TIME_LIMIT_80

set -e ### stops bash script if line ends with error

echo ${HOSTNAME} ${SLURM_ARRAY_TASK_ID}

module purge

ml load GCCcore/8.3.0 \
        Trim_Galore/0.6.5-Java-11.0.2-Python-3.7.4

curPath=$(pwd)

sraArray=(SRR6507311
          SRR6507312
          SRR6507313
          SRR6507300
          SRR6507301
          SRR6507302)

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

# Hard trim to 50bp to keep consistent across samples
trimLen=50
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
    --readFilesIn ${sraNum}_1.fastq.gz \
                  ${sraNum}_2.fastq.gz \
    --outFileNamePrefix ${curPath}/output/aligned/drKras/${sraNum}

samtools index ${curPath}/output/aligned/drKras/${sraNum}Aligned.sortedByCoord.out.bam

rm ${sraNum}_[12]_val_[12].fq \
    ${sraNum}_[12].fastq \
    ${sraNum}_[12].fastq_trimming_report.txt \
    ${sraNum}_[12].fastq.gz

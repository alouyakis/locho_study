#!/usr/bin/env bash
#SBATCH --cpus-per-task=10
#SBATCH --mem=10G
#SBATCH --job-name=trim
#SBATCH --output=scripts/logs/trim_%A_%a.out
#SBATCH --array=1-106%5
#SBATCH --export=ALL
#SBATCH -p defq
#SBATCH -N 1
hostname; pwd; date

source /home/artemisl/miniforge3/etc/profile.d/conda.sh
conda activate mtenv

RUN=${SLURM_ARRAY_TASK_ID}

export rawseqs_mg="/data/Microbiome/CosmosID/TOBI_Study"
export filtseqs="filtered"
export nohost="nohost"
export dog="/data/databases/kneaddata_2023/dog/canFam4"

INPUT_PATH=$(ls ${rawseqs_mg}/*R1*.fastq.gz | sed -n "${RUN} p")
sample=$(basename -a ${INPUT_PATH} | cut -f 1 -d '_')

trimmomatic PE \
  ${rawseqs_mg}/${sample}_*R1*.fastq.gz ${rawseqs_mg}/${sample}_*R2*.fastq.gz \
  ${filtseqs}/mg_${sample}_trimmed-pair_R1.fastq.gz ${filtseqs}/mg_${sample}_trimmed-single_R1.fastq.gz \
  ${filtseqs}/mg_${sample}_trimmed-pair_R2.fastq.gz ${filtseqs}/mg_${sample}_trimmed-single_R2.fastq.gz \
  ILLUMINACLIP:/data/databases/adapters/TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36

bowtie2 -p 10 -x ${dog} \
  -1 ${filtseqs}/mg_${sample}_trimmed-pair_R1.fastq.gz \
  -2 ${filtseqs}/mg_${sample}_trimmed-pair_R2.fastq.gz \
  -S ${nohost}/mg_${sample}_mapped_and_unmapped_dog.sam;
samtools view -@ 10 -bS ${nohost}/mg_${sample}_mapped_and_unmapped_dog.sam > ${nohost}/mg_${sample}_mapped_and_unmapped_dog.bam;
samtools view -@ 10 -b -f 12 -F 256 ${nohost}/mg_${sample}_mapped_and_unmapped_dog.bam > ${nohost}/mg_${sample}_unmapped_dog.bam;
samtools sort -n ${nohost}/mg_${sample}_unmapped_dog.bam > ${nohost}/mg_${sample}_unmapped_sorted_dog.bam;
bedtools bamtofastq -i ${nohost}/mg_${sample}_unmapped_sorted_dog.bam \
  -fq ${nohost}/mg_${sample}_rmhost_dog_r1.fq \
  -fq2 ${nohost}/mg_${sample}_rmhost_dog_r2.fq

rm ${nohost}/mg_${sample}_mapped_and_unmapped_dog.sam
rm ${nohost}/mg_${sample}_mapped_and_unmapped_dog.bam
rm ${nohost}/mg_${sample}_unmapped_dog.bam
rm ${nohost}/mg_${sample}_unmapped_sorted_dog.bam

done


#!/usr/bin/env bash
#SBATCH --cpus-per-task=10
#SBATCH --mem=50G
#SBATCH --job-name=assembly
#SBATCH --output=scripts/logs/spades_%A_%a.out
#SBATCH --array=1-106%5
#SBATCH --export=ALL
#SBATCH -p defq
#SBATCH -N 1
hostname; pwd; date

source /home/artemisl/miniforge3/etc/profile.d/conda.sh
conda activate mtenv

export rawseqs_mg="/data/Microbiome/CosmosID/TOBI_Study"
export nohost="nohost"
export spades="spades"

RUN=${SLURM_ARRAY_TASK_ID}

INPUT_PATH=$(ls ${rawseqs_mg}/*R1*.fastq.gz | sed -n "${RUN} p")
sample=$(basename -a ${INPUT_PATH} | cut -f 1 -d '_')

echo -e "\nRun ID: ${RUN}"
echo -e "\nSample: ${sample}"

spades.py --meta -t 10 -m 50 \
  -1 ${nohost}/mg_${sample}*_rmhost_dog_r1.fq \
  -2 ${nohost}/mg_${sample}*_rmhost_dog_r2.fq \
  -o ${spades}/nohost_${sample}

date


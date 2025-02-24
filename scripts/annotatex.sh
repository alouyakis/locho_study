#!/bin/sh
#SBATCH --job-name=annotatex
#SBATCH --qos=defq
#SBATCH --cpus-per-task=25
#SBATCH --mem=125gb
#SBATCH -t 21-00:00:00
#SBATCH -o scripts/annotatex_%j.out
pwd; hostname; date

conda activate mtenv

export trinity_out="/home/artemisl/locho_study/trinity_out"
export annotations="/home/artemisl/locho_study/annotations"

export uniprot="/data/databases/uniprot"

## run simulatanous to transdecoder
blastx -query ${trinity_out}/Trinity.fasta \
  -db ${uniprot}/uniprot_sprot.fasta \
  -out ${annotations}/blastx.outfmt6 \
  -evalue 1e-3 \
  -num_threads 25 \
  -outfmt 6

date

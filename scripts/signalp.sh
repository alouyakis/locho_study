#!/bin/sh
#SBATCH --job-name=signalp
#SBATCH --qos=defq
#SBATCH --cpus-per-task=10
#SBATCH --mem=20gb
#SBATCH -t 21-00:00:00
#SBATCH -o scripts/signalp_%j.out
pwd; hostname; date

conda activate mtenv

export annotations="/home/artemisl/locho_study/annotations"
export transdecoder="${annotations}/transdecoder"

signalp \
  -format short \
  -plot png \
  -org euk \
  -tmp tmp \
  -prefix ${annotations}/signalp_euk.out \
  -fasta ${transdecoder}/Trinity.fasta.transdecoder.pep

signalp \
  -format short \
  -plot png \
  -org gram- \
  -tmp tmp \
  -prefix ${annotations}/signalp_gramn.out \
  -fasta ${transdecoder}/Trinity.fasta.transdecoder.pep

signalp \
  -format short \
  -plot png \
  -org gram+ \
  -tmp tmp \
  -prefix ${annotations}/signalp_gramp.out \
  -fasta ${transdecoder}/Trinity.fasta.transdecoder.pep

date

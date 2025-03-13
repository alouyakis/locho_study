#!/bin/sh
#SBATCH --job-name=trinotate
#SBATCH --qos=defq
#SBATCH --cpus-per-task=10
#SBATCH --mem=20gb
#SBATCH -t 21-00:00:00
#SBATCH -o scripts/trinotate_%j.out
pwd; hostname; date

conda activate mtenv

export trinity_out="/home/artemisl/locho_study/trinity_out"
export annotations="/home/artemisl/locho_study/annotations"
export transdecoder="${annotations}/transdecoder"

Trinotate ${annotations}/Trinotate.sqlite init \
  --gene_trans_map ${trinity_out}/Trinity.fasta.gene_trans_map \
  --transcript_fasta ${trinity_out}/Trinity.fasta \
  --transdecoder_pep ${transdecoder}/Trinity.fasta.transdecoder.pep

Trinotate ${annotations}/Trinotate.sqlite LOAD_swissprot_blastx ${annotations}/filter_blastx.fmt6
Trinotate ${annotations}/Trinotate.sqlite LOAD_swissprot_blastp ${annotations}/filter_blastp.fmt6
Trinotate ${annotations}/Trinotate.sqlite LOAD_pfam ${annotations}/TrinotatePFAM.out
Trinotate ${annotations}/Trinotate.sqlite LOAD_tmhmm ${annotations}/tmhmm.out

Trinotate ${annotations}/Trinotate.sqlite report > ${annotations}/trinotate_annotation_report.xls

date

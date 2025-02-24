#!/bin/sh
#SBATCH --job-name=annotatep
#SBATCH --qos=defq
#SBATCH --cpus-per-task=25
#SBATCH --mem=125gb
#SBATCH -t 21-00:00:00
#SBATCH -o scripts/annotatep_%j.out
pwd; hostname; date

conda activate mtenv

export trinity_out="/home/artemisl/locho_study/trinity_out"
export annotations="/home/artemisl/locho_study/annotations"
export transdecoder="${annotations}/transdecoder"

export uniprot="/data/databases/uniprot"
export pfam="/data/databases/pfam"

## identify protein coding regions
TransDecoder.LongOrfs -t ${trinity_out}/Trinity.fasta \
  --output_dir ${transdecoder}
TransDecoder.Predict -t ${trinity_out}/Trinity.fasta \
  --output_dir ${transdecoder}
mv Trinity.fasta.transdecoder.* ${transdecoder}

## run blastp on protein coding regions
blastp -query ${transdecoder}/Trinity.fasta.transdecoder.pep \
  -db ${uniprot}/uniprot_sprot.fasta \
  -out ${annotations}/blastp.outfmt6 \
  -evalue 1e-3 \
  -num_threads 25 \
  -outfmt 6

## filter blast output for top hit
mkdir tmp
sort -T tmp -g -k11,11 ${annotations}/blastp.outfmt6 | sort -T tmp -u -k1,1 > ${annotations}/filter_blastp.fmt6
rm -r tmp/

## identify protein families
hmmscan \
  --cpu 25 \
  --domtblout ${annotations}/TrinotatePFAM.out \
  ${pfam}/Pfam-A.hmm \
  ${transdecoder}/Trinity.fasta.transdecoder.pep

## identify transmembrane proteins
tmhmm --short \
  ${transdecoder}/Trinity.fasta.transdecoder.pep > \
  ${annotations}/tmhmm.out

date

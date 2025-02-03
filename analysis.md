**project: locho metatranscriptome**
**date: 2025.01.31**
**author:**
- name: Artemis S. Louyakis

Sample number: 106
sequence path: /data/Microbiome/NextSeq2000/250115_VH00394_26_2225LTJNX/Analysis/1/Data/fastq
analysis path: /home/artemisl/locho_study

set up:
```bash
## clone git repo
git clone git@github.com:alouyakis/locho_study.git

## activate environment
screen -R locho_study
conda activate mtenv

## create output dirs and variables
rawseqs="/data/Microbiome/NextSeq2000/250115_VH00394_26_2225LTJNX/Analysis/1/Data/fastq"
mkdir -p quality && export quality="quality"
mkdir -p filtered && export filtseqs="filtered"
mkdir -p sortmerna && export sortmerna="sortmerna"
mkdir -p nohost && export nohost="nohost"
mkdir -p tables && export tables="tables"
mkdir -p alignments && export alignments="alignments"
mkdir -p annotations && export annotations="annotations"
mkdir -p ${annotations}/transdecoder && export transdecoder="${annotations}/transdecoder"
export trinity_out="trinity_out"

export rrna_db="/data/databases/sortmerna/data/rRNA_databases"
export dog="/data/databases/kneaddata_2023/dog/canFam4"
export cat="/data/databases/kneaddata_2023/cat/fca126"
export uniprot="/data/databases/uniprot"
export pfam="/data/databases/pfam"

## ignore undetermined - may need sudo
sudo mv ${rawseqs}/Undetermined_S0_R1_001.fastq.gz ${rawseqs}/Undetermined_S0_R1_001.fastq.gz.ignore
sudo mv ${rawseqs}/Undetermined_S0_R2_001.fastq.gz ${rawseqs}/Undetermined_S0_R2_001.fastq.gz.ignore
```

quality trim:
```bash
parallel -j 10 \
  trimmomatic PE \
    ${rawseqs}/{1}_*R1*.fastq.gz ${rawseqs}/{1}_*R2*.fastq.gz \
    ${filtseqs}/{1}_trimmed-pair_R1.fastq.gz ${filtseqs}/{1}_trimmed-single_R1.fastq.gz \
    ${filtseqs}/{1}_trimmed-pair_R2.fastq.gz ${filtseqs}/{1}_trimmed-single_R2.fastq.gz \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36 ::: $( basename -a ${rawseqs}/*R1*.fastq.gz | cut -f 1 -d '_')
```

remove rRNA from trimmed reads:
```bash
parallel -j 10 \
  sortmerna --ref ${rrna_db}/rfam-5.8s-database-id98.fasta \
    --ref ${rrna_db}/rfam-5s-database-id98.fasta \
    --ref ${rrna_db}/silva-arc-16s-id95.fasta \
    --ref ${rrna_db}/silva-arc-23s-id98.fasta \
    --ref ${rrna_db}/silva-bac-16s-id90.fasta \
    --ref ${rrna_db}/silva-bac-23s-id98.fasta \
    --ref ${rrna_db}/silva-euk-18s-id95.fasta \
    --ref ${rrna_db}/silva-euk-28s-id98.fasta \
    --reads ${filtseqs}/{1}_trimmed-pair_R1.fastq.gz \
    --reads ${filtseqs}/{1}_trimmed-pair_R2.fastq.gz \
    --aligned ${sortmerna}/rrna_{1} \
    --other ${sortmerna}/mrna_{1} --out2 \
    --paired_in --fastx --num_alignments 1 \
    --workdir ${sortmerna}/{1}_wd \
    --threads 5 -m 5000 ::: $( basename -a ${rawseqs}/*R1*.fastq.gz | cut -f 1 -d '_')
```

sequence counts:
```bash
## raw counts
## QC - counts should match:
zgrep -c "^+$" ${rawseqs}/*R1*.fastq.gz > ${tables}/raw_fwd.csv
zgrep -c "^+$" ${rawseqs}/*R2*.fastq.gz > ${tables}/raw_rev.csv

## filtered counts
## QC - counts should match:
zgrep -c "^+$" ${filtseqs}/*pair_R1*.fastq.gz > ${tables}/filt_paired_fwd.csv
zgrep -c "^+$" ${filtseqs}/*pair_R2*.fastq.gz > ${tables}/filt_paired_rev.csv
## singleton counts:
zgrep -c "^+$" ${filtseqs}/*single_R1*.fastq.gz > ${tables}/filt_singles_fwd.csv
zgrep -c "^+$" ${filtseqs}/*single_R2*.fastq.gz > ${tables}/filt_singles_rev.csv

## rrna vs mrna counts


## mrna with host removed counts

```








slurm script template:
```bash
#!/bin/sh
#SBATCH --job-name=align
#SBATCH --cpus-per-task=10
#SBATCH --mem=50gb
#SBATCH -t 100:00:00
#SBATCH -o sortme_%A_%a.out
#SBATCH --array=1-111%10
pwd; hostname; date

conda activate mtenv

RUN=${SLURM_ARRAY_TASK_ID}
nohost="rm_host_mt"
sortmerna="sortmerna_nohost"
mkdir -p ${sortmerna}
db="/data/databases/sortmerna/data/rRNA_databases"

INPUT_PATH=$(ls ${nohost}/*_rmhost_r1.fqq | sed -n "${RUN} p" | sed 's/_rmhost_r1.fq//' )
SAMPLE=$(basename ${INPUT_PATH})

echo -e "\nRun ID: ${RUN}"
echo -e "\nSample: ${SAMPLE}"
```

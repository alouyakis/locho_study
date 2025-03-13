**project: locho mg and mt analysis**
**date: 2025.03.13**
**author:**
- name: Artemis S. Louyakis

rerunning analysis using metagenome data. for mt steps, will use processed reads from analysis.md. combine pipelines after comparisons.
- compare annotations from assembled mg to biobakery ref based annotations
- compare de novo mt assembly annotations to mg assist

Sample number: 106
sequence path: /data/Microbiome/NextSeq2000/250115_VH00394_26_2225LTJNX/Analysis/1/Data/fastq
analysis path: /home/artemisl/locho_study
metagenome path: /data/Microbiome/CosmosID/TOBI_Study
mg google drive: https://drive.google.com/drive/folders/1DOPGLeTZeCOzhImquL7inTqiR6rrSfGt

set up:
```bash
## clone git repo
git clone git@github.com:alouyakis/locho_study.git

## activate environment
screen -R locho_study_mg
conda activate mgenv

## create output dirs and variables
export rawseqs_mt="/data/Microbiome/NextSeq2000/250115_VH00394_26_2225LTJNX/Analysis/1/Data/fastq"
export rawseqs_mg="/data/Microbiome/CosmosID/TOBI_Study"
export data="data"
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
    ${rawseqs_mg}/{1}_*R1*.fastq.gz ${rawseqs_mg}/{1}_*R2*.fastq.gz \
    ${filtseqs}/mg_{1}_trimmed-pair_R1.fastq.gz ${filtseqs}/mg_{1}_trimmed-single_R1.fastq.gz \
    ${filtseqs}/mg_{1}_trimmed-pair_R2.fastq.gz ${filtseqs}/mg_{1}_trimmed-single_R2.fastq.gz \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36 ::: $( basename -a ${rawseqs_mg}/*R1*.fastq.gz | cut -f 1 -d '_')
```

remove host from mRNA:  
```bash
## align to dog and output fq for unaligned (i.e. non-host)
parallel -j 10 \
  'bowtie2 -p 5 -x ${dog} \
    -1 ${sortmerna}/mrna_{1}_fwd.fq.gz \
    -2 ${sortmerna}/mrna_{1}_rev.fq.gz \
    -S ${nohost}/mrna_{1}_mapped_and_unmapped_dog.sam;
  samtools view -@ 5 -bS ${nohost}/mrna_{1}_mapped_and_unmapped_dog.sam > ${nohost}/mrna_{1}_mapped_and_unmapped_dog.bam;
  samtools view -@ 5 -b -f 12 -F 256 ${nohost}/mrna_{1}_mapped_and_unmapped_dog.bam > ${nohost}/mrna_{1}_unmapped_dog.bam;
  samtools sort -n ${nohost}/mrna_{1}_unmapped_dog.bam > ${nohost}/mrna_{1}_unmapped_sorted_dog.bam;
  bedtools bamtofastq -i ${nohost}/mrna_{1}_unmapped_sorted_dog.bam \
    -fq ${nohost}/mrna_{1}_rmhost_dog_r1.fq \
    -fq2 ${nohost}/mrna_{1}_rmhost_dog_r2.fq' ::: $( basename -a ${rawseqs}/*R1*.fastq.gz | cut -f 1 -d '_')

## output reads aligned to host - not run - TODO: add to unaligned and remove sam/bam when finished
# parallel -j 5 \
#   'gunzip ${nohost}/mrna_{1}_mapped_and_unmapped_dog.bam.gz;
#   samtools view -@ 5 -b -f 12 -F 4 ${nohost}/mrna_{1}_mapped_and_unmapped_dog.bam > ${nohost}/mrna_{1}_mapped_dog.bam;
#   samtools sort -n ${nohost}/mrna_{1}_mapped_dog.bam > ${nohost}/mrna_{1}_mapped_sorted_dog.bam;
#   bedtools bamtofastq -i ${nohost}/mrna_{1}_mapped_sorted_dog.bam \
#     -fq ${nohost}/mrna_{1}_dog_r1.fq \
#     -fq2 ${nohost}/mrna_{1}_dog_r2.fq' ::: $( basename -a ${rawseqs}/*R1*.fastq.gz | cut -f 1 -d '_')

## cleanup - may just delete
gzip ${nohost}/*.sam
gzip ${nohost}/*.bam
```












create yaml for spades if co-assembling
```bash
cd ${rawseqs}
printf "[\n  {\n    orientation: \"fr\",\n    type: \"paired-end\",\n    right reads: [\n" > ${wd}/input.yaml
for i in *_R1_*; do
  r1=$(realpath ${i})
  printf "      \"${r1}\"\n" >> ${wd}/input.yaml.tmp;
done
sed -i '$!s/$/,/' ${wd}/input.yaml.tmp
cat ${wd}/input.yaml.tmp >> ${wd}/input.yaml && rm ${wd}/input.yaml.tmp
printf "    ],\n    left reads: [\n" >> ${wd}/input.yaml
for i in *_R2_*; do
  r2=$(realpath ${i})
  printf "      \"${r2}\"\n" >> ${wd}/input.yaml.tmp;
done
sed -i '$!s/$/,/' ${wd}/input.yaml.tmp
cat ${wd}/input.yaml.tmp >> ${wd}/input.yaml && rm ${wd}/input.yaml.tmp
printf "    ]\n  }\n]" >> ${wd}/input.yaml
cd ${wd}

```

scratch slurm script construction
TODO: write pipeline to create slurm scripts
```
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



assembly="/home/artemisl/haircare/phase-ii/assembly"

INPUT_PATH=$(ls ${data_in}/rmhost_*_R1_001.fastq | sed -n "${RUN} p" | sed 's/_R1_001.fastq//' )
sample=$(basename ${INPUT_PATH})

echo -e "\nRun ID: ${RUN}"
echo -e "\nSample: ${sample}"

spades.py --meta -t 10 -m 50 \
  -1 ${data_in}/${sample}_R1_001.fastq \
  -2 ${data_in}/${sample}_R2_001.fastq \
  -o ${assembly}/nohost_${sample}

date
```

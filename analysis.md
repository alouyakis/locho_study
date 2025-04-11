**project: locho metatranscriptome**
**date: 2025.01.31**
**author:**
- name: Artemis S. Louyakis

Sample number: 106
sequence path: /data/Microbiome/NextSeq2000/250115_VH00394_26_2225LTJNX/Analysis/1/Data/fastq
analysis path: /home/artemisl/locho_study
metagenome path: /data/Microbiome/CosmosID/TOBI_Study
mg google drive: https://drive.google.com/drive/folders/1DOPGLeTZeCOzhImquL7inTqiR6rrSfGt

questions for matt and dayakar about the analysis -
1. use the metagenomes from the same samples to improve annotating taxa and functions?
2. will this study include the metabolite data with the unknowns?

considerations:
 - crossover study, so look into mixed effects model for repeated measures in mg/mt data analysis


set up:
```bash
## clone git repo
git clone git@github.com:alouyakis/locho_study.git

## activate environment
screen -R locho_study
conda activate mtenv

## create output dirs and variables
export rawseqs="/data/Microbiome/NextSeq2000/250115_VH00394_26_2225LTJNX/Analysis/1/Data/fastq"
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
zgrep -c "^+$" ${sortmerna}/mrna_*fwd*.fq.gz > ${tables}/mrna_paired_fwd.csv
zgrep -c "^+$" ${sortmerna}/mrna_*rev*.fq.gz > ${tables}/mrna_paired_rev.csv
zgrep -c "^+$" ${sortmerna}/rrna_*fwd*.fq.gz > ${tables}/rrna_paired_fwd.csv
zgrep -c "^+$" ${sortmerna}/rrna_*rev*.fq.gz > ${tables}/rrna_paired_rev.csv

## mrna with host removed counts
grep -c "^+$" ${nohost}/mrna_*_rmhost_dog_r1.fq > ${tables}/mrna_nohost_fwd.csv
grep -c "^+$" ${nohost}/mrna_*_rmhost_dog_r2.fq > ${tables}/mrna_nohost_rev.csv
```

quality check:
```bash
## fastqc
parallel -j 10 \
  fastqc ${rawseqs}/{1}*R1*.fastq.gz ${rawseqs}/{1}*R2*.fastq.gz \
    -o ${quality} -t 5 ::: $( basename -a ${rawseqs}/*R1*.fastq.gz | cut -f 1 -d '_')

## fastx toolkit
parallel -j 5 \
  'gunzip -c ${rawseqs}/{1}*R1*.fastq.gz | fastx_quality_stats \
  -o ${quality}/fastx_{1}_r1_stats.csv' ::: $( basename -a ${rawseqs}/*R1*.fastq.gz | cut -f 1 -d '_')
parallel -j 5 \
  gunzip -c ${rawseqs}/{1}*R2*.fastq.gz | fastx_quality_stats \
  -o ${quality}/fastx_{1}_r1_stats.csv ::: $( basename -a ${rawseqs}/*R1*.fastq.gz | cut -f 1 -d '_')

## TODO - try RSeQC https://rseqc.sourceforge.net/
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

### two methods for running trinity - create input file or dynamically create command from input directory
method 1 - generate command with fwd/rev seq lists for trinity inputs:
```bash
## create command variables
cmd1="Trinity \
        --seqType fq \
        --max_memory 300G \
        --normalize_reads \
        --CPU 75 \
        --output ${trinity_out} \
        --left"

cmd2="        --right"

## append with sequence file names
first=true
for file in $(find "${nohost}" -name "mrna_*_rmhost_dog_r1.fq"); do
  fwd="${file}"
  rev="${file/_r1/_r2}"

  if $first; then
    cmd1+=" ${fwd}"
    cmd2+=" ${rev}"
    first=false
  else
    cmd1+=",${fwd}"
    cmd2+=",${rev}"
  fi
done

## combine into single command
cmd="$cmd1 $cmd2"

# execute the command
eval $cmd
```

method 2 - create trinity input file:
this method is less reproducible
```bash
## filter sample data by sample ids
awk -F":" '{print $1}' ${tables}/filt_paired_fwd.csv | sed 's/filtered\///g;s/_trimmed-pair_R1.fastq.gz//g' > ${data}/sampleids.tmp
grep -Ff data/sampleids.tmp ${data}/locho_105536_sample_data.tsv > ${data}/filtered_locho_105536_sample_data.tsv

## create temp files for each column
awk -F"\t" '{OFS="\t"}{print $4,$28,$29}' ${data}/filtered_locho_105536_sample_data.tsv | sort | sed 's/Pre-Feed\t/prefeed_/g;s/Treatment (test)\t/treatment_/g' | awk -F"\t" '{print $2}' > ${data}/col1.tmp
awk -F"\t" '{OFS="\t"}{print $4,$28,$29}' ${data}/filtered_locho_105536_sample_data.tsv | sort | sed 's/Pre-Feed\t/prefeed_/g;s/Treatment (test)\t/treatment_/g' | awk -F"\t" '{split($2, arr, "_"); key = arr[1]"_"arr[2]; count[key]++; print $1"\t"$2"_"count[key]}' | sort | awk -F"\t" '{print $2}' > ${data}/col2.tmp
awk -F":" '{print $1}' ${tables}/mrna_paired_fwd.csv | sort > ${data}/col3.tmp
awk -F":" '{print $1}' ${tables}/mrna_paired_rev.csv | sort > ${data}/col4.tmp
## paste together
paste -d '\t' ${data}/col1.tmp ${data}/col2.tmp ${data}/col3.tmp ${data}/col4.tmp > ${data}/trinity_input.txt

## remove tmp files
rm ${data}/*.tmp

## run trinity
Trinity \
  --seqType fq \
  --max_memory 300G \
  --normalize_reads \
  --CPU 75 \
  --output ${trinity_out} \
  --samples_file ${data}/trinity_input.txt
```

assembly housekeeping:
```
~/.conda/pkgs/trinity-2.8.5-h8b12597_5/bin/get_Trinity_gene_to_trans_map.pl ${trinity_out}/Trinity.fasta > ${trinity_out}/Trinity.fasta.gene_trans_map
TrinityStats.pl ${trinity_out}/Trinity.fasta > ${tables}/TrinityStats_out.txt

for i in ${trinity_out}/read_partitions/Fb_*/CBin_*; do rm -r ${i}; done
```

align reads to assembly:
```bash
bowtie2-build ${trinity_out}/Trinity.fasta ${trinity_out}/Trinity.fasta
parallel -j 10 \
  'bowtie2 -p 10 -q --no-unal -k 20 -x ${trinity_out}/Trinity.fasta \
    -1 ${nohost}/mrna_{1}_rmhost_dog_r1.fq \
    -2 ${nohost}/mrna_{1}_rmhost_dog_r2.fq \
    2>${tables}/align_stats_{1}.txt | samtools view -@10 -Sb \
    -o ${alignments}/bowtie2_{1}.bam' ::: $( basename -a ${rawseqs}/*R1*.fastq.gz | cut -f 1 -d '_')

cat 2>&1 ${tables}/align_stats_*.txt | less

## alternative
parallel -j 10 \
  'bowtie2 -p 10 --local -q --no-unal -x ${trinity_out}/Trinity.fasta \
    -1 ${nohost}/mrna_{1}_rmhost_dog_r1.fq \
    -2 ${nohost}/mrna_{1}_rmhost_dog_r2.fq \
    2>${tables}/alt_align_stats_{1}.txt | samtools view -@10 -Sb - | \
    samtools sort -o ${alignments}/alt_bowtie2_{1}.bam' ::: $( basename -a ${rawseqs}/*R1*.fastq.gz | cut -f 1 -d '_')

bowtie2 -p 10 --local -q --no-unal -x ${trinity_out}/Trinity.fasta \
  -1 ${nohost}/mrna_F-50179058-003S_rmhost_dog_r1.fq \
  -2 ${nohost}/mrna_F-50179058-003S_rmhost_dog_r2.fq \
  2>${tables}/alt_align_stats_F-50179058-003S.txt | samtools view -@10 -Sb - | \
  samtools sort -o ${alignments}/alt_bowtie2_F-50179058-003S.bam
```

alignment counts:
```bash
## prep reference
align_and_estimate_abundance.pl \
  --transcripts ${trinity_out}/Trinity.fasta \
  --est_method RSEM \
  --aln_method bowtie2 \
  --trinity_mode \
  --output_dir ${trinity_out} \
  --prep_reference

parallel -j 5 \
  align_and_estimate_abundance.pl \
    --transcripts ${trinity_out}/Trinity.fasta \
    --seqType fq \
    --left ${nohost}/mrna_{1}_rmhost_dog_r1.fq \
    --right ${nohost}/mrna_{1}_rmhost_dog_r2.fq \
    --est_method RSEM \
    --output_dir ${alignments}/mrna_{1} \
    --aln_method bowtie2 \
    --thread_count 10 \
    --gene_trans_map ${trinity_out}/Trinity.fasta.gene_trans_map \
    --output_prefix mrna_{1} ::: $( basename -a ${rawseqs}/*R1*.fastq.gz | cut -f 1 -d '_')

## merge into counts table
module load R/4.4.0
# base command
cmd="abundance_estimates_to_matrix.pl \
      --est_method RSEM \
      --cross_sample_norm TMM \
      --out_prefix ${alignments}/iso \
      --name_sample_by_basedir"
# find all isoforms.results files in the directory and append them to the command
for file in $(find "${alignments}" -name "*.isoforms.results"); do
  cmd+=" $file"
done
# execute
eval $cmd
```

annotate assembly:
```bash
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

## run simulatanous to transdecoder
blastx -query ${trinity_out}/Trinity.fasta \
  -db ${uniprot}/uniprot_sprot.fasta \
  -out ${annotations}/blastx.outfmt6 \
  -evalue 1e-3 \
  -num_threads 25 \
  -outfmt 6

### filter for top hit
mkdir tmp
sort -T tmp -g -k11,11 ${annotations}/blastx.outfmt6 | sort -T tmp -u -k1,1 > ${annotations}/filter_blastx.fmt6
rm -r tmp/

## signalp on bacteria and eukaryotes
## NOT RUNNING - LIKELY RUNNING OUT OF SPACE
# signalp \
#   -format short \
#   -plot png \
#   -org euk \
#   -tmp tmp \
#   -prefix ${annotations}/signalp_euk.out \
#   -fasta ${transdecoder}/Trinity.fasta.transdecoder.pep
#
# signalp \
#   -format short \
#   -plot png \
#   -org gram- \
#   -tmp tmp \
#   -prefix ${annotations}/signalp_gramn.out \
#   -fasta ${transdecoder}/Trinity.fasta.transdecoder.pep
#
# signalp \
#   -format short \
#   -plot png \
#   -org gram+ \
#   -tmp tmp \
#   -prefix ${annotations}/signalp_gramp.out \
#   -fasta ${transdecoder}/Trinity.fasta.transdecoder.pep


## TO DO - ADD GZIP TO LARGE FILES
```

align reads to transdecoder predicted proteins
```bash
align_and_estimate_abundance.pl \
  --transcripts ${transdecoder}/Trinity.fasta.transdecoder.cds \
  --est_method RSEM \
  --aln_method bowtie2 \
  --trinity_mode \
  --output_dir ${transdecoder} \
  --prep_reference

parallel -j 5 \
  align_and_estimate_abundance.pl \
    --transcripts ${transdecoder}/Trinity.fasta.transdecoder.cds \
    --seqType fq \
    --left ${nohost}/mrna_{1}_rmhost_dog_r1.fq \
    --right ${nohost}/mrna_{1}_rmhost_dog_r2.fq \
    --est_method RSEM \
    --output_dir ${alignments}/mrna2prot_{1} \
    --aln_method bowtie2 \
    --thread_count 10 \
    --gene_trans_map 'none' \
    --output_prefix mrna_{1} ::: $( basename -a ${rawseqs}/*R1*.fastq.gz | cut -f 1 -d '_')

## not working due to gene_trans_map file

# parallel -j 5 \
#   rsem-calculate-expression \
#     --paired-end ${nohost}/mrna_{1}_rmhost_dog_r1.fq ${nohost}/mrna_{1}_rmhost_dog_r2.fq \
#     ${transdecoder}/Trinity.fasta.transdecoder.cds {1} ::: $( basename -a ${rawseqs}/*R1*.fastq.gz | cut -f 1 -d '_')

```


TODO: add Zhang et al 2025 co-expression analysis
TODO: add maaslin3 https://huttenhower.sph.harvard.edu/maaslin3/
TODO: add multimedia analysis using mg/mt/mb https://github.com/krisrs1128/multimedia
TODO: add map coloring for kegg https://www.kegg.jp/kegg/webapp/color_url.html







slurm script template:
```bash
#!/bin/sh
#SBATCH --job-name=align
#SBATCH --qos=defq
#SBATCH --cpus-per-task=10
#SBATCH --mem=50gb
#SBATCH -t 21-00:00:00
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



housekeeping
```bash
for i in *; do du -sh ${i}; done
# 325G	alignments
rm alignments/*.bam.gz
# 16K	analysis.md
# 31G	annotations
gzip annotations/blast*
# 1.1T	bb4_wmgx_wmtx_out
# 15M	data
# 216G	filtered
# 36K	LICENSE
# 692K	locho_study.nb.html
# 8.0K	locho_study.Rmd
# 4.0K	locho_study.Rproj
# 12K	mgmt_analysis.md
# 1.3T	nohost
# 4.0K	plots
# 131M	quality
# 4.0K	README.md
# 24M	scripts
# 411G	sortmerna
# 323G	spades
# 660K	tables
# 0	TMHMM_194952
# 0	tmp
# 126G	trinity_out
```


END

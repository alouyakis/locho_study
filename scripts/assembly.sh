#!/bin/sh
#SBATCH --job-name=assembly
#SBATCH --qos=defq
#SBATCH --cpus-per-task=75
#SBATCH --mem=325gb
#SBATCH -t 21-00:00:00
#SBATCH -o assembly_%j.out
pwd; hostname; date

conda activate mtenv

export trinity_out="/home/artemisl/locho_study/trinity_out"
export nohost="/home/artemisl//locho_study/nohost"

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

date

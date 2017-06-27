#!/bin/bash
# request Bourne shell as shell for job
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -m e
#$ -pe whole_nodes 8
#$ -M klkim@mit.edu
###############################

###----- PARAMETERS -----
ID=160701
sam=6211
#sam1=160701Ess\_D16-6211\_1
#sam2=160701Ess\_D16-6211\_2
###----- END -----

###----- RECOMMENDED DEFAULTS -----
linesPerFile=5000000

###----- END -----


#len1=$(wc -l $sam1\_sequence.fastq) #237821488 for 11
#len2=$(wc -l $sam2\_sequence.fastq) #237821488 for 11
len=$(wc -l $sam1\_sequence.fastq) # since should be equal
numFiles=$(($len / $linesPerFile))
numFiles=${numFiles%.*}
numFiles=$(($numFiles + 1))

split -l $linesPerFile data/$sam1\_sequence.fastq data/$sam1\_ #48 files total
split -l $linesPerFile data/$sam2\_sequence.fastq data/$sam2\_

splits=()
count=1
for l1 in {a..z}; do
  for l2 in {a..z}; do
    if [ $count -gt $numFiles ]; then
      break
    fi
    splits+=($l1$l2)
    count=$(($count + 1))
  done
done

for sam in $sam1 $sam2; do
  for suf in "${splits[@]}"; do
    mv data/$sam\_$suf data/$sam\_$suf.sequence.fastq
  done
done

mkdir pipeline
mkdir pipeline/01

seq1=''
seq2=''
for l in "${splits[@]}"; do
  python DCS-2.00/tag\_to\_header.py --infile1 data/$sam1\_$l.sequence.fastq --infile2 data/$sam2\_$l.sequence.fastq --outfile1 pipeline/01/$sam1\_$l.sequence.fq.smi --outfile2 pipeline/01/$sam2\_$l.sequence.fq.smi
  seq1+=pipeline/01/$sam1\_$l.sequence.fq.smi' '
  seq2+=pipeline/01/$sam2\_$l.sequence.fq.smi' '
done

cat $seq1 > pipeline/01/$sam1\_sequence.fq.smi
cat $seq2 > pipeline/01/$sam2\_sequence.fq.smi

activate=/net/bmc-pub7/data2/essigmannlab/users/klkim/dependencies/anaconda3/bin/activate;
dcs2=/net/bmc-pub7/data2/essigmannlab/users/klkim/dependencies/anaconda3/envs/dcs2/;

export PATH=/home/klkim/dependencies/anaconda3/bin:$PATH

source $activate $dcs2;

snakemake --latency-wait 90;

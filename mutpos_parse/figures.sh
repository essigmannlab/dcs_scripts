#!/bin/bash
# Counter for trinucleotide contexts using Jellyfish
# Input:  fasta file
# Outputs:  text file with context counts
# Dependencies: Python, Jellyfish, dependencies of mutpos_parse_normalized.py

# Default parameters
min_clon=0
max_clon=0.2
min_depth=100
notation='pyrimidine'

# Help message
help="usage: figures.sh -r REFERENCE_FASTA -m MUTPOS_FILE -o OUT_FILE [-c MIN_CLONALITY] [-l MAX_CLONALITY] [-d MIN_DEPTH] [-n NOTATION] [-h HELP]"

while getopts r:m:o:c::l::d::n::h option
do
  case "${option}"
  in
  r) REF_FILE=${OPTARG};;
  m) MUTPOS_FILE=${OPTARG};;
  o) OUT_FILE=${OPTARG};;
  c) min_clon=${OPTARG};;
  l) max_clon=${OPTARG};;
  d) min_depth=${OPTARG};;
  n) notation=${OPTARG};;
  h) echo "$help"
     exit;;
  esac
done

###----- PARAMETERS -----
#REF_FILE='EG10_custom.fasta' # Assumes it has a header line
#MUTPOS_FILE='/media/sf_share/triplicate_analysis/x6614-16.mutpos' # Assumes no header?
###----- END -----

declare -a CONTEXT=('ACA' 'ACC' 'ACG' 'ACT'
                    'ATA' 'ATC' 'ATG' 'ATT'
                    'CCA' 'CCC' 'CCG' 'CCT'
                    'CTA' 'CTC' 'CTG' 'CTT'
                    'GCA' 'GCC' 'GCG' 'GCT'
                    'GTA' 'GTC' 'GTG' 'GTT'
                    'TCA' 'TCC' 'TCG' 'TCT' 
					'TTA' 'TTC' 'TTG' 'TTT');

JF_FILE='jf_counter.jf'
JF_OUT='jf_counted.txt'
rm $JF_FILE
rm $JF_OUT
touch $JF_OUT

jellyfish count -m 3 -s 10M -C $REF_FILE -o $JF_FILE

for i in "${CONTEXT[@]}"
do
  COUNT=$(jellyfish query $JF_FILE $i | cut -d ' ' -f2)
  echo $i $COUNT >> $JF_OUT
  echo $i 'has been counted'
  echo $COUNT
done

python mutpos_parse_normalized.py -r $JF_OUT -a $REF_FILE -m $MUTPOS_FILE -c $min_clon -C $max_clon -d $min_depth -n $notation -o $OUT_FILE

#Rscript trinucleotide-graphs.r $JF_OUT --save

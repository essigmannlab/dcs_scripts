#!/bin/bash
# Counter for trinucleotide contexts using Jellyfish
# Input:  FASTA reference, mutpos file(s) to plot
# Outputs:  text file with context counts; four mutational spectra (combination of total/unique muts, frequencies/proportions); Excel file with substitution counts
# Dependencies: Python, Jellyfish, dependencies of mutpos_update.py

###----- PARAMETERS -----
REF_FILE='data/ref/EG10_custom.fasta' # Assumes it has a header line
declare -a files=('data/a10_ucm.mutpos') # Allows for multiple mutpos files at once
###----- END -----

declare -a CONTEXT=('ACA' 'ACC' 'ACG' 'ACT'
                    'ATA' 'ATC' 'ATG' 'ATT'
                    'CCA' 'CCC' 'CCG' 'CCT'
                    'CTA' 'CTC' 'CTG' 'CTT'
                    'GCA' 'GCC' 'GCG' 'GCT'
                    'GTA' 'GTC' 'GTG' 'GTT'
                    'TCA' 'TCC' 'TCG' 'TCT'
                    'TTA' 'TTC' 'TTG' 'TTT'
		    'AAA' 'AAC' 'AAG' 'AAT'
                    'AGA' 'AGC' 'AGG' 'AGT'
                    'CGA' 'CGC' 'CGG' 'CGT'
                    'CAA' 'CAC' 'CAG' 'CAT'
                    'GGA' 'GGC' 'GGG' 'GGT'
                    'GAA' 'GAC' 'GAG' 'GAT'
                    'TGA' 'TGC' 'TGG' 'TGT'
                    'TAA' 'TAC' 'TAG' 'TAT');

JF_FILE='jf_counter.jf'
OUT_FILE='ref_counts.txt'

touch $OUT_FILE

jellyfish count -m 3 -s 10M -C $REF_FILE -o $JF_FILE

for i in "${CONTEXT[@]}"
do
  COUNT=$(jellyfish query $JF_FILE $i | cut -d ' ' -f2)
  echo $i $COUNT >> $OUT_FILE
  echo $i 'has been counted'
  echo $COUNT
done

for fi in ${files[@]}; do
  python mutpos_update.py -k $OUT_FILE -a $REF_FILE -m $fi -u "unique" -p "frequencies"
  python mutpos_update.py -k $OUT_FILE -a $REF_FILE -m $fi -u "total" -p "frequencies"
  python mutpos_update.py -k $OUT_FILE -a $REF_FILE -m $fi -u "unique" -p "proportions"
  python mutpos_update.py -k $OUT_FILE -a $REF_FILE -m $fi -u "total" -p "proportions"
done

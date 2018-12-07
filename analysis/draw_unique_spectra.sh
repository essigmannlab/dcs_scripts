#!/bin/bash
# Converts mutpos files to extract only unique mutations
# Generate mutational spectra for total and unique mutations

###----- Parameters to change -----

declare -a files=('data/x8896.mutpos'
		  'data/x8897.mutpos'
		  'data/x8898.mutpos');

ref_file='data/ref/EG10_custom.fasta'

###------- End parameters ---------

# Count trinucleotide contexts in reference with Jellyfish
declare -a context=('ACA' 'ACC' 'ACG' 'ACT'
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

jf_file='jf_counter.jf'
count_file='ref_counts.txt'

rm $count_file
touch $count_file

jellyfish count -m 3 -s 10M -C $ref_file -o $jf_file

echo "Counting trinucleotide contexts in the reference sequence..."

for i in "${context[@]}"; do
  count=$(jellyfish query $jf_file $i | cut -d ' ' -f2)
  echo $i $count >> $count_file
done

# Draw mutational spectra for total mutations
echo "Drawing mutational spectra for total mutations..."

for fi in ${files[@]}; do
  python draw_spectra.py -k $count_file -a $ref_file -m $fi -p "frequencies" -n "purine"
  python draw_spectra.py -k $count_file -a $ref_file -m $fi -p "proportions" -n "purine"
done

# Create mutpos files for unique mutations
echo "Creating mutpos files for unique mutations..."

declare -a unique_files=();
for fi in ${files[@]}; do
  Rscript unique-mutpos.R $fi
  char_len=${#fi}
  chars=$((char_len-7))
  unique_files+=(${fi:0:$chars}'-unique.mutpos')
done

# Draw mutational spectra for unique mutations
echo "Drawing mutational spectra for unique mutations..."

for ufi in ${unique_files[@]}; do
  python draw_spectra.py -k $count_file -a $ref_file -m $ufi -p "frequencies" -n "purine"
  python draw_spectra.py -k $count_file -a $ref_file -m $ufi -p "proportions" -n "purine"
done


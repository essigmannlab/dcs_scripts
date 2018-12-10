#!/bin/bash
# Converts mutpos files to extract only unique mutations
# Generate mutational spectra for total and unique mutations

###----- Parameters to change -----

declare -a files=('data/x8896.mutpos'
		  'data/x8897.mutpos'
		  'data/x8898.mutpos');

ref_file='data/ref/EG10_custom.fasta'
group_name='afb1m_25wk'

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

jf_file=$group_name/'jf_counter.jf'
count_file=$group_name/'ref_counts.txt'
total_list=''
unique_list=''

# Set up Conda environment
activate=/mnt/d/software/anaconda3/bin/activate;
spectra=/mnt/d/software/anaconda3/envs/spectra/;
export PATH=/mnt/d/software/anaconda3/bin:$PATH
export PATH=/bin:$PATH
source $activate $spectra;
echo "Virtual environment activated!"

rm -rf $group_name
mkdir $group_name

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
  python draw_spectra.py -k $count_file -a $ref_file -m $fi -o $group_name -p "frequencies" -n "purine"
  python draw_spectra.py -k $count_file -a $ref_file -m $fi -o $group_name -p "proportions" -n "purine"
done

# Create mutpos files for unique mutations
echo "Creating mutpos files for unique mutations..."

declare -a unique_files=();
for fi in ${files[@]}; do
  Rscript unique-mutpos.R $fi
  char_len=${#fi}
  chars=$((char_len-7))
  unique_files+=(${fi:0:$chars}'-unique.mutpos')
  total_list+=$fi','
done

# Draw mutational spectra for unique mutations
echo "Drawing mutational spectra for unique mutations..."

for ufi in ${unique_files[@]}; do
  python draw_spectra.py -k $count_file -a $ref_file -m $ufi -o $group_name -p "frequencies" -n "purine"
  python draw_spectra.py -k $count_file -a $ref_file -m $ufi -o $group_name -p "proportions" -n "purine"
  unique_list+=$ufi','
done

# Combine mutpos files to create combined mutational spectra
# for both total and unique mutations
echo "Combining mutpos files in the sample group..."

total_list=${total_list:0:${#total_list}-1}
unique_list=${unique_list:0:${#unique_list}-1}

python combine_mutpos.py -i $total_list -o $group_name/$group_name'-total'
python combine_mutpos.py -i $unique_list -o $group_name/$group_name'-unique'

# Draw mutational spectra for combined samples
echo "Drawing mutational spectra for combined samples..."

python draw_spectra.py -k $count_file -a $ref_file -m $group_name/$group_name'-total'.mutpos -o $group_name -p "frequencies" -n "purine"
python draw_spectra.py -k $count_file -a $ref_file -m $group_name/$group_name'-total'.mutpos -o $group_name -p "proportions" -n "purine"

python draw_spectra.py -k $count_file -a $ref_file -m $group_name/$group_name'-unique'.mutpos -o $group_name -p "frequencies" -n "purine"
python draw_spectra.py -k $count_file -a $ref_file -m $group_name/$group_name'-unique'.mutpos -o $group_name -p "proportions" -n "purine"

 

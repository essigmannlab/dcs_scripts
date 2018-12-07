#!/bin/bash
# Counter for trinucleotide contexts using Jellyfish
# Input: FASTA reference, mutpos file(s) to plot
# Outputs: text file with context counts; four mutation spectra (combination of total/unique muts, frequencies/proportions); CSV file with substitution counts
# Dependencies: Python, Jellyfish, dependencies of mutpos_update

###----- PARAMETERS -----
ref_file='data/ref/EG10_custom.fasta' # Assumes it has a header line
#declare -a files=('data/a10_ucm.mutpos')
declare -a files=('data/pre-process/x2430.mutpos'
		  'data/pre-process/x2431.mutpos'
		  'data/pre-process/x2432.mutpos'
		  'data/pre-process/x2433.mutpos'
		  'data/pre-process/x2434.mutpos'
		  'data/pre-process/x2435.mutpos'
		  'data/pre-process/x8896.mutpos'
		  'data/pre-process/x8897.mutpos'
		  'data/pre-process/x8898.mutpos'
		  'data/pre-process/x5148.mutpos'
		  'data/pre-process/x5149.mutpos'
		  'data/pre-process/x5150.mutpos');
declare -a combos=('data/pre-process/m_hep_tcpo_25wk.mutpos'
		   'data/pre-process/m_hep_tcpodmso_25wk.mutpos'
		   'data/pre-process/m_hep_afb1_25wk.mutpos'
		   'data/pre-process/m_hep_afb1tcpo_25wk.mutpos');
###----- END -----

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

touch $count_file

jellyfish count -m 3 -s 10M -C $ref_file -o $jf_file

for i in "${context[@]}"
do
  count=$(jellyfish query $jf_file $i | cut -d ' ' -f2)
#  echo $i $count >> $count_file
#  echo $i 'has been counted'
#  echo $count
done

declare -a u_files=(); 
for fi in ${files[@]}; do
  Rscript unique-mutpos.R $fi
  var_len=${#fi}
  chars=$((var_len-7))
  u_files+=(${fi:0:$chars}'-unique.mutpos')
done

for fi in ${files[@]}; do
  python mutpos_update.py -k $count_file -a $ref_file -m $fi -u "total" -p "frequencies" -n "purine"
  python mutpos_update.py -k $count_file -a $ref_file -m $fi -u "total" -p "proportions" -n "purine"
done

# Mutpos files are already unique anyway, go for total
for ufi in ${u_files[@]}; do
  python mutpos_update.py -k $count_file -a $ref_file -m $ufi -u "total" -p "frequencies" -n "purine"
  python mutpos_update.py -k $count_file -a $ref_file -m $ufi -u "total" -p "proportions" -n "purine"
done

fpath=data/pre-process
python data/pre-process/mutpos_combine.py -i $fpath/x2430-unique.mutpos,$fpath/x2431-unique.mutpos,$fpath/x2432-unique.mutpos -o $fpath/m_hep_tcpo_25wk.mutpos
python data/pre-process/mutpos_combine.py -i $fpath/x2433-unique.mutpos,$fpath/x2434-unique.mutpos,$fpath/x2435-unique.mutpos -o $fpath/m_hep_tcpodmso_25wk.mutpos
python data/pre-process/mutpos_combine.py -i $fpath/x8896-unique.mutpos,$fpath/x8897-unique.mutpos,$fpath/x8898-unique.mutpos -o $fpath/m_hep_afb1_25wk.mutpos
python data/pre-process/mutpos_combine.py -i $fpath/x5148-unique.mutpos,$fpath/x5149-unique.mutpos,$fpath/x5150-unique.mutpos -o $fpath/m_hep_afb1tcpo_25wk.mutpos

# Total mutations from unique mutpos files, to not get rid of unique mutations occurring at same pos, but different mice
for cfi in ${combos[@]}; do
  python mutpos_update.py -k $count_file -a $ref_file -m $cfi -u "total" -p "frequencies" -n "purine"
  python mutpos_update.py -k $count_file -a $ref_file -m $cfi -u "total" -p "proportions" -n "purine"
done 

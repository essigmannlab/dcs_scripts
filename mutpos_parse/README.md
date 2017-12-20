# DCS Post-Processing
### _Author_:
Lina Kim, klkim [at] mit dot edu; Essigmann Lab, MIT
### _Version_: 
2.5
### _Date_:
2.5: 12/20/17
2.0: 08/29/17

```
usage: bash figures.sh [-h] -r REFERENCE_FASTA -m MUTPOS_FILE -o OUT_FILE [-c MIN_CLONALITY] [-l MAX_CLONALITY] [-d MIN_DEPTH] [-n NOTATION]

Find the normalized trinucleotide mutation frequencies from parsing a FASTA reference and .mutpos file as outputted from the duplex sequencing pipeline.

optional arguments:
  -h, --help            Show this help message and exit.
  -r REF_FASTA, REF_FILE
                        The reference genome in FASTA format.
  -m MUTPOS_FILE, MUTPOS_FILE
                        The mutpos file.
  -o OUT_FILE, OUT_FILE
                        The output plot name.
  -c MIN_CLONALITY, --minClonality MIN_CLONALITY
                        The minimum clonality for counting a mutation [0].
  -l MAX_CLONALITY, --maxClonality MAX_CLONALITY
                        The maximum clonality for counting a mutation
                        (inclusive) [100].
  -d MIN_DEPTH, --minDepth MIN_DEPTH
                        The minimum depth of reads for each location
                        (inclusive) [100].
  -n NOTATION, --notation NOTATION
                        This is useful for labelling figures as purines or
                        pyrimidines. [pyrimidine]
```
### _Changes from 2.0_:
* Removed 'lab' option, making Essigmann Lab the default
* Changed the default max_clonality value from 0.2 to 100, to account for total mutations
* Calculations for unique mutation placement are made later in the script, by position instead of clonality threshold
* Normalization for context frequencies changed to be made in the spectrum dictionary
* Added functionality: plotting with purine labels in the same order as pyrimidine labels (a la Stratton signatures)
* Context labels moved 0.25 units to the left, for better alignment with bars

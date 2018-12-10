# Mutational Analysis

Draws mutational spectra for given .mutpos files (outputs of DCS pipeline); generates the following files:
* 2 spectra for each sample's total mutations: normalized ("-prop") and unnormalized ("-freq")
* 2 spectra for each sample's unique mutations: normalized ("-prop") and unnormalized ("-freq")

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
                        (inclusive) [1].
  -d MIN_DEPTH, --minDepth MIN_DEPTH
                        The minimum depth of reads for each location
                        (inclusive) [100].
  -n NOTATION, --notation NOTATION
                        This is useful for labelling figures as purines or
                        pyrimidines. [pyrimidine]
```

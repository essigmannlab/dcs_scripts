# DCS Post-Processing
_Author_: Lina Kim, klkim [at] mit dot edu; Essigmann Lab, MIT

_Version_: 2.0

_Date_: 08/29/17

NOTE: first step is **source activate mutpos_parse**

```
usage: bash figures.sh [-h] -r REFERENCE_FASTA -m MUTPOS_FILE -o OUT_FILE [-c MIN_CLONALITY] [-l MAX_CLONALITY] [-d MIN_DEPTH] [-n NOTATION]

Find the normalized trinucleotide mutation frequencies from parsing a FASTA reference and .mutpos file as outputted from the duplex sequencing pipeline.

optional arguments:
  -h, --help            show this help message and exit
  -r REF_FASTA, REF_FILE
                        The reference genome in FASTA format.
  -m MUTPOS_FILE, MUTPOS_FILE
                        The mutpos file.
  -o OUT_FILE, OUT_FILE
                        The output plot.
  -c MIN_CLONALITY, --minClonality MIN_CLONALITY
                        The minimum clonality for counting a mutation [0]
  -l MAX_CLONALITY, --maxClonality MAX_CLONALITY
                        The maximum clonality for counting a mutation
                        (inclusive) [0.2]
  -d MIN_DEPTH, --minDepth MIN_DEPTH
                        The minimum depth of reads for each location
                        (inclusive) [100].
  -n NOTATION, --notation NOTATION
                        This is useful for labelling figures aspurines or
                        pyrimidines. [pyrimidine]
```

```bash
usage: mcs.py [-h] -r REF_FILE -m MUTPOS_FILE [-c MIN_CLONALITY]
              [-C MAX_CLONALITY] [-d MIN_DEPTH] [-n NOTATION] [-i DPI]
              [-f FORMAT] [-y YMAX] [-s SAVE] [-t TITLE] [-l LAB]

Find the non-normalized trinucleotide mutation frequencies from parsing a
FASTA reference and .mutpos file as outputted from the duplex sequencing
pipeline.

optional arguments:
  -h, --help            show this help message and exit
  -r REF_FILE, --ref_file REF_FILE
                        The reference genome in FASTA format.
  -m MUTPOS_FILE, --mutpos_file MUTPOS_FILE
                        The mutpos file.
  -c MIN_CLONALITY, --minClonality MIN_CLONALITY
                        The minimum clonality for counting a mutation [0]
  -C MAX_CLONALITY, --maxClonality MAX_CLONALITY
                        The maximum clonality for counting a mutation
                        (inclusive) [0.2]
  -d MIN_DEPTH, --minDepth MIN_DEPTH
                        The minimum depth of reads for each location
                        (inclusive) [100].
  -n NOTATION, --notation NOTATION
                        This is useful for labelling figures aspurines or
                        pyrimidines. [pyrimidine]
  -i DPI, --dpi DPI     The dots per inch of the final saved figure [320].
  -f FORMAT, --format FORMAT
                        The format of the output image file. Popular options
                        are png or svg [png].
  -y YMAX, --ymax YMAX  Sets the Y_MAX value in figure. Must be 0-100.
  -s SAVE, --save SAVE  Whether to save both, ratio, or total [both].
  -t TITLE, --title TITLE
                        Plot title to be overlayed if supplied. The default
                        option creates an informative title. Type None for no
                        title.
  -l LAB, --lab LAB     Laboratory running program in [loeb].
```

blastfind
==============

Copyright 2014 Clint Valentine

GitHub, Inc. Repository https://github.com/clintval/blastfind

Download: BioPython, SciPy, & NCBI BLAST+ Executables

Open this readme *RAW* or through Notepad

This program is designed to compare a query sequence with a subject sequence
using BLASTn algorithms (somewhat similar sequences) to identify all types
of mutations in the subject sequence. Two files are generated from this code:

1. Comma Separated Value Raw data BLASTn output of all mutations
2. Comma Separated Value Transition & Transversion Analysis output

This program is specific to the Massachusetts Institute of Technology Essigmann
Lab GPT Protocol where sequences are analyzed considering the categories gender,
sample cell type, and treatment.

All input .FASTA or .seq files should be formatted in the following way:

1st number group indicates mouse number
        |
        10H-1.seq
           ||
           | 2nd number group indicates mutant sequenced
          (-) indicates break between mouse number and mutant number

Cell type letters can be found grouped anywhere in the filename. Cell types that are
currently supported are H (hepatocyte), NH (non-hepatocyte), & (T) tumor.

1. Install Python 2.7
2. Install NumPy
3. Install SciPy
4. Install BioPython
5. Install Blast+ Executable
6. Make sure Python and Blast+ are appended to PATH
7. Download repository (all files) into 'working folder'. Test with sample file TF1-1_Plasmid.seq.
8. Open terminal run->cmd
9. Edit queryseq file to change subject sequence for alignment. Currently GPT transgene fragment.
10. Change directory to folder with blastfind scripts and supporting files (ex: cd C:/scripts)
11. Type: python blastfind.py and hit _Enter_
12. Mutation_Results.csv will be created in directory
13. Mutations_Analysis.csv will be created in directory
14. Save all .csv to Excel files before manipulating data.
15. You can run multiple sets of .seq files into the same Mutation_Results.csv. They will merely append. Be sure to delete already ran .seq files from the directory so they are not run twice.
16. Mutation_Analysis.csv will always be overwritten.

Sample Sequence Results (TF1-1_Plasmid.seq)
A deletion at position 12
G to A transition at position 110

###License

See LICENSE.md

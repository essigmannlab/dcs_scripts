#!/bin/env python3.4
# Copyright 2014 Clint Valentine
# GitHub, Inc. Repository https://github.com/clintval/blastfind
# Download: Python 2.7, BioPython, & NCBI BLAST+ Executables

"""
This script is designed to compare a query sequence with a subject
sequence using BLASTn algorithms (somewhat similar sequences) to
identify all types of mutations in the query sequence. Two files are
generated from this code:

1. Comma Separated Value Raw data BLASTn output of all mutations
2. Comma Separated Value Transition & Transversion Analysis output

This program is specific to the Massachusetts Institute of Technology
Essigmann Lab GPT Protocol where sequences are analyzed considering the
categories gender, sample cell type, and toxins. No other categories are
supported at this time. To accomodate backwards compatibilty with
previously saved sequence data, only mouse number, sequence number, and
celltype are parsed from the filename.

All input .FASTA or .seq files should be formatted in the following way:

1st number group indicates mouse number
    |
    10H-1.seq/fasta/fsa/fa
       ||
       | 2nd number group indicates mutant sequenced
      (-) indicates break between mouse number and mutant number

XML Blast output files are also valid data inputs and are recommended
due to the time processing of the BLASTn algorithms to create them.

Cell type letters can be found grouped anywhere in the filename.
Cell types that are currently supported are (H) hepatocyte, (NH) non-
hepatocyte, & (T) tumor.
"""

import sys
import os
import re
import glob
import csv

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Blast import NCBIXML as xmlread
    from Bio.Blast.Applications import NcbiblastnCommandline as blastn
except:
    print ("You must download Biopython & NCBI Blast Executables",
           "Append both to PATH before running blastfind.py script")
    raise SystemExit

__version__ = '08.09.2014'

QNAME = "GPT_TRANSGENE_FRAGMENT"
QDESCRIPTION = "Essigmann Lab MIT GPT Protocol Query"
QSEQ = ('atgagcgaaaaatacatcgtcacctgggacatgttgcagatccatgcacgtaaactcgca'
        'agccgactgatgccttctgaacaatggaaaggcattattgccgtaagccgtggcggtctg'
        'gtaccgggtgcgttactggcgcgtgaactgggtattcgtcatgtcgataccgtttgtatt'
        'tccagctacgatcacgacaaccagcgcgagcttaaagtgctgaaacgcgcagaaggcgat'
        'ggcgaaggcttcatcgttattgatgacctggtggataccggtggtactgcggttgcgatt'
        'cgtgaaatgtatccaaaagcgcactttgtcaccatcttcgcaaaaccggctggtcgtccg'
        'ctggttgatgactatgttgttgatatcccgcaagatacctggattgaacagccgtgggat'
        'atgggcgtcgtattcgtcccgccaatctccggtcgctaa')

def cls():
    """Clears terminal"""
    os.system('cls' if os.name == 'nt' else 'clear')

def meter(process, progress, length=30):
    """Provides visual indicator of progress"""
    block = int(round(length * progress))
    disp = "\r{} [{}] {}%".format(process,
                                  '=' * block + '>' + ' ' * (length - block),
                                  int(progress * 100))
    sys.stdout.write(disp)
    sys.stdout.flush()

def createquery(query, qid, description):
    """Creates .fasta file in input directory for BLAST executables"""
    queryseq = SeqRecord(Seq(query), id=qid, description=description)
    SeqIO.write(queryseq, "queryseq", "fasta")

def xmlspitter():
    """Creates XML blast output files from sequence data in cwd"""
    for seq in os.listdir('.'):
        if re.match(r'.*\.(seq|fasta|fsa|fa)$', seq, re.IGNORECASE):
            name, _ = seq.split('.')
            blastn(query="queryseq", subject=seq, outfmt=5,
                   out=str(name.upper()) + ".xml")()

def naming(reg={}):
    """
    Identifies categorical information for each sample. Directory example:
    {filename: [idmouse, idsample, gender, celltype, toxins]}
    The reg dict keeps track of samples to eliminate redundant queries
    """
    cls()
    directory = {key: ['' for _ in xrange(5)] for key in glob.iglob('*.xml')}

    def sortby(tuplesort):
        """Sort directory dictionary by idmouse then by idsample"""
        return (tuplesort[1][0], tuplesort[1][1])

    for key, val in directory.items():
        idmouse, idsample = autoname(key)
        val[0] = idmouse
        val[1] = idsample
        val[3] = autocelltype(key)
    for i, (key, val) in enumerate(sorted(directory.items(), key=sortby)):
        meter("Input Manual Data: ", (i + 1) / float(len(directory)))
        lookup = reg[val[0]] = reg[val[0]] if val[0] in reg else ask(val[0])
        val[2], val[4] = lookup
    return directory

def autoname(name):
    """Parses sample idmouse + idsample from filename"""
    first, second = name.split('-')
    return (int(re.search(r'\d+', first).group()),
            int(re.search(r'\d+', second).group()))

def autocelltype(name):
    """Parses celltype from filename"""
    if 'NH' in name.upper():
        return 'NH'
    elif 'H' in name.upper():
        return 'H'
    elif 'T' in name.upper():
        return 'T'

def ask(name):
    """Queries user for manual assignment of categories"""
    gender, toxins = ('.', '.')
    while gender not in 'FfMm':
        gender = raw_input('\nGender F/M of: {}\n\n--> '.format(name))
    while toxins not in 'YyNn':
        toxins = raw_input('\nToxin Y/N of: {}\n\n--> '.format(name))
    toxins = 'Dosed' if toxins in 'Yy' else 'Control'
    cls()
    return (gender.upper(), toxins)

def mutfinder(directory, mutlist=[]):
    """Uses BLAST commandline to find point mutations in xml files"""
    delta = {key: 0 for key in directory}
    for i, blast in enumerate(directory):
        meter("Finding Mutations: ", (i + 1) / float(len(delta)))
        with open(blast) as data:
            parsed = xmlread.parse(data)
            for record in parsed:
                for alig in record.alignments:
                    for hsp in alig.hsps:
                        for pos, mask in enumerate(hsp.match):
                            if mask == ' ':
                                mutline = list(directory[blast])
                                mutline.extend([str(hsp.query[pos]),
                                                str(hsp.sbjct[pos]),
                                                str(pos + delta[blast] + 1)])
                                mutwriter(mutline)
                                mutlist.append(mutline)
                                delta[blast] += inadjuster(mutline[5])
    return mutlist

def inadjuster(query):
    """Returns -1 integer for scanning shift due to insertion only"""
    return 0 if query == '-' else -1

def mutationsummer(mutlist):
    """Bins instances of point mutation type per category of samples"""
    cls()
    categories = ('FHDosed FHControl FNHDosed FNHControl FTDosed FTControl '
                  'MHDosed MHControl MNHDosed MNHControl MTDosed MTControl'
                  .split())
    muttypes = ('TRI_1 TRI_2 TRV_1 TRV_2 TRV_3 TRV_4 Ins_1 Del_1'.split())
    masterdict = {key: {key: 0 for key in categories} for key in muttypes}
    dmap = {
        ('G', 'A'): 'TRI_1', ('C', 'T'): 'TRI_1',
        ('A', 'G'): 'TRI_2', ('T', 'C'): 'TRI_2',
        ('G', 'T'): 'TRV_1', ('C', 'A'): 'TRV_1',
        ('G', 'C'): 'TRV_2', ('C', 'G'): 'TRV_2',
        ('A', 'T'): 'TRV_3', ('T', 'A'): 'TRV_3',
        ('A', 'C'): 'TRV_4', ('T', 'G'): 'TRV_4',
        ('-', 'A'): 'Ins_1', ('-', 'C'): 'Ins_1',
        ('-', 'G'): 'Ins_1', ('-', 'T'): 'Ins_1',
        ('A', '-'): 'Del_1', ('C', '-'): 'Del_1',
        ('G', '-'): 'Del_1', ('T', '-'): 'Del_1'}
    for i, mutline in enumerate(mutlist):
        key = ''.join([_ for _ in mutline[2:5]])
        query = mutline[5]
        subject = mutline[6]
        masterdict[dmap[query, subject]][key] += 1
        meter("Analyzing Frequencies: ", (i + 1) / len(mutlist))
    cls()
    analysiswriter(masterdict, categories)

def mutwriter(mutline):
    """Writes mutation results to a .txt file in directory of script"""
    with open('Mutation_Results.txt', 'ab') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(mutline)

def analysiswriter(masterdict, categories):
    """Writes analysis results to a .txt file in directory of script"""
    with open('Mutation_Analysis.txt', 'wb') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(', G:C to A:T, A:T to G:C, G:C to T:A'
                        ', G:C to C:G, A:T to T:A, A:T to C:G'
                        ', Ins_1, Del_1'.split(', '))
        for key in categories:
            writer.writerow([key, masterdict['TRI_1'][key],
                             masterdict['TRI_2'][key],
                             masterdict['TRV_1'][key],
                             masterdict['TRV_2'][key],
                             masterdict['TRV_3'][key],
                             masterdict['TRV_4'][key],
                             masterdict['Ins_1'][key],
                             masterdict['Del_1'][key]])

if __name__ == '__main__':
    createquery(QSEQ, QNAME, QDESCRIPTION)
    xmlspitter()
    mutationsummer(mutfinder(naming()))
    raise SystemExit

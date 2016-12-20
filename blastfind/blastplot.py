#!/usr/bin/python
# -*- coding: utf-8 -*-
# Copyright 2014 Clint Valentine
# Disable Pylint Numpy Errors:
# pylint: disable=E1101
# pylint: disable=E1103
"""
This script is designed to plot a spectral profile of all the mutations
found from the output source of blastfind.py.

Development in progress...

To do
1a. Consider color stacked histogram relative to mutation type
1b. Show legend for stacked mutation types
2. Incorporate responsive toggles for expand/minimzing selection on gene
3. Find clever solution to illustrate indel events
4a. Output a Monte Carlo Bayesian report on control vs. dosed
4b. Hypergeometric test?
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

__version__ = '2014.08.25'

def framefromcsv():
    """Reads csv file"""
    frame = pd.read_csv(r'Mutation_Results.csv',
                        header=0, sep=r',', index_col=[0, 1],
                        skipinitialspace=True)
    frame['categories'] = frame.apply(category_return, axis=1)
    return frame

def datasort(frame, category, case):
    """Returns a subset of the main dataframe based on case for category"""
    genelength = frame.position.max()
    data = frame[frame[category] == case]
    data = pd.crosstab(data['position'], data['categories'])
    data = data.reindex([_ + 1 for _ in range(genelength)], fill_value=0)
    return data

def category_return(dtf):
    """Provide query and subject values to dmap for category type return"""
    dmap = {('G', 'A'): 'TRI_1', ('C', 'T'): 'TRI_1',
            ('A', 'G'): 'TRI_2', ('T', 'C'): 'TRI_2',
            ('G', 'T'): 'TRV_1', ('C', 'A'): 'TRV_1',
            ('G', 'C'): 'TRV_2', ('C', 'G'): 'TRV_2',
            ('A', 'T'): 'TRV_3', ('T', 'A'): 'TRV_3',
            ('A', 'C'): 'TRV_4', ('T', 'G'): 'TRV_4',
            ('-', 'A'): 'Ins', ('-', 'C'): 'Ins',
            ('-', 'G'): 'Ins', ('-', 'T'): 'Ins',
            ('A', '-'): 'Del', ('C', '-'): 'Del',
            ('G', '-'): 'Del', ('T', '-'): 'Del'}
    return dmap[dtf['query'], dtf['subject']]

def histoplot(xpos, xneg, case1, case2):
    """Plots histogram for mutation spectra"""
    _, axes = plt.subplots()
    cplot = axes.bar(xneg.index, -(xneg.sum(axis=1)), width=1,
                     color='k', edgecolor="none")
    dplot = axes.bar(xpos.index, xpos.sum(axis=1), width=1,
                     color='k', edgecolor="none")

    def autolabel(rects, delta=1):
        """attach some text labels"""
        for rect in rects:
            position = rect.get_x()
            height = (rect.get_height() * delta)
            if height > 8:
                axes.text(position, height + 1, '%d'%int(position),
                          rotation=90, ha='center', va='bottom')
            elif height < -8:
                axes.text(position, height - len(str(position)) - 1,
                          '%d'%int(position), rotation=90,
                          ha='center', va='bottom')

    plt.title('Mutational Spectra of GPT Gene')
    leg1 = plt.legend([case1], loc=1)
    plt.legend([case2], loc=4)
    plt.gca().add_artist(leg1)
    plt.xlabel('Position')
    plt.ylabel('Frequency')
    plt.xticks(np.arange(0, 500, 50))
    plt.yticks(np.arange(-20, 61, 10),
               ('-20', '-10', '0', '10', '20', '30', '40', '50', '60'),)
    plt.axhline(0, color='r')
    autolabel(dplot)
    autolabel(cplot, -1)
    spines_to_remove = ['top', 'right']
    for spine in spines_to_remove:
        axes.spines[spine].set_visible(False)
    axes.xaxis.set_ticks_position('bottom')
    axes.yaxis.set_ticks_position('left')
    plt.minorticks_on()
    axes.grid(which='both', color='0.08', linestyle=':')
    plt.show()

def hypergprep(dosed, control):
    """Prepares dat file for Hyperg analysis"""
    statframe = pd.DataFrame()
    statframe['dsum'] = dosed.sum(axis=1)
    statframe['csum'] = control.sum(axis=1)
    statframe = statframe[(statframe['dsum'] != 0) | (statframe['csum'] != 0)]
    statframe.to_csv('data.dat', columns=False, index=False, sep='\t')
    return statframe

if __name__ == '__main__':
    CATEGORY = 'treatment'
    CASE1 = 'Dosed'
    CASE2 = 'Control'

    DATA = framefromcsv()
    XPOS = datasort(DATA, CATEGORY, CASE1)
    XNEG = datasort(DATA, CATEGORY, CASE2)
    #FRAME = hypergprep(XPOS, XNEG)
    histoplot(XPOS, XNEG, CASE1, CASE2)

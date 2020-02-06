#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Generate bins from sorted confspace index

Info on module.
"""

##################################################
# Imports
##################################################
import argparse
import os
import sys
import re
import itertools

##################################################
# Info
##################################################
__author__ = "Graham Holt"
__version__ = "0.0.0"
__status__ = "Prototype"

##################################################
# Globals
##################################################

##################################################
# Classes / Functions
##################################################

# Main
def main(args):

    pattern = r"[eE]\+(\d+)"
    def grouper(x):
        n = args.binsize
        match = re.search(pattern, x)
        exp = int(match.group(1))
        if exp == 0:
            return 0
        else:
            return exp//n
    # Set up output dir
    name = './bins_'+str(args.binsize)
    try:
        os.mkdir(name)
    except:
        print('Error! dir already exists')
        raise OSError
    # Setup
    with open(args.index, 'r') as f:
        lines = f.readlines()

    # Split the list into groups based on exponent
    bins = [list(g) for k, g in itertools.groupby(lines, grouper)]
    for e in bins:
        first = re.search(pattern, e[0]).group(1)
        last = re.search(pattern, e[-1]).group(1)
        with open(name+'/bin_'+str(first)+'_'+str(last)+'.csv', 'w') as f:
            for line in e:
                f.write(line)


##################################################
# Conditional main
##################################################

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='This description is shown when -h or --help are passed as arguments.')

    parser.add_argument('index',
            type = str,
            default = None,
            help = 'The file containing sorted cfs index')

    parser.add_argument('binsize',
            type = int,
            default = None,
            help = 'The exponent size to bin by')
    args = parser.parse_args()

    main(args)

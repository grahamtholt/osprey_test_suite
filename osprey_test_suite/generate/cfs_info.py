#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Populate a confspace index with information

Info on module.
"""

##################################################
# Imports
##################################################
import argparse
import os
import sys
import re

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
    # Setup
    pattern = r"_(\d+)res_(.+)\.cfs"
    with open(args.filename, 'r') as f:
        lines = f.readlines()

    info_list = []
    with open(args.filename+".info", 'w') as f:
        for line in lines:
            m = re.search(pattern, line)
            f.write(','.join([line.strip(),m.group(1),m.group(2)])+'\n')

    with open(args.filename+".info.noname", 'w') as f:
        for line in lines:
            m = re.search(pattern, line)
            f.write(','.join([m.group(1),m.group(2)])+'\n')

##################################################
# Conditional main
##################################################

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='This description is shown when -h or --help are passed as arguments.')

    parser.add_argument('filename',
            type = str,
            default = None,
            help = 'The file containing .cfs file names')

    args = parser.parse_args()

    main(args)

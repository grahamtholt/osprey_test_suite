#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module Docstring

Info on module.
"""

##################################################
# Imports
##################################################
import argparse
import os
import sys
import glob
from shutil import move
import pymol

import pymol_tools

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
    orig_dir = os.getcwd()
    dir_name, pdb_name = os.path.split(os.path.abspath(args.pdb))
    print(dir_name, pdb_name)
    out_dir = os.path.join(orig_dir, str(pdb_name)[:4].lower())
    os.mkdir(out_dir)
    os.chdir(dir_name)

    # Launch PyMol
    # Because we might be using an old PyMol version
    pymol.pymol_argv = ['pymol','-qcK']
    pymol.finish_launching()

    # load pdb
    pymol.cmd.load(pdb_name)
    pymol_tools.test_gen()

    # Copy .cfs files from current directory to pdb directory
    for f in glob.glob("*.cfs"):
        move(f,os.path.join(out_dir,f))

    # Return to original directory
    os.chdir(orig_dir)
    # return exit code
    sys.exit()

##################################################
# Conditional main
##################################################

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='This description is shown when -h or --help are passed as arguments.')

    parser.add_argument('pdb',
            type = str,
            default = None,
            help = 'The pdb file from which to generate test designs')

    args = parser.parse_args()

    main(args)

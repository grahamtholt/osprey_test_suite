#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""test_wrapper.py

Script used to submit design tests using the slurm batch system
"""

##################################################
# Version
#   0.0
# Author
#   Graham Holt
##################################################

#TODO: Make directories configurable, not hardcoded

##################################################
# Imports
##################################################
import argparse
from datetime import date
import os
from shutil import copy

##################################################
# Globals
##################################################
FRAM_DIR = "/usr/project/dlab/Users/gth/projects/osprey_test_suite/scripts"

DEBUG_DIR = "/usr/project/dlab/Users/gth/projects/SHARKStar/debug_tests"

PROD_DIR = "/usr/project/dlab/Users/gth/projects/SHARKStar/production_tests"

CFS_DIR = "/usr/project/dlab/Users/gth/projects/SHARKStar/cfs_lists"
##################################################
# Classes / Functions
##################################################

def init(dir_name, cfs_name, debug):
    """Initialize directory structure and read files
    """
    try:
        if debug:
            the_dir = os.path.join(DEBUG_DIR, dir_name)
        else:
            the_dir = os.path.join(PROD_DIR, dir_name)
        os.makedirs(the_dir)
    except OSError as e:
        raise

    copy(os.path.abspath(cfs_name), the_dir)
    os.chdir(the_dir)

# Main
def main(args):
    # Get the number of designs
    num_lines = sum(1 for line in open(os.path.abspath(args.cfsfile)))
    init(args.output, args.cfsfile, args.debug)
    # Call out to a slurm bash script
    os.system("sbatch --array=1-%d%%5 %s/slurm_array.sh %s %d" \
              % ( num_lines,
                  FRAM_DIR,
                  os.path.split(args.cfsfile)[1],
                  int(args.bbk)
                ))

##################################################
# Conditional main
##################################################

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='This description is shown when -h or --help are passed as arguments.')

    parser.add_argument('cfsfile',
                        type = str,
                        default = None,
                        help = "The CFS file from which to run tests")
    parser.add_argument('-d', "--debug",
                        help = "switch output to debug_tests directory",
                        action="store_true"
                       )
    parser.add_argument('-o', "--output",
                        help = "specify output directory name",
                        type = str,
                        default = '_'.join((date.today().strftime("%y%m%d"),
                                            "tests"))
                       )
    parser.add_argument('-b', "--bbk",
                        help = "use the BBK* design algorithm",
                        action="store_true"
                       )


    args = parser.parse_args()

    main(args)

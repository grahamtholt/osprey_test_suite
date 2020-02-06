#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module Docstring

Info on module.
"""

##################################################
# Version
#   0.0
# Author
#   Graham Holt
##################################################

##################################################
# Imports
##################################################
import argparse
from os import path
from os import rename
import re
import pymol_tools
import imp
import sys

##################################################
# Globals
##################################################

##################################################
# Classes / Functions
##################################################

# Main
def main(args):
    # import confspace file
    fn = path.abspath(args.cfsfile)
    if fn[-4:] != ".cfs":
        sys.exit()
    design = imp.load_source('design', fn)

    # calculate confspace size
    design_size = pymol_tools.get_confspace_size(design.strand_flex)

    #format output
    size_string = "_%.3E" % design_size

    # generate new name
    file_dir = path.dirname(fn)
    name = path.basename(fn)
    size_pattern = r"_\d\.\d+E\+\d+"
    new_name = re.sub(size_pattern, size_string, name)

    new_fn = path.join(file_dir, new_name)

    # rename file
    rename(fn, new_fn)
    # return exit code

##################################################
# Conditional main
##################################################

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='This description is shown when -h or --help are passed as arguments.')

    parser.add_argument('cfsfile',
            type = str,
            default = None,
            help = "The confspace file for which to correct the size")

    args = parser.parse_args()

    main(args)

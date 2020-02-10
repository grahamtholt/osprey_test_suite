#!/usr/project/dlab/Users/gth/code/osprey/virtual_envs/sharkstar/bin/python
# -*- coding: utf-8 -*-
"""run_sharkstar.py

Script to run the SHARK* algorithm on a design project.

Takes as input a confspace (.cfs) file.

Outputs a JSON file containing the sequences, K* scores, and upper and lower
bounds on the partition functions for each state.
"""

##################################################
# Version
#   0.0
# Author
#   Graham Holt
##################################################

import argparse
from osprey_test_suite.design import design
import socket
import sys
import osprey
from math import floor
import json
import os


# Main
def main(args):
    # Record relevant parameters
    data = dict.fromkeys(design.HEADERS)
    print args
    data["status"] = design.STATUS.ESTIMATING.value
    data["slurm out"] = args.slurm_out
    data["slurm err"] = args.slurm_err
    data["design name"] = design.get_emat_name(args.cfsfile)
    data["host"] = socket.gethostname()
    data["memory (GB)"] = design.HEAP_GB
    data["epsilon"] = args.epsilon

    data["algorithm"] = design.ALGO_LIST[args.algo]

    # Print out misc information
    print("Output_file: %s" % args.output)
    print("Emat id for this design: %s" % data["design name"])
    print("Running on host: %s" % socket.gethostname())
    sys.stdout.flush()

    osprey.start(heapSizeMiB=floor(design.HEAP_GB * 953.674)) # Convert to MiB
    conf_spaces = design.make_confspace(args.cfsfile, data)
    shark = design.setup_design(conf_spaces, args.epsilon, args.numseqs,
                                args.algo, data)

    # Write out JSON file before beginning computation
    with open(os.path.abspath(args.output), 'w') as f:
        f.write(json.dumps(data, indent = 4, sort_keys=True)+
                '\n')

    design.run(shark, data)

    # Print out JSON output
    with open(os.path.abspath(args.output) , 'w') as f:
        f.write(json.dumps(data ,indent = 4, sort_keys=True)+'\n')

    print("Finished writing")
    # return exit code

##################################################
# Conditional main
##################################################
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='This description is shown when -h or --help are passed as arguments.')

    parser.add_argument('cfsfile',
                        type = str,
                        default = None,
                        help = "The input confspace file")
    parser.add_argument('output',
                        type = str,
                        default=None,
                        help = 'the output json file')
    parser.add_argument('-e', '--epsilon',
                        type = float,
                        default = 0.68,
                        help = "epsilon")
    parser.add_argument('-s', '--numseqs',
                        type = int,
                        default = 5,
                        help = "number of sequences")
    parser.add_argument('--slurm-out',
                        type=str,
                        default=None,
                        help = "slurm standard output")
    parser.add_argument('--slurm-err',
                        type=str,
                        default=None,
                        help = "slurm standard error")
    parser.add_argument('-a', "--algo",
                        help = "choose algorithm: <0> SHARK* <1> MARK* <2> BBK*",
                        type = int,
                        default = 0,
                       )

    args = parser.parse_args()

    main(args)

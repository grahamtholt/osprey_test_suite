#!/usr/project/dlab/Users/gth/code/osprey/virtual_envs/sharkstar_refactor/bin/python
######!/usr/project/dlab/Users/gth/code/osprey/virtual_envs/sharkstar/bin/python
# -*- coding: utf-8 -*-
"""run_design_single_pfunc.py

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
from datetime import datetime
import jpype


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

    osprey.start(heapSizeMiB=floor(design.HEAP_GB * 953.674),
                                  allowRemoteManagement=False) # Convert to MiB
    conf_spaces = design.make_confspace(args.cfsfile, data)

    # Write out JSON file before beginning computation
    with open(os.path.abspath(args.output), 'w') as f:
        f.write(json.dumps(data, indent = 4, sort_keys=True)+
                '\n')

    print("\nRunning Algorithm")
    design_start = datetime.now()
    data['start_time'] = design_start.strftime(design.DATETIME_FORMAT)

    try:
        pfunc = design.make_complex_pfunc(
            args.num_cores,
            conf_spaces,
            args.epsilon,
            args.numseqs,
            args.algo, data)

        if args.algo == 2:
            pfunc.setReportProgress(True)
            pfunc.compute(2147483647)
        elif args.algo==1:
            pfunc.compute()

        data["status"] = design.STATUS.FINISHED.value

    except (jpype.JavaException, AttributeError) as e:
        data["status"] = design.STATUS.ERROR.value
        print str(e)
        print e.stacktrace()
        data['errmsg'] = str(e)

    design_end = datetime.now()
    data['end_time'] = design_end.strftime(design.DATETIME_FORMAT)
    data['runtime (s)'] = (design_end - design_start).total_seconds()
    print("Algorithm Runtime (s): %d" % data['runtime (s)'])

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
    parser.add_argument('-c', "--num-cores",
                        help = "the number of cpus",
                        type=int,
                        default=20,
                       )

    args = parser.parse_args()

    main(args)

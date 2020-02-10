import csv
import json
from os import listdir
from os.path import isfile, join, splitext
import re
import sys
from osprey_test_suite.analyze.result import Result as R
import matplotlib.pyplot as plt
from osprey_test_suite.design.design import STATUS

TEST_FILE = "../tmp/test_results.json"
TEST_DIR = "../191031_tests/outputs"
CFS_SIZE_PATTERN = re.compile(r"_(\d\.\d\d\dE\+\d\d)\.\S*$")

def load_results( dirname='./' ):
    """Loads json test cases from a directory into memory.
    Returns a list of Result objects
    """
    files = [ join(dirname,fp) for fp in listdir(dirname)
                 if isfile(join(dirname, fp)) and splitext(fp)[-1] == '.json']
    data = [ R.from_file(open(f)) for f in files ]
    return data

def report ( results_list ):
    """Print a basic report on results in the input list
    """
    finished = [e for e in results_list if e.status == "FINISHED"]
    estimating = [e for e in results_list if e.status == "ESTIMATING"]
    errors = [e for e in results_list if e.status == "ERROR"]

    print("""Out of %d total tests:
                %d are finished running
                %d are still running (or hit the memory limit)
                %d threw an error
          """\
          % (len(results_list), len(finished), len(estimating), len(errors)))

def print_errors( results_list ):
    """Print any error messages
    """
    errors = [e for e in results_list if e.status == "ERROR"]
    errlist = [ '\t'.join((e.__dict__.get('slurm out'), e.design_name, e.errmsg)) for e in errors]
    print('\n'.join(errlist))

def group( results, f=lambda x: x.design_name ):
    """Groups designs by property

    Takes in a list of results and groups them into tuples by the specified
    accessor function.

    Args:
        results: An iterable containing Results objects
        f: A function returning the property to group by

    Returns:
        A list of tuples. Each tuple contains all results that are equal in the
        specified property.
    """
    keys = set([ f(e) for e in results ])
    grouped = dict.fromkeys(keys)

    for r in results:
        if grouped[f(r)] is None:
            grouped[f(r)] =[r]
        else:
            grouped[f(r)] = grouped[f(r)] + [r]

    return [ tuple(e) for e in grouped.values() ]

def plot_all( results, x_accessor=lambda x: x.get_cfssize(), y_accessor=lambda x: x.runtime, group_func=lambda x: x.algorithm ):
    """Plot data from all algorithms against a design property
    """
    SYMBOL_LIST = ['go', 'bo', 'mo', 'gx', 'bx', 'mx']
    STATUS_SYMBOLS = {'ERROR': 'r^', 'ESTIMATING': 'y^'}
    MARKER_SIZE = 7.5

    grouped = group(results, f=group_func)
    print("there are %d groups" % (len(grouped)))

    counter = 0
    for g in grouped:
        # Plot different statuses differently
        for status in STATUS:
            f = [e for e in g if e.status == status.value]
            if status.value == 'FINISHED':
                plt.loglog([x_accessor(e) for e in f],
                           [y_accessor(e) for e in f],
                           SYMBOL_LIST[counter],
                           label=group_func(f[0]),
                           markersize=MARKER_SIZE
                          )
            else:
                plt.loglog([x_accessor(e) for e in f],
                           [y_accessor(e) for e in f],
                           STATUS_SYMBOLS[status.value],
                           markersize=MARKER_SIZE
                          )
                # plotting again for differentiation
                plt.loglog([x_accessor(e) for e in f],
                           [y_accessor(e) for e in f],
                           SYMBOL_LIST[counter],
                           markersize=0.5*MARKER_SIZE
                          )


        counter = counter + 1

    plt.legend(loc='upper left')
    plt.show()

def plot_vs( results,
            y_accessor=lambda x: x.runtime,
            group_func=lambda x: x.design_name,
            x_selector= lambda x: x.algorithm == 'BBK'):
    """Plot data against reference data
    """
    ## not currently implemented correctly
    raise NotImplementedError()
    ##
    SYMBOL_LIST = ['go', 'bo', 'mo', 'gx', 'bx', 'mx']
    STATUS_SYMBOLS = {'ERROR': 'r^', 'ESTIMATING': 'y^'}
    MARKER_SIZE = 7.5

    def selector(tup, prop_selector):
        x = [ e for e in tup if prop_selector(e)]
        return x

    algorithms = ([e[0] for e in group(results, lambda x: x.algorithm)])
    grouped = group(results, group_func)
    x_results = [ selector(g, x_selector) for g in grouped ]

    counter = 0
    for a in algorithms:
        if not x_selector(a):
            y_results = [ selector(g, lambda x: x.algorithm == a.algorithm) for g in grouped ]

            for status in STATUS:
                f = [e for e in y_results if e.status == status.value]
                if status.value == 'FINISHED':
                    plt.loglog([y_accessor(e) for e in x_results],
                               [y_accessor(e) for e in y_results],
                               SYMBOL_LIST[counter],
                               label=a.algorithm,
                               markersize=MARKER_SIZE
                              )
                else:
                    plt.loglog([y_accessor(e) for e in x_results],
                               [y_accessor(e) for e in y_results],
                               STATUS_SYMBOLS[status.value],
                               markersize=MARKER_SIZE
                              )
                    # plotting again for differentiation
                    plt.loglog([y_accessor(e) for e in x_results],
                               [y_accessor(e) for e in y_results],
                               SYMBOL_LIST[counter],
                               markersize=MARKER_SIZE*0.5
                              )
            counter = counter + 1



def print_csv( data, filename ):
    """Prints a list of iterables to file
    """
    try:
        with open(filename, 'w') as f:
            writer = csv.writer(f)
            writer.writerows(data)
    except Exception as e:
        print(e.message)



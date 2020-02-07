import csv
import json
from os import listdir
from os.path import isfile, join, splitext
import re
import sys
from testytools.result import Result as R

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
    errlist = [ '\t'.join((e.__dict__.get('slurm out'), e.errmsg)) for e in errors]
    print('\n'.join(errlist))

def group_designs( results, f=lambda x: x.design_name ):
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
            grouped[f(r)] = list(r)
        grouped[f(r)] = grouped[f(r)] + list(r)

    return [ tuple(e) for e in grouped.values() ]




def print_csv( data, filename ):
    """Prints a list of iterables to file
    """
    try:
        with open(filename, 'w') as f:
            writer = csv.writer(f)
            writer.writerows(data)
    except Exception as e:
        print(e.message)



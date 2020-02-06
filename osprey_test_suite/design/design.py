# -*- coding: utf-8 -*-
"""design.py

Module containing metyhods to run design algorithms in OSPREY

"""

##################################################
# Version
#   0.0
# Author
#   Graham Holt
##################################################

# TODO: Make the JSON file print out at the beginning and then be overwritten
# upon completion

##################################################
# Imports
##################################################
import argparse
from datetime import datetime
import errno
from enum import Enum
import fcntl
import imp
from itertools import izip_longest
import jpype
import json
from math import floor
import os
import osprey
import re
import socket
import sys
import time

##################################################
# Globals
##################################################
CFS_DIR = \
"/usr/project/dlab/Users/gth/projects/osprey_test_suite/updated_full_mut"
""" The directory in which confspace files are kept
"""

XTMP_DIR = "/usr/xtmp/gth/designScratch"
""" The directory in which to store large, unused files
"""

SCORE_PATTERN = re.compile(r"(^\S+)\s+in\s+\[(\S+)\s*,\s*(\S+)\]\s+\(log10\)")
""" The pattern for regexing K* scores
"""

HEADERS = ['cfs',
           'design name',
           'host',
           'pdb',
           'results',
           'runtime (s)',
           'epsilon',
           'memory (GB)',
           'status',
           'slurm out',
           'slurm err',
           'algorithm',
           'errmsg',
          ]
""" The field names for data that we are interested in collecting

cfs:    The full path to the design .cfs file
design name:    The name of the design
host:   The host on which the design was run
pdb:    The PDB code for the design
results:    A dictionary containing the sequences and their K* scores
epsilon:    The epsilon for the design
memory: The total RAM available to OSPREY
"""

DATETIME_FORMAT = "%Y-%m-%dT%H:%M:%S:%f"

# Here we store selected OSPREY Settings
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
HEAP_GB = 150
""" The Java heap size in GB"""

##################################################
# Classes / Functions
##################################################
class STATUS(Enum):
    """STATUS
    Contains the possible design status
    """
    ESTIMATING = "ESTIMATING"
    FINISHED = "FINISHED"
    ERROR = "ERROR"

def make_confspace (cfs_file, data):
    """make_confspace

    Generate conf spaces from input file

    Takes as input a cfs_file that contains chain definitions, flexible
    residues, and mutable residues.

    Returns a dictionary containing forcefield parameters and OSPREY
    confspace objects for:
        protein
        ligand
        complex
    """
    # Get the pdb code for this design
    data["pdb"] = cfs_file[:4]
    print("pdb: %s" % data["pdb"])

    # Load the design conformation space information
    data["cfs"] = os.path.join(CFS_DIR, cfs_file)
    confspace = imp.load_source('confspace',data["cfs"])
    print("CFS File: %s" % data["cfs"])

    ffparams = osprey.ForcefieldParams()
    mol = osprey.readPdb(confspace.mol)
    template_library = osprey.TemplateLibrary(ffparams.forcefld)

    # Make sure we don't have 3 strands (some cfs files do)
    assert len(confspace.strand_defs)==2

    # Define the protein strand
    protein = osprey.Strand(mol,
                            templateLib=template_library,
                            residues=confspace.strand_defs["strand0"]
                           )
    for resi, res_allowed in confspace.strand_flex["strand0"].iteritems():
        protein.flexibility[resi]\
                .setLibraryRotamers(osprey.WILD_TYPE, *res_allowed)\
                .setContinuous()\
                .addWildTypeRotamers()

    # Define the ligand strand
    ligand = osprey.Strand(mol,
                           templateLib=template_library,
                           residues=confspace.strand_defs["strand1"]
                          )
    for resi, res_allowed in confspace.strand_flex["strand1"].iteritems():
        ligand.flexibility[resi]\
                .setLibraryRotamers(osprey.WILD_TYPE,*res_allowed)\
                .setContinuous()\
                .addWildTypeRotamers()

    # Build spaces
    return {'protein' : osprey.ConfSpace(protein),
            'ligand' : osprey.ConfSpace(ligand),
            'complex': osprey.ConfSpace([protein,ligand]),
            'ffparams' : ffparams
           }

def get_emat_name(fileName):
    """get_emat_name

    Get a unique identifier for the energy matrices based on the filename
    """
    m = re.search('/(\\w+)_[^/]*res\\_(\\d+\\.\\d+E\\+\\d+).cfs',fileName)
    return m.group(1)+'_'+m.group(2)

def setup_design(conf_spaces, eps, num_seqs, use_shark, data):
    """setup_design

    Set up the classes required to run the design

    conf_spaces: A dictionary containing protein, ligand, and complex conf
    spaces along with force-field parameters
    num_seqs:   The number of top sequences we want (usually 5)

    Returns a BBKStar instance, ready to run
    """

    parallelism = osprey.Parallelism(cpuCores=20)

    # how should we compute energies of molecules?
    minimizingEcalc = osprey.EnergyCalculator(conf_spaces['complex'],
                                    conf_spaces['ffparams'],
                                    parallelism=parallelism,
                                   isMinimizing=True
                                   )

    # configure basic inputs for SHARKStar
    design = osprey.BBKStar(
        conf_spaces['protein'],
        conf_spaces['ligand'],
        conf_spaces['complex'],
        numBestSequences=num_seqs,
        epsilon=eps, # you proabably want something more precise in your real designs
        showPfuncProgress=True,
        maxSimultaneousMutations=9,
    )

    # configure SHARK*/BBK* inputs for each conf space
    configure_bbk(design, minimizingEcalc, "shark", data["design name"],
                  use_shark)

    return design

def configure_bbk(instance, minimizingEcalc, type_string, id_obj, useSHARKStar):
    """Configure the energy calculators and info for BBK* instance"""

    for info in instance.confSpaceInfos():
        # Compute reference energies
        eref = osprey.ReferenceEnergies(info.confSpace, minimizingEcalc)
        #Create a minimizing energy calculator
        info.confEcalcMinimized = osprey.ConfEnergyCalculator(
            info.confSpace,
            minimizingEcalc,
            referenceEnergies=eref)
 
        # BBK* needs rigid energies too
        rigidEcalc = osprey.SharedEnergyCalculator(
            minimizingEcalc,
            isMinimizing=False)
        rigidConfEcalc = osprey.ConfEnergyCalculatorCopy(
            info.confEcalcMinimized,
            rigidEcalc)
        info.confEcalcRigid = rigidConfEcalc

        # Specify the input for the partition functions. Providing the confUpperBoundcalc turns on SHARK*
        if useSHARKStar :
            info.pfuncFactory = osprey.PartitionFunctionFactory(
                info.confSpace,
                info.confEcalcMinimized,
                info.id,
                confUpperBoundcalc=rigidConfEcalc)
        else :
            info.pfuncFactory = osprey.PartitionFunctionFactory(info.confSpace, info.confEcalcMinimized, info.id)

        # Set cache pattern
        info.pfuncFactory.setCachePattern('%s/emat.%s.%s'
                                          % (XTMP_DIR, info.id, id_obj))
        print('Cache pattern: %s/emat.%s.%s' % (XTMP_DIR, info.id, id_obj))

        # compute the energy matrices
        info.ematMinimized = info.pfuncFactory.getOrMakeEmat(
            info.confEcalcMinimized,
            'minimized')

        info.ematRigid = info.pfuncFactory.getOrMakeEmat(
            info.confEcalcRigid,
            'rigid')

        # Updating energy matrix?
        info.ematCorrected =\
            osprey.c.ematrix.UpdatingEnergyMatrix(
                info.confSpace,
                info.ematMinimized)

        # how should confs be ordered and searched? (don't forget to capture emat by using a defaulted argument)
        def makeAStar_min(rcs, emat=info.ematMinimized):
            return osprey.AStarTraditional(emat, rcs, showProgress=True)

        info.confSearchFactoryMinimized =\
                osprey.KStar.ConfSearchFactory(makeAStar_min)

        def makeRigidAStar(rcs, emat=info.ematRigid):
            return osprey.AStarTraditional(emat, rcs, showProgress=True)

        info.confSearchFactoryRigid =\
                osprey.KStar.ConfSearchFactory(makeRigidAStar)

def parse_result(sequence):
    """parse_result

    Try to get design information out of an OSPREY sequence result
    """
    keys = ['sequence', 'kscore', 'lowerbound', 'upperbound']
    seqname = str(sequence.sequence)
    values = [seqname] +\
        list(re.search(SCORE_PATTERN, str(sequence.score)).groups())
    return dict(zip(keys,values))

def run(design, data):
    """run

    Runs the design instance
    """
    print("\nRunning SHARK*")
    design_start = datetime.now()
    data['start_time'] = design_start.strftime(DATETIME_FORMAT)
    try:
        scoredSequences = design.run()
        data['results'] = [ parse_result(seq) for seq in scoredSequences ]
        data["status"] = STATUS.FINISHED.value
        for seq in scoredSequences:
            print("%s: %s",(str(seq.sequence),str(seq.score)))
    except (jpype.JavaException, AttributeError) as e:
        data["status"] = STATUS.ERROR.value
        print str(e)
        print e.stacktrace()
        data['errmsg'] = str(e)

    design_end = datetime.now()
    data['end_time'] = design_end.strftime(DATETIME_FORMAT)
    data['runtime (s)'] = (design_end - design_start).total_seconds()
    print("SHARKStar Runtime (s): %d" % data['runtime (s)'])


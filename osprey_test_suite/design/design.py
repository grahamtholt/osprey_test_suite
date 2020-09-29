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

# TODO: Consolidate code. Currently setup_design, make_complex* and
# configure_bbk all have some overlap

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
#
DEFAULT_CFS_DIR = \
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

ALGO_LIST = ['SHARK', 'MARK', 'BBK']

BBK_BATCH_SIZE = 8
""" The number of conformations to compute before switching in BBK*. This has a
HUGE impact on performance, and that impact depends on how many cores you are using.
"""

# Here we store selected OSPREY Settings
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
HEAP_GB = 240
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

def make_confspace (cfs_file, data, cfs_dir=DEFAULT_CFS_DIR):
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
    data["cfs"] = os.path.join(cfs_dir, cfs_file)
    confspace = imp.load_source('confspace',data["cfs"])
    print("CFS File: %s" % data["cfs"])

    # Get the sequences if they exist
    try:
        data["sequences"] = confspace.sequences
    except AttributeError:
        data["sequences"] = [None];


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
    m = re.search('(\\w+)_[^/]*res\\_(\\d+\\.\\d+E\\+\\d+).cfs',fileName)
    return m.group(1)+'_'+m.group(2)

def setup_design(numcores, conf_spaces, eps, num_seqs, algo_index, data):
    """setup_design

    Set up the classes required to run the design

    conf_spaces: A dictionary containing protein, ligand, and complex conf
    spaces along with force-field parameters
    num_seqs:   The number of top sequences we want (usually 5)

    Returns a BBKStar instance, ready to run
    """

    parallelism = osprey.Parallelism(cpuCores=numcores)
    data['numCpus'] = numcores

    # how should we compute energies of molecules?
    minimizingEcalc = osprey.EnergyCalculator(conf_spaces['complex'],
                                    conf_spaces['ffparams'],
                                    parallelism=parallelism,
                                   isMinimizing=True
                                   )
    # record the seq tree file name
    data['seq_tree_name'] = data["design name"]+"seqTree.tree"

    # configure basic inputs for SHARKStar
    design = osprey.BBKStar(
        conf_spaces['protein'],
        conf_spaces['ligand'],
        conf_spaces['complex'],
        numBestSequences=num_seqs,
        epsilon=eps, # you proabably want something more precise in your real designs
        showPfuncProgress=True,
        maxSimultaneousMutations=9,
        numConfsPerBatch=BBK_BATCH_SIZE,
        maxNumConfsPerBatch=8,
        #printSequenceTree=data["seq_tree_name"]
    )

    # configure SHARK*/BBK* inputs for each conf space
    configure_bbk(design, minimizingEcalc, "shark", data["design name"],
                  algo_index)

    return design

def configure_bbk(instance, minimizingEcalc, type_string, id_obj, algo_index):
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
        if ALGO_LIST[algo_index] == 'SHARK':
            impt_ecalc = rigidConfEcalc
            choose_markstar=False
        elif ALGO_LIST[algo_index] == 'MARK':
            impt_ecalc = rigidConfEcalc
            choose_markstar=True
        else:
            impt_ecalc = None
            choose_markstar=False

        info.pfuncFactory = osprey.PartitionFunctionFactory(
            info.confSpace,
            info.confEcalcMinimized,
            info.id,
            confUpperBoundcalc=impt_ecalc,
            useMARK=choose_markstar
        )

        # Set cache pattern
        info.pfuncFactory.setCachePattern('%s/emat.%s.%s'
                                          % (XTMP_DIR,
                                             info.id,
                                             id_obj))
        print('Cache pattern: %s/emat.%s.%s'
              % (XTMP_DIR,
                 info.id,
                 id_obj))

        # compute the energy matrices
        info.ematMinimized = info.pfuncFactory.getOrMakeEmat(
            info.confEcalcMinimized,
            'minimized')

        info.ematRigid = info.pfuncFactory.getOrMakeEmat(
            info.confEcalcRigid,
            'rigid')

        # Updating energy matrix?
        #info.ematCorrected =\
            #osprey.c.ematrix.UpdatingEnergyMatrix(
                #info.confSpace,
                #info.ematMinimized,
                #info.confEcalcMinimized
            #)

        # how should confs be ordered and searched? (don't forget to capture emat by using a defaulted argument)
        def makeAStar_min(rcs, emat=info.ematMinimized):
            return osprey.AStarTraditional(emat, rcs, showProgress=True)

        info.confSearchFactoryMinimized =\
                osprey.KStar.ConfSearchFactory(makeAStar_min)

        def makeRigidAStar(rcs, emat=info.ematRigid):
            return osprey.AStarTraditional(emat, rcs, showProgress=True)

        info.confSearchFactoryRigid =\
                osprey.KStar.ConfSearchFactory(makeRigidAStar)

def make_complex_pfunc(numcores, conf_spaces, eps, num_seqs, algo_index, data):
    """Make a single partition function for the flexible part of the complex
    confspace. This will be useful for comparing the core of SHARK* to the core
    of MARK*
    """
    # make the flexible, complex confspace
    flex_complex_space = conf_spaces['complex'].makeFlexibleCopy()
    # make the factory
    pfunc_factory = make_pfunc_factory(flex_complex_space,
                                     conf_spaces['ffparams'],
                                     numcores,
                                     eps,
                                     algo_index,
                                     "flexible_complex_confspace",
                                     data)
    # make the pfunc for the wild-type sequence
    return make_pfunc_for_sequence(flex_complex_space, pfunc_factory, eps, data,
                                   seq_list=None)

def make_pfunc_factory(conf_space, ffparams, numcores, eps, algo_index, emat_cache_pattern ,data):
    """make_pfunc_factory

    Return a PartitionFunctionFactory object for the input confspace, which we can use to make
    partition functions for various sequences.
    """
    parallelism = osprey.Parallelism(cpuCores=numcores)
    data['numCpus'] = numcores

    # how should we compute energies of molecules?
    minimizingEcalc = osprey.EnergyCalculator(conf_space,
                                    ffparams,
                                    parallelism=parallelism,
                                   isMinimizing=True
                                   )
    # Compute reference energies
    eref = osprey.ReferenceEnergies(conf_space, minimizingEcalc)

    #Create a minimizing energy calculator
    confEcalcMinimized = osprey.ConfEnergyCalculator(
        conf_space,
        minimizingEcalc,
        referenceEnergies=eref)

    # we need rigid energies too for many algorithms
    rigidEcalc = osprey.SharedEnergyCalculator(
        minimizingEcalc,
        isMinimizing=False)
    rigidConfEcalc = osprey.ConfEnergyCalculatorCopy(
        confEcalcMinimized,
        rigidEcalc)
    confEcalcRigid = rigidConfEcalc

    # Specify the type of partitionFunction
    if ALGO_LIST[algo_index] == 'SHARK': # using SHARK*
        impt_ecalc = rigidConfEcalc
        choose_markstar=False
    elif ALGO_LIST[algo_index] == 'MARK': # using MARK*
        impt_ecalc = rigidConfEcalc
        choose_markstar=True
    else: # using Gradient descent pfunc
        impt_ecalc = None
        choose_markstar=False

    pfuncFactory = osprey.PartitionFunctionFactory(
        conf_space,
        confEcalcMinimized,
        emat_cache_pattern,
        confUpperBoundcalc=impt_ecalc,
        useMARK=choose_markstar
    )
    # Set cache pattern
    pfuncFactory.setCachePattern('%s/emat.%s.%s'
                                      % (XTMP_DIR,
                                         emat_cache_pattern,
                                         data["design name"]))
    print('Cache pattern: %s/emat.%s.%s'
          % (XTMP_DIR,
             emat_cache_pattern,
             data["design name"]))

    return pfuncFactory

def make_pfunc_for_sequence(conf_space, pfunc_factory, epsilon, data, seq_list=None):
    """make_pfunc_for_sequence

    Make a partition function for the given sequence and conf space from
        the given pfunc_factory.

    Args:
        conf_space: SimpleConfSpace
        pfunc_factory: PartitionFunctionFactory
        epsilon: The target approximation error
        data: A dictionary for recording test data -- side-effects
        seq_list: List of 3-letter residue code strings, or None.
            Residue code strings must be in the same order as positions in
            conf_space
    Returns:
        Object inheriting the PartitionFunction interface
    """
    # make the sequence object
    if seq_list is None:
        sequence = conf_space.makeWildTypeSequence()
    else:
        # We have to use this datastructure to match lists
        myList = jpype.java.util.ArrayList()
        for s in seq_list:
            myList.add(s)
        sequence = conf_space.seqSpace.makeSequence(myList)
    # make the rcs
    rcs = sequence.makeRCs(conf_space)
    data['numconfs'] = rcs.getNumConformations().toString() # record the size
    # make the pfunc and return it
    return pfunc_factory.makePartitionFunctionFor(rcs,
                                                  rcs.getNumConformations(),
                                                  epsilon,
                                                  sequence
                                                 )

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
    print("\nRunning Algorithm")
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
    print("Algorithm Runtime (s): %d" % data['runtime (s)'])
    if design.complexSHARK is not None and design.proteinSHARK is not None and\
        design.ligandSHARK is not None:
        data['precompute_complex_time'] =\
            design.complexSHARK.precomputedFlexComputeTime +\
                design.proteinSHARK.precomputedFlexComputeTime +\
                design.ligandSHARK.precomputedFlexComputeTime
        data['make_pfunc_time'] =\
            design.complexSHARK.pfuncCreationTime +\
                design.proteinSHARK.pfuncCreationTime +\
                design.ligandSHARK.pfuncCreationTime
    else:
        data['precompute_complex_time'] = 0
        data['make_pfunc_time'] = 0


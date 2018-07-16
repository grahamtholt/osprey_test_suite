#!/usr/bin/env python
# -*- coding: <utf-8> -*-
"""Generate a set of reasonable OSPREY designs from pdb

"""

from pymol import cmd, CmdException, stored
import os
import itertools

import atom
import residue
import interface
import gen_shell
import design

__author__ = "Graham Holt"
__version__ = "0.0.0"
__status__ = "Prototype"

def test_gen(pdb):
    """Generate all reasonable designs from pdb

    Args:
        pdb (string): The filepath to a pdb file

    Raises:
        CmdException: Error in a PyMol command

    """

    # Load pymol file
    #try:
       #cmd.load(pdb)
    #except:
        #raise CmdException

    # Generate all interfaces
    #interface.interface()
    # Note - this command seems to work well for picking only
    # sidechain-sidechain contacts
    interface.interface("","",3,1)
    try:
        inter_list = cmd.get_names("selections", all)
        print inter_list
    except:
        raise CmdException
        inter_list = [ e for e in inter_list if e[0:2]=="int" ]

    # Find potential mutable residues
    for inter in inter_list:
        mut_list.append(*gen_mut(inter))

    # For each mutable residue set, generate a flexible shell and print design
    for mut in mut_list:
        #Select muts
        #gen_shell(1)
        #design()
        pass


def gen_mut(sele_name):
    """Select combinations of mutable residues

    Selects all possible combinations of good mutable residues within an
    interface.

    Args:
        sele_name (string): The name of a PyMol selection that defines a PPI

    Returns:
        A list of lists of Residue objects that are good design choices
    """
    mut_list = []
    res_list = residue.Residue.sele_to_res(sele_name)

    # Add all single mutations
    mut_list.append(*[ [res] for res in res_list ])

    # Construct basis for close pairs
    close_pair_basis = [ e for e in itertools.combinations(mut_list, 2)
                        if e[0].min_sc_sc_distance(e[1]) < 3 ]
    print close_pair_basis

# Make this runnable as a command from pymol
cmd.extend( "test_gen", test_gen )

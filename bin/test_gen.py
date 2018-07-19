#!/usr/bin/env python
# -*- coding: <utf-8> -*-
"""Generate a set of reasonable OSPREY designs from pdb

"""

from pymol import cmd, CmdException, stored
import pymol
import os
import itertools
import argparse

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
    try:
       cmd.load(pdb)
    except:
        raise CmdException

    # Generate all interfaces
    #interface.interface()
    # Note - this command seems to work well for picking only
    # sidechain-sidechain contacts
    interface.interface(None,None,3,1)
    try:
        inter_list = cmd.get_names("selections", 0, "all")
    except:
        raise CmdException

    inter_list = [ e for e in inter_list if e[0:3]=="int" ]

    # Find potential mutable residues
    mut_list = []
    for inter in inter_list:
        chains = cmd.get_chains(inter)
        inter_1 = inter+' and chain '+chains[0]
        inter_2 = inter+' and chain '+chains[1]
        print inter_1
        mut_list.extend(gen_mut(inter_1))
        print inter_2
        mut_list.extend(gen_mut(inter_2))


    # For each mutable residue set, generate a flexible shell and print design
    counter = 0
    for mut in mut_list:
        #Select muts
        obj_name = cmd.get_names('objects',0,'all')[0]
        resi_list = [ str(res.res_seq)+res.i_code for res in mut ]
        chain = next(iter(mut)).chain_id
        selection = "br. chain "+chain+" and resi "+'+'.join(resi_list)
        cmd.select("mut", selection)

        #gen_shell
        gen_shell.gen_shell(1,0,"mut",5)

        #design()
        design.design("mut","flex",('_'.join([obj_name, chain,
                                              "confsize", str(counter)]))+".cfs")
        counter = counter+1

    gen_scenes(mut_list)

def gen_mut(sele_name):
    """Select combinations of mutable residues

    Selects all possible combinations of good mutable residues within an
    interface.

    Args:
        sele_name (string): The name of a PyMol selection that defines a PPI

    Returns:
        A set of sets of Residue objects that are good design choices
    """
    mut_list = set()
    res_list = residue.Residue.sele_to_res(sele_name)

    # Add all single mutations
    mut_list.update(frozenset([res]) for res in res_list)

    # Construct close pairs
    close_pairs = set( frozenset(e) for e in itertools.combinations(
        [ res for res in res_list ], 2)
                        if e[0].min_sc_sc_distance(e[1]) < 3 )

    # Add close pairs to muts
    mut_list.update(close_pairs)

    # Iteratively union close pairs until cannot any longer
    new_union = True
    old_set = close_pairs
    while new_union == True:
        new_union = False
        new_mut_set = set()
        for e in old_set:
            for pair in close_pairs:
                if len(e.intersection(pair)) == 1:
                    new_mut_set.add(e.union(pair))
                    new_union = True
        mut_list.update(new_mut_set)
        old_set = new_mut_set

    for mut in mut_list:
        print ', '.join(str(e) for e in mut)
    return mut_list

def is_globular(mut_set, dist):
    max_dist = 0
    for pair in itertools.combinations(iter(mut_set), 2):
        if pair[0].max_sc_sc_distance(pair[1]) < dist:
            return False

    return True

def gen_scenes(mut_list):

    # create scenes for each mut
    for index, mut in enumerate(mut_list):
        obj_name = cmd.get_names('objects',0,'all')[0]

    # modify individual states
        resi_list = [ str(res.res_seq)+res.i_code for res in mut ]
        chain = next(iter(mut)).chain_id
        selection = "br. chain "+chain+" and resi "+'+'.join(resi_list)
        cmd.create(obj_name, selection, 1, index+2)

        cmd.set("stick_color","red", "all", state=index)
        cmd.show("sticks",selection)

    # generate a reference
    cmd.create("reference", "all", 1, 1)

# Make this runnable as a command from pymol
cmd.extend( "test_gen", test_gen )

def main(args):
    pymol.pymol_argv = ['pymol','-qcK']
    pymol.finish_launching()

    test_gen(args.pdb)

##################################################
# Conditional main
##################################################

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='This description is shown when -h or --help are passed as arguments.')

    parser.add_argument('pdb',
            type = str,
            default = None,
            help = 'The pdb file from which to generate mutations')

    args = parser.parse_args()

    main(args)

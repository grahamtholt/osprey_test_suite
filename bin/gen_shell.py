#!/usr/bin/env python
"""Generate flexible shell around mutable residues

"""

import re
from pymol import cmd, CmdException, stored

import atom, residue

__author__ = "Graham Holt"
__version__ = "0.0.0"
__status__ = "Prototype"

def gen_shell(prune_sele=0, aggr=0, mut="mut", r=4):
    """Select a flexible shell around a mutable selection.

    @param prune_sele: Do we change the selection? If not, we just color AAs
    @param aggr: Integer setting for how aggresively we prune residues
    @param mut: The PyMol selection defining the mutable residues
    @param r: Radius cutoff for flexible residues in angstroms

    """

    # Define the flexible shell name
    if mut == "mut":
        flex_name = "flex"
    else:
        flex_name = mut+"_flex"

    # Select the full possible flexible shell
    selection = "br. (all and not "+mut+") within "+str(r)+" of "+mut
    try:
        flex_sele = cmd.select(flex_name, selection)
    except:
        raise CmdException

    mut_list = res_from_sele(mut)
    flex_list = res_from_sele(flex_name)


    prune_gly_pro(flex_name, mut_list, flex_list, aggr, r, prune_sele)
    prune_distance(flex_name, mut_list, flex_list, aggr, r, prune_sele)

def res_from_sele(sele_name):
    """Generate a list of residue objects from a PyMol selection

    @param sele_name: PyMol selection name
    """

    # Get the list of all atoms in the selection
    stored.atoms = []
    try:
        cmd.iterate(sele_name,
                    "stored.atoms.append([type, ID, name, alt, resn, chain,\
                    resi,\"\",(0,0,0), q, b, elem, \"\"])")
        # Get coordinates
        # NOTE: PYMOL MIGHT HAVE ROUNDING ERRORS FOR COORDINATES!
        coords_list = cmd.get_coords(sele_name, 1)
    except:
        raise CmdException

    # Write coordinates, extract icode, make atom lists, reslist
    res_list = []
    for e, coord in zip(stored.atoms, coords_list):
        e[8] = tuple(coord)
        m = re.match(r"(\d+)(\w?)", e[6])
        e[6] = int(m.group(1))
        e[7] = m.group(2)

        #make residues
        a = atom.Atom(*e)
        if not res_list or not res_list[-1].can_contain(a):
            res_list.append(residue.Residue([a]))
        else:
            res_list[-1].add_atom(a)

    return res_list

def prune_distance(flex_sele_name, mut_list, flex_list, aggr, r, prune_sele):
    """Remove residues from the flexible selection based on criterion

    This method either removes or colors residues in the flexible shell that
    are evaluated by the given criterion function as false for ALL in mut_list.

    @param flex_sele_name: The Pymol selection from which to prune
    @param mut_list: A stored list of mutable Residue objects
    @param flex_list: A stored list of flexible Residue objects
    @param aggr: Integer setting for how aggresively we prune residues
    @param r: Radius cutoff for flexible residues in angstroms
    @param prune_sele: Do we change the selection? If not, we just color AAs

    """

    for flex in flex_list:
        to_prune = True
        for mut in mut_list:
            if flex.min_sc_distance(mut) < (r-aggr*0.2):
                to_prune = False
        if to_prune:
            do_pruning(flex_sele_name, flex, prune_sele)

def prune_gly_pro(flex_sele_name, mut_list, flex_list, aggr, r, prune_sele):
    """Remove residues from the flexible selection based on residue type

    This method either removes or colors residues in the flexible shell that
    are glycines or prolines.

    @param flex_sele_name: The Pymol selection from which to prune
    @param mut_list: A stored list of mutable Residue objects
    @param flex_list: A stored list of flexible Residue objects
    @param aggr: Integer setting for how aggresively we prune residues
    @param r: Radius cutoff for flexible residues in angstroms
    @param prune_sele: Do we change the selection? If not, we just color AAs

    """

    for flex in flex_list:
        if flex.res_name == "GLY" or flex.res_name == "PRO":
            do_pruning(flex_sele_name, flex, prune_sele)

def do_pruning(sele_name, res, prune_sele):
    """Remove or color residue in PyMol selection

    @param sele_name
    @param res
    """

    if prune_sele:
        selection = ( sele_name+" and not (resi "+
                     str(res.res_seq)+res.i_code+
                     " and chain "+res.chain_id +")")
        try:
            cmd.select(sele_name, selection)
        except:
            raise CmdException
    else:
        selection = ( sele_name+" and (name C*) and (resi "+
                     str(res.res_seq)+res.i_code+
                     " and chain "+res.chain_id +")")
        try:
            cmd.color("gray50",selection)
        except:
            raise CmdException

# Make this runnable as a command from pymol
cmd.extend( "gen_shell", gen_shell )
#!/usr/bin/env python
"""Select protein-protein interface in PyMol

This module creates PyMol selections from the interface between
1) Two given selections, or 2) All polymer chains in the pdb.
"""

from pymol import cmd, CmdException, stored

__author__ = "Graham Holt"
__version__ = "0.0.0"
__status__ = "Prototype"

def select_flex(prune_sele=0, aggr=0, mut="mut", r=4):
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
    selection = "br. (all and not "+mut+") within "+r+" of "+mut
    try:
        flex_sele = cmd.select(flex_name, selection)
    except:
        raise CmdException

    # Get the residues in the flexible selection
    stored.flex = []
    try:
        cmd.iterate(flex_sele+" and name ca",
                    "stored.flex.append(\"\"+chain+resi)")
    except:
        raise CmdException

def prune_distance(mut_list, flex_list, prune_sele, aggr, r):
    """Remove residues from the flexible selection based on distance

    This method either removes or colors residues in the flexible shell that
    have no side chain atoms within the specified radius of any mutable atom.
    Scales with aggressiveness, up to r-1.

    @param mut_list: A stored list of mutable residues
    @param flex_list: A stored list of flexible residues
    @param prune_sele: Do we change the selection? If not, we just color AAs
    @param aggr: Integer setting for how aggresively we prune residues
    @param r: Radius cutoff for flexible residues in angstroms

    """
    pass


#!/usr/bin/env python
"""Print out a design config file from PyMol

This module reads one mutable and one flexible PyMol selection, and converts
them into variables convenient for defining an OSPREY 3.0 configuration space.
"""

import os
from pymol import cmd, CmdException, stored

__author__ = "Graham Holt"
__version__ = "0.0.0"
__status__ = "Production"

STRAND_NAME = "strand"
"""
The prefix to write to output for strand definitions
"""

DEFAULT_MUTATIONS_ALL = [
'ALA', 'VAL', 'LEU', 'ILE', 'PHE', 'TYR', 'TRP', 'CYS', 'MET', 'SER',
'THR', 'LYS', 'ARG', 'HIP', 'HIE', 'ASP', 'GLU', 'ASN', 'GLN', 'GLY'
]
"""list: List of amino acid templates to allow by default for mutable residues.

"""

HEADER = """#The following information is sufficient to generate an OSPREY 3.0
#configuration space"""
"""
The header to write out to the config file
"""

def design(mut="mut", flex="flex", out="design.cfg"):
    """
    Print out the OSPREY config information for a design to file

    @param mut: PyMol selection name for the mutable design residues
    @param flex: PyMol selection name for the flexible design residues
    @param out: System path to output file
    """
    # Find the input pdb file
    pdb_file_name = find_pdb_file(mut, flex)

    # Collect the chain ID and residue number for mutable and flexible res
    stored.mut_list = []
    stored.flex_list = []

    cmd.iterate(mut+" and name ca", "stored.mut_list.append(\"\"+chain+resi)")
    cmd.iterate(flex+" and name ca", "stored.flex_list.append(\"\"+chain+resi)")

    if set(stored.mut_list).intersection(set(stored.flex_list)):
        print ("You forgot to remove the mutable residues from the"+
            "flexible ones!")
        raise CmdException

    # Collect the chains required for design
    stored.chain_list = []

    cmd.iterate(mut+" and name ca", "stored.chain_list.append(\"\"+chain)")
    cmd.iterate(flex+" and name ca", "stored.chain_list.append(\"\"+chain)")
    chains = set( stored.chain_list )

    # For each chain, get the beginning and end residues
    strand_defs = {}
    strand_defs_all = {}
    for counter, chain in enumerate(chains):
        # First, get all residues in chain
        stored.all_res_list = []
        cmd.iterate("chain "+chain+" and name ca",
                    "stored.all_res_list.append(\"\"+chain+resi)")
        # Store only the first and last
        strand_defs[STRAND_NAME+str(counter)] = [ stored.all_res_list[0],
                                                 stored.all_res_list[-1] ]
        # Also store the whole thing for further use
        strand_defs_all[STRAND_NAME+str(counter)] = stored.all_res_list

    # Define strand flexibilities / mutabilities
    strand_flex_all = {}
    for strand in strand_defs_all:
        this_strand_flex = {}
        # For mutable residues in strand, add all mutations by default
        for el in set(strand_defs_all[strand]).intersection(set(stored.mut_list)):
            this_strand_flex[el] = DEFAULT_MUTATIONS_ALL
        # For flexible residues in strand, add all mutations by default
        for el in set(strand_defs_all[strand]).intersection(set(stored.flex_list)):
            this_strand_flex[el] = []
        strand_flex_all[strand] = this_strand_flex

    # Print out variables
    with open(out, 'w') as f:
        f.write(HEADER+"\n")
        f.write("mol = \""+str(pdb_file_name)+"\"\n")
        for strand in strand_defs:
            f.write(strand+" = "+str(strand_defs[strand])+"\n")
        for strand in strand_flex_all:
            f.write(strand+"_flex = "+str(strand_flex_all[strand])+"\n")

def find_pdb_file(mut, flex):
    """
    Find the pdb file of a pymol selection

    This function checks the current working directory for a pdb file with the
    same name as the selection model. This will ONLY work if the file is in the
    cwd.

    @param mut: PyMol selection name for the mutable design residues
    @param flex: PyMol selection name for the flexible design residues
    """
    mut_obj_list = cmd.get_names("objects", 0, mut)
    flex_obj_list = cmd.get_names("objects", 0, flex)

    assert len(mut_obj_list) == 1
    assert len(flex_obj_list) == 1

    mut_obj = mut_obj_list[0]
    flex_obj = flex_obj_list[0]

    assert mut_obj == flex_obj

    pdb_rel = mut_obj+".pdb"

    if os.path.isfile(pdb_rel):
        # Return the absolute path to the pdb
        return os.path.abspath(pdb_rel)
        # Return the relative path to the pdb
        # return pdb_rel
    else:
        print "ERROR! Could not find pdb file. Please set manually"
        return "UNSET"


# Make this executable as a command from pymol
cmd.extend( "design", design )

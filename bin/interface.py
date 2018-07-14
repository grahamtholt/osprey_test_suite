#!/usr/bin/env python
"""Select protein-protein interface in PyMol

This module creates PyMol selections from the interface between
1) Two given selections, or 2) All polymer chains in the pdb.
"""

from itertools import combinations
from pymol import cmd, CmdException, stored

__author__ = "Graham Holt"
__version__ = "0.0.0"
__status__ = "Prototype"

def interface( sele1="", sele2="", d=4):
    """Select and name the interfaces between selections.

    Defaults to combinations between all chains.

    @param sele1: First PyMol selection
    @param sele2: Second PyMol selection
    @param d: Depth of interface

    """
    d = str(d)
    stored.chains = []

    if sele1 == "" or sele2 == "":
        chain_names = cmd.get_chains('polymer')
        stored.chains = [ "chain "+s for s in chain_names ]

    else:
        stored.chains.append(sele1)
        stored.chains.append(sele2)

    for a, b in combinations(stored.chains, 2):
        selection = ( "(br. "+a+" within "+d+" of "+b+")"+
                     "or (br. "+b+" within "+d+" of "+a+")"
                    )
        # Make nicer names for the default case
        a_name = "".join(a.split()).replace("chain","")
        b_name = "".join(b.split()).replace("chain","")
        sel_name = "int_"+a_name+b_name
        try:
            sele = cmd.select(sel_name, selection)
        except:
            print "Error in selection: "+selection
            raise CmdException

        # Remove empty selections
        stored.test = []
        try:
            cmd.iterate(sel_name+" and name ca",
                        "stored.test.append(\"\"+resi)")
        except:
            print "Cannot iterate through "+sel_name
            raise CmdException
        if not stored.test:
            try:
                cmd.delete(sel_name)
            except:
                print "Cannot delete "+sel_name
                raise CmdException

# Make this runnable as a command from pymol
cmd.extend( "interface", interface )

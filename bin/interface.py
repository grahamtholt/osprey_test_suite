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

def interface( sele1=None, sele2=None, d=4, sc_only=0):
    """Select and name the interfaces between selections.

    Defaults to combinations between all chains.

    Args:
        sele1 (str): First PyMol selection
        sele2 (str): Second PyMol selection
        d (int): Depth of interface
        sc_only (int): Filter by side-chain contacts only, instead of all atoms

    """
    d = str(d)
    stored.chains = []

    if "\"\"" in [sele1, sele2] or None in [sele1, sele2]:
        chain_names = cmd.get_chains('polymer')
        stored.chains = [ "chain "+s for s in chain_names ]

    else:
        stored.chains.append(sele1)
        stored.chains.append(sele2)

    for a, b in combinations(stored.chains, 2):
        if int(sc_only) == 1:
            selection = ( "(br. ("+a+" and not name c+ca+o+n+h+ha) "+
                         "within "+d+" of ("+b+" and not name c+ca+o+n+h+ha)) "+
                         " or (br. ("+b+" and not name c+ca+o+n+h+ha) "+
                         "within "+d+" of ("+a+" and not name c+ca+o+n+h+ha)) "
                        )
        else:
            selection = ( "(br. "+a+" within "+d+" of "+b+")"+
                         " or (br. "+b+" within "+d+" of "+a+")"
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

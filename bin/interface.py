from pymol import cmd, CmdException, stored
from itertools import combinations

def interface( sele1="", sele2="", d=4):
    """
Select and name the interfaces between selections.
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
        try:
            sele = cmd.select( "int_"+a_name+b_name, selection)
        except:
            print "Error in selection: "+selection
            raise CmdException

# Make this runnable as a command from pymol
cmd.extend( "interface", interface )

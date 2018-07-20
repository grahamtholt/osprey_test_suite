"""Methods to manipulate PyMol sessions and PDB files

Todo:
    * Improve gen_shell residue pruning. This is the most obvious way to
    improve things
"""
from pymol import cmd, CmdException, stored

import re
import os
from itertools import combinations
from itertools import chain as iterchain
import argparse
import imp

import atom
import residue

__author__ = "Graham Holt"
__version__ = "0.0.0"
__status__ = "Prototype"

STRAND_NAME = "strand"
"""
The prefix to write to output for strand definitions
"""

DEFAULT_MUTATIONS_ALL = [
'ALA', 'VAL', 'LEU', 'ILE', 'PHE', 'TYR', 'TRP', 'CYS', 'MET', 'SER',
'THR', 'LYS', 'ARG', 'HIP', 'HIE', 'HID', 'ASP', 'GLU', 'ASN', 'GLN', 'GLY'
]
"""list: List of amino acid templates to allow by default for mutable residues.

"""

HEADER = """#The following information is sufficient to generate an OSPREY 3.0
#configuration space"""
"""
The header to write out to the config file
"""

def interface(sele1=None, sele2=None, d=4, sc_only=0):
    """Select and name the interfaces between selections.

    Defaults to combinations between all chains.

    Args:
        sele1 (str): First PyMol selection
        sele2 (str): Second PyMol selection
        d (int): Depth of interface
        sc_only (int): Filter by side-chain contacts only, instead of all atoms

    Returns:
        Nothing

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

def gen_shell(prune_sele=0, aggr=0, mut="mut", r=4):
    """Select a flexible shell around a mutable selection.

    Args:
        prune_sele (int): Do we change the selection? If not, we just color AAs
        aggr (int): how aggressively we prune residues
        mut (str): The PyMol selection defining the mutable residues
        r i(int): Radius cutoff for flexible residues in angstroms

    Returns:
        Nothing

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

    # modify view
    try:
        cmd.hide("everything", "all")
        cmd.show("sticks", mut)
        cmd.show("lines", flex_name)
        cmd.hide("everything","hydrogens")
        cmd.orient(flex_name)
    except:
        raise CmdException

    mut_list = residue.Residue.sele_to_res(mut)
    flex_list = residue.Residue.sele_to_res(flex_name)


    prune_gly_pro(flex_name, mut_list, flex_list, aggr, r, prune_sele)
    prune_distance(flex_name, mut_list, flex_list, aggr, r, prune_sele)
    prune_angle(flex_name, mut_list, flex_list, aggr, r, prune_sele)

def prune_distance(flex_sele_name, mut_list, flex_list, aggr, r, prune_sele):
    """Remove residues from the flexible selection based on distance

    This method either removes or colors residues in the flexible shell that
    are not within the given distance for ALL in mut_list.

    Args:
        flex_sele_name (str): The Pymol selection from which to prune
        mut_list (list<Residue>): A stored list of mutable Residue objects
        flex_list (list<Residue>): A stored list of flexible Residue objects
        aggr (int): How aggresively we prune residues
        r (int): Radius cutoff for flexible residues in angstroms
        prune_sele (int): Do we change the selection? If not, we just color AAs

    Returns:
        Nothing

    """

    for flex in flex_list:
        to_prune = True
        for mut in mut_list:
            #TODO: Make a decision on this. MAybe sc_any for bb flex?
            #if flex.min_sc_any_distance(mut) < (r-aggr*0.2):
            if flex.min_sc_sc_distance(mut) < (r-aggr*0.2):
                to_prune = False
        if to_prune:
            do_pruning(flex_sele_name, flex, prune_sele)

def prune_angle(flex_sele_name, mut_list, flex_list, aggr, r, prune_sele):
    """Prune by angle between calpha-cbeta vectors

    This method either removes or colors residues in the flexible shell that
    do not point sufficiently toward the mutable residue

    Args:
        flex_sele_name (str): The Pymol selection from which to prune
        mut_list (list<Residue>): A stored list of mutable Residue objects
        flex_list (list<Residue>): A stored list of flexible Residue objects
        aggr (int): How aggresively we prune residues
        r (int): Radius cutoff for flexible residues in angstroms
        prune_sele (int): Do we change the selection? If not, we just color AAs

    """

    for flex in flex_list:
        to_prune = True
        for mut in mut_list:
            if not flex.angles_diverge(mut, 90):
                to_prune = False
        if to_prune:
            do_pruning(flex_sele_name, flex, prune_sele)

def prune_gly_pro(flex_sele_name, mut_list, flex_list, aggr, r, prune_sele):
    """Remove residues from the flexible selection based on residue type

    This method either removes or colors residues in the flexible shell that
    are glycines or prolines.

    Args:
        flex_sele_name (str): The Pymol selection from which to prune
        mut_list (list<Residue>): A stored list of mutable Residue objects
        flex_list (list<Residue>): A stored list of flexible Residue objects
        aggr (int): How aggresively we prune residues
        r (int): Radius cutoff for flexible residues in angstroms
        prune_sele (int): Do we change the selection? If not, we just color AAs

    """

    for flex in flex_list:
        if flex.res_name == "GLY" or flex.res_name == "PRO":
            do_pruning(flex_sele_name, flex, prune_sele)

def do_pruning(sele_name, res, prune_sele):
    """Remove or color residue in PyMol selection

    Args:
        sele_name (str): A PyMol selection name
        res (Residue): A residue object to prune
        prune_sele (int): Do we change the selection? If not, we just color AAs

    Returns:
        Nothing

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

def test_gen():
    """Generate all reasonable designs from the current pdb

    Args:

    Raises:
        CmdException: Error in a PyMol command

    """

    # Generate all interfaces
    # Note - this command seems to work well for picking only
    # sidechain-sidechain contacts
    interface(None,None,3,1)
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


    print 'Generated '+str(len(mut_list))+' sets of mutations'
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
        gen_shell(1,0,"mut",5)

        #print output
        print_design("mut","flex",('_'.join([obj_name, chain,
                                              "confsize", str(counter),
                                             str(len(mut))+'muts']))+".cfs")
        counter = counter+1

    #optional, but slow
    #gen_scenes(mut_list)

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
    close_pairs = set( frozenset(e) for e in combinations(
        [ res for res in res_list ], 2)
                        if e[0].min_sc_sc_distance(e[1]) < 3 )

    # Add close pairs to muts
    mut_list.update(close_pairs)

    # Iteratively union close pairs until cannot any longer
    # Or, if we reach a maximum of 6 mutable residues
    max_num_mut = 6
    mut_counter = 2
    new_union = True
    old_set = close_pairs
    while new_union == True and mut_counter < max_num_mut:
        new_union = False
        new_mut_set = set()
        for e in old_set:
            for pair in close_pairs:
                if len(e.intersection(pair)) == 1:
                    new_mut_set.add(e.union(pair))
                    new_union = True
        mut_list.update(new_mut_set)
        old_set = new_mut_set
        mut_counter = mut_counter + 1

    return mut_list

def is_globular(mut_set, dist):
    max_dist = 0
    for pair in combinations(iter(mut_set), 2):
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

def print_design(mut="mut", flex="flex", out="design.cfs"):
    """
    Print out the OSPREY config information for a design to file

    Args:
        mut (str): PyMol selection name for the mutable design residues
        flex (str): PyMol selection name for the flexible design residues
        out (str): System path to output file
    """
    # Find the input pdb file
    pdb_file_name = find_pdb_file(mut, flex)

    # Collect residue lists from selections
    mut_list = residue.Residue.sele_to_res(mut)
    flex_list = residue.Residue.sele_to_res(flex)

    # Make sure we don't have the mutable residues in the flexible selection
    if set(mut_list).intersection(set(flex_list)):
        print ("You forgot to remove the mutable residues from the"+
            "flexible ones!")
        raise CmdException

    # Collect the chains required for design
    # Note that the "chain" function here comes from itertools
    chains = set([ res.chain_id for res in iterchain(mut_list, flex_list) ])

    # TODO: Consider changing the format to be easier to read and write...
    # For each chain, process muts and flex
    strand_defs = {}
    strand_flex_all = {}
    for counter, chain in enumerate(chains):
        strand_key = STRAND_NAME+str(counter)

        # First, get all residues in chain
        stored.all_res_list = []
        cmd.iterate("chain "+chain+" and name ca",
                    "stored.all_res_list.append(\"\"+chain+resi)")
        # Store only the first and last
        strand_defs[strand_key] = [ stored.all_res_list[0],
                                                 stored.all_res_list[-1] ]

        # Define strand flexibilities / mutabilities
        this_strand_flex = {}
        # For mutable residues in strand, add all mutations by default
        chain_muts = (res for res in mut_list if res.chain_id == chain)
        for res in chain_muts:
            this_strand_flex[res.chain_id+str(res.res_seq)+res.i_code] = \
                    DEFAULT_MUTATIONS_ALL
        # For flexible residues in strand, add all mutations by default
        chain_muts = (res for res in flex_list if res.chain_id == chain)
        for res in chain_muts:
            this_strand_flex[res.chain_id+str(res.res_seq)+res.i_code] = \
                    [res.res_name]

        strand_flex_all[strand_key] = this_strand_flex

    # Print out variables
    with open(out, 'w') as f:
        f.write(HEADER+"\n")
        f.write("mol = \""+str(pdb_file_name)+"\"\n")
        f.write("strand_defs = "+str(strand_defs)+"\n")
        f.write("strand_flex = "+str(strand_flex_all)+"\n")
        #for strand in strand_defs:
            #f.write(strand+" = "+str(strand_defs[strand])+"\n")
        #for strand in strand_flex_all:
            #f.write(strand+"_flex = "+str(strand_flex_all[strand])+"\n")

def find_pdb_file(mut, flex):
    """
    Find the pdb file of a pymol selection

    This function checks the current working directory for a pdb file with the
    same name as the selection model. This will ONLY work if the file is in the
    cwd.

    Args:
        mut (str): PyMol selection name for the mutable design residues
        flex(str): PyMol selection name for the flexible design residues
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

def load_cfs(cfs_file):
    """Load a configuration space file as a visual PyMol design

    Args:
        cfs_file (str): The cfs file that contains the design confspace
    """

    # import confspace file
    conf_space = imp.load_source('conf_space', cfs_file)
    # get flexible and mutable residues in nice format
    mut_dict = dict()
    flex_dict = dict()
    for d in conf_space.strand_flex.values():
        mut_list = []
        flex_list = []
        chain = next(iter(d.keys()))[:1]
        for res, mut_allowed in d.iteritems():
            # Check to make sure all chains match
            this_chain = res[:1]
            assert chain == this_chain

            resi = res[1:]
            if len(mut_allowed) > 1:
                mut_list.append(resi)
            else:
                flex_list.append(resi)

        # Assign resis to dicts
        mut_dict[chain] = mut_list
        flex_dict[chain] = flex_list
    print(str(mut_dict))
    # Make pymol session
    try:
        cmd.load(conf_space.mol)
    except:
        raise CmdException
        print "Error! Cannot load pdb from "+conf_space.mol
    try:
        cmd.select("mut", "none")
        cmd.select("flex", "none")
    except:
        raise CmdException
        print "Error! Cannot select empty selections"
    # Select mutable residues
    for chain, resi_list in mut_dict.iteritems():
        if resi_list:
            mut_selection = "mut or (chain "+chain+" and resi "+\
                    "+".join(resi_list)+")"
            try:
                cmd.select("mut", mut_selection)
            except:
                raise CmdException
    # Select flexible residues
    for chain, resi_list in flex_dict.iteritems():
        if resi_list:
            flex_selection = "flex or (chain "+chain+" and resi "+\
                    "+".join(resi_list)+")"
            try:
                cmd.select("flex", flex_selection)
            except:
                raise CmdException
    # Show sticks and/or lines
    try:
        cmd.show("sticks", "mut")
        cmd.show("lines", "flex")
        cmd.hide("everything", "hydrogens")
    except:
        raise CmdException

# Make this executable as a command from pymol
cmd.extend( "print_design", print_design )
cmd.extend( "test_gen", test_gen )
cmd.extend( "gen_shell", gen_shell )
cmd.extend( "interface", interface )
cmd.extend( "load_cfs", load_cfs )

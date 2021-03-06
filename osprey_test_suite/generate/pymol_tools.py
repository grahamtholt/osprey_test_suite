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

SIMILAR_MUTATIONS_1 = {
    'ALA': ['GLY', 'SER', 'THR', 'VAL', 'ALA'],
    'VAL': ['VAL', 'LEU', 'ILE', 'PHE', 'TYR'],
    'LEU': ['VAL', 'LEU', 'ILE', 'PHE', 'TYR'],
    'ILE': ['VAL', 'LEU', 'ILE', 'PHE', 'TYR'],
    'PHE': ['VAL', 'LEU', 'ILE', 'PHE', 'TYR'],
    'TYR': ['VAL', 'LEU', 'ILE', 'PHE', 'TYR'],
    'TRP': ['TRP', 'LEU', 'ILE', 'PHE', 'TYR'],
    'CYS': ['ASP', 'GLU', 'ASN', 'GLN', 'HID'],
    'MET': ['VAL', 'LEU', 'ILE', 'PHE', 'TYR'],
    'SER': ['GLY', 'SER', 'THR', 'VAL', 'ALA'],
    'THR': ['GLY', 'SER', 'THR', 'VAL', 'ALA'],
    'LYS': ['LYS', 'ARG', 'GLN', 'ASN', 'HID', 'THR'],
    'ARG': ['LYS', 'ARG', 'GLN', 'ASN', 'HID', 'THR'],
    'HIP': ['ASP', 'GLU', 'ASN', 'GLN', 'HID', 'THR'],
    'HIE': ['ASP', 'GLU', 'ASN', 'GLN', 'HID', 'THR'],
    'HID': ['ASP', 'GLU', 'ASN', 'GLN', 'HID', 'THR'],
    'ASP': ['ASP', 'GLU', 'ASN', 'GLN', 'HID', 'THR'],
    'GLU': ['ASP', 'GLU', 'ASN', 'GLN', 'HID', 'THR'],
    'ASN': ['ASP', 'GLU', 'ASN', 'GLN', 'HID', 'THR'],
    'GLN': ['ASP', 'GLU', 'ASN', 'GLN', 'HID', 'THR'],
    'GLY': ['GLY', 'SER', 'THR', 'VAL', 'ALA', 'THR']
}
"""dict: Dictionary of amino acid templates by AA. Goal is to provide a smaller
set of mutations for designs. Current limit is 5 mutations.

"""

MODULE_LOC = os.path.dirname((os.path.realpath(__file__)))

DEFAULT_ROT_FILE = MODULE_LOC + "/../../resources/datafiles/LovellRotamer.dat"
"""The default rotamer file
"""

HEADER = """#The following information is sufficient to generate an OSPREY 3.0
#configuration space"""
"""The header to write out to the config file
"""

class TooManyChainsException(Exception):
    pass

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
            print("Error in selection: "+selection)
            raise CmdException

        # Remove empty selections
        stored.test = []
        try:
            cmd.iterate(sel_name+" and name ca",
                        "stored.test.append(\"\"+resi)")
        except:
            print("Cannot iterate through "+sel_name)
            raise CmdException
        if not stored.test:
            try:
                cmd.delete(sel_name)
            except:
                print("Cannot delete "+sel_name)
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
            if flex.min_sc_sc_distance(mut) < (int(r)-int(aggr)*0.2):
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

    if int(prune_sele):
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

def test_gen(max_num_mut=6, target_num_mut=None, flex_only=False,
             limit_seqspace=False):
    """Generate all reasonable designs from the current pdb

    Args:
        max_num_mut (int): The maximum number of mutations to allow
        target_num_mut (int): Only return mutable sets of this size
        flex_only (boolean): Allow mutations?

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
        print(inter_1)
        mut_list.extend(gen_mut(inter_1, max_num_mut, target_num_mut))
        print(inter_2)
        mut_list.extend(gen_mut(inter_2, max_num_mut, target_num_mut))


    print('Generated '+str(len(mut_list))+' sets of mutations')
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
        name = '_'.join([obj_name[:4], chain])+".cfs"
        print_design("mut", "flex", name, flex_only, limit_seqspace)
        counter = counter+1

    #optional, but slow
    #gen_scenes(mut_list)

def gen_mut(sele_name, max_num_mut, target_num_mut=None):
    """Select combinations of mutable residues

    Selects all possible combinations of good mutable residues within an
    interface.

    Args:
        sele_name (string): The name of a PyMol selection that defines a PPI
        max_num_mut (int): The maximum number of mutations to allow
        target_num_mut (int): Only return mutable sets of this size

    Returns:
        A set of sets of Residue objects that are good design choices
    """
    mut_list = set()
    # Added condition that we don't consider prolines
    res_list = [res for res in residue.Residue.sele_to_res(sele_name)
                if res.res_name != "PRO" and res.res_name != "CYX"]
    # Since we are working with sets, we'll just work with indices
    res_dict = dict(zip(list(range(0,len(res_list))),
                        res_list))

    # Add all single mutation indices
    mut_list.update(frozenset((i,)) for i in res_dict.keys())

    # Determine how many mutations we want
    if target_num_mut is not None:
        max_num_mut = target_num_mut

    if max_num_mut < 2:
        return mut_list

    # Construct close pairs
    close_pairs = set( frozenset(e) for e in combinations(res_dict.keys(), 2)
                        if res_dict[e[0]].min_sc_sc_distance(res_dict[e[1]]) < 3 )

    # Add close pairs to muts
    mut_list.update(close_pairs)

    # Iteratively union close pairs until cannot any longer
    # Or, if we reach a maximum of 6 mutable residues
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

    #convert back from indices to residue objects
    def recover_res(tup):
        return tuple([res_dict[e] for e in tup])

    if target_num_mut is not None:
        return [recover_res(e) for e in mut_list if len(e)==target_num_mut]
    else:
        return [recover_res(e) for e in mut_list]

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

def print_design(mut="mut", flex="flex", out="design.cfs", flex_only=False,
                 limit_seqspace = False):
    """
    Print out the OSPREY config information for a design to file

    Args:
        mut (str): PyMol selection name for the mutable design residues
        flex (str): PyMol selection name for the flexible design residues
        out (str): System path to output file
    """
    try:
        ## Find the input pdb file
        #pdb_file_name = find_pdb_file(mut, flex)

        # Collect residue lists from selections
        mut_list = residue.Residue.sele_to_res(mut)
        flex_list = residue.Residue.sele_to_res(flex)

        # Find the input pdb file
        pdb_file_name = find_pdb_file(mut, flex)

        # Make sure we don't have the mutable residues in the flexible selection
        if set(e.str_long() for e in mut_list).intersection(
            set(e.str_long() for e in flex_list)):
            raise CmdException("ERROR!: You forgot to remove the mutable residues from the"+
                " flexible ones!")

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
                allowed_mutations = [res.res_name]
                if not flex_only:
                    if limit_seqspace:
                        allowed_mutations = SIMILAR_MUTATIONS_1[res.res_name]
                    else:
                        allowed_mutations = DEFAULT_MUTATIONS_ALL

                this_strand_flex[res.chain_id+str(res.res_seq)+res.i_code] =\
                    allowed_mutations

            # For flexible residues in strand, add only res name
            chain_muts = (res for res in flex_list if res.chain_id == chain)
            for res in chain_muts:
                this_strand_flex[res.chain_id+str(res.res_seq)+res.i_code] = \
                        [res.res_name]

            strand_flex_all[strand_key] = this_strand_flex

        conf_size = get_confspace_size(strand_flex_all)
        conf_size_str = '%.3E' % conf_size
        num_res = str(get_num_res(strand_flex_all))+"res"
        # Modify out with confsize
        basename, ext = os.path.splitext(out)
        new_name = '_'.join([basename, num_res, conf_size_str])+ext
        print("Writing "+new_name)
        # Print out variables
        with open(new_name, 'w') as f:
            f.write(HEADER+"\n")
            f.write("mol = \""+str(pdb_file_name)+"\"\n")
            f.write("strand_defs = "+str(strand_defs)+"\n")
            f.write("strand_flex = "+str(strand_flex_all)+"\n")
            #for strand in strand_defs:
                #f.write(strand+" = "+str(strand_defs[strand])+"\n")
            #for strand in strand_flex_all:
                #f.write(strand+"_flex = "+str(strand_flex_all[strand])+"\n")
    except Exception as e:
        print(e)


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

    assert len(mut_obj_list) == 1, "ERROR!: mut_obj_list has length %d, not 1." %len(mut_obj_list)
    assert len(flex_obj_list) == 1, "ERROR!: flex_obj_list has length %d, not 1." %len(flex_obj_list)

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
        print("WARNING! Could not find pdb file. Please set manually")
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
        for res, mut_allowed in d.items():
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
        print("Error! Cannot load pdb from "+conf_space.mol)
    try:
        cmd.select("mut", "none")
        cmd.select("flex", "none")
    except:
        raise CmdException
        print("Error! Cannot select empty selections")
    # Select mutable residues
    for chain, resi_list in mut_dict.items():
        if resi_list:
            mut_selection = "mut or (chain "+chain+" and resi "+\
                    "+".join(resi_list)+")"
            try:
                cmd.select("mut", mut_selection)
            except:
                raise CmdException
    # Select flexible residues
    for chain, resi_list in flex_dict.items():
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

def get_confspace_size(all_flex_dict):
    """
    """
    rots = read_rotamer_info(DEFAULT_ROT_FILE)

    confspace_size = 1

    for strand_dict in all_flex_dict.values():
        strand_cfs_size = 1 # for each strand
        for res_allowed in strand_dict.values():
            pos_rots = 1 # for wild type rotamer
            for res_type in res_allowed:
                pos_rots = pos_rots + rots[res_type]
            strand_cfs_size = strand_cfs_size * pos_rots
        confspace_size = confspace_size * strand_cfs_size
    return confspace_size

def get_num_res(all_flex_dict):
    num_res = 0
    for strand_dict in all_flex_dict.values():
        num_res = num_res + len(strand_dict.keys())
    return num_res

def read_rotamer_info(rotamer_file):
    """Read rotamer dat file and return a dict
    """
    rots = dict()
    with open(rotamer_file, 'r') as f:
        data = f.read()
        m = re.findall(r"([A-Za-z]{3}) \d+ (\d+)", data)
        for match in m:
            rots[match[0]] = int(match[1])
            if match[0] == "ALA" or match[0] == "GLY":
                rots[match[0]] = rots[match[0]]+1

    return rots

# Make this executable as a command from pymol
cmd.extend( "print_design", print_design )
cmd.extend( "test_gen", test_gen )
cmd.extend( "gen_shell", gen_shell )
cmd.extend( "interface", interface )
cmd.extend( "load_cfs", load_cfs )

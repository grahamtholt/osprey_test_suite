#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Defines a class and tools to handle ATOM records in a PDB file.

Atom class contains all fields in ATOM record.
read_pdb_line method parses a single string into fields and returns an atom
object containing those fields.
"""

import math
import re
from pymol import cmd, CmdException, stored
##################################################
# Version
#   0.0.0 -- First version
#   0.1.0 -- Added methods for distances
# Author
#   Graham Holt
##################################################

__author__ = "Graham Holt"
__version__ = "0.1.0"
__status__ = "Prototype"

class Atom(object):
    """An atom based on the PDB file format. Atoms have the following
    properties:

    Attributes:
        prefix(str): ATOM or HETATM, depending on record type.
        serial(int): Atom serial number. Restricted to 5 digits.
        name(str): Atom name. Restricted to 4 characters
        alt_loc(str): Alternate location. Restricted to 1 character.
        res_name(str): Residue name. Restricted to 3 characters.
        chain_id(char): Chain identifier.
        res_seq(int): Residue number: Restricted to 4 characters.
        i_code(char): Insertion code. Distinguishes residues of same number
        loc(3-tuple): Cartesian coordinates of Atom in 3-space. Each coord is 7
            digit float (signed)
        occ(float): Occupancy of atom. Restricted to 5 digit float (signed)
        temp(float): Temperature factor. Restricted to 5 digit float (signed)
        ele(str): Element symbol, right justified. Restricted to 2 characters.
        charge(str): Atom charge. Restricted to 2 characters.

        """

    def __init__(self, prefix, serial, name, alt_loc, res_name, chain_id, res_seq, i_code, loc,
                occ, temp, ele, charge):
        self.prefix= prefix
        self.serial = serial
        self.name = name
        self.alt_loc = alt_loc
        self.res_name = res_name
        self.chain_id = chain_id
        self.res_seq = res_seq
        self.i_code = i_code
        self.loc = loc
        self.occ = occ
        self.temp = temp
        self.ele = ele
        self.charge = charge

    def __str__(self):
        return '{:6s}{:5d} {:4s}{:1s}{:3s} {:1s}{:4d}{:1s}'\
                '   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}'\
                .format(self.prefix, self.serial, self.name, self.alt_loc,\
                self.res_name, self.chain_id, self.res_seq, self.i_code,\
                self.loc[0], self.loc[1], self.loc[2], self.occ, self.temp,\
                self.ele, self.charge)

    # prefix
    @property
    def prefix(self):
        return self._prefix

    @prefix.setter
    def prefix(self, value):
        if type(value) is not str:
            raise TypeError('prefix must be a string')
        elif len(value) <= 6:
            self._prefix = value
        else:
            raise ValueError('prefix is out of bounds')

    # serial
    @property
    def serial(self):
        return self._serial

    @serial.setter
    def serial(self, value):
        if type(value) is not int:
            raise TypeError('serial must be an integer')
        elif len(str(value)) <= 5 and value >= 0:
            self._serial = value
        else:
            raise ValueError('serial is out of bounds')

    # name
    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        if type(value) is not str:
            raise TypeError('name must be a string')
        elif len(value) <= 4:
            self._name = value
        else:
            raise ValueError('name is out of bounds')

    # alt_loc
    @property
    def alt_loc(self):
        return self._alt_loc

    @alt_loc.setter
    def alt_loc(self, value):
        if type(value) is not str:
            raise TypeError('alt_loc must be a string')
        elif len(value) <= 1:
            self._alt_loc = value
        else:
            raise ValueError('alt_loc is out of bounds')

    # res_name
    @property
    def res_name(self):
        return self._res_name

    @res_name.setter
    def res_name(self, value):
        if type(value) is not str:
            raise TypeError('res_name must be a string')
        elif len(value) <= 3:
            self._res_name = value
        else:
            raise ValueError('res_name is out of bounds')

    # chain_id
    @property
    def chain_id(self):
        return self._chain_id

    @chain_id.setter
    def chain_id(self, value):
        if type(value) is not str:
            raise TypeError('chain_id must be a string')
        elif len(value) <= 1:
            self._chain_id = value
        else:
            raise ValueError('chain_id is out of bounds')

    # res_seq
    @property
    def res_seq(self):
        return self._res_seq

    @res_seq.setter
    def res_seq(self, value):
        if type(value) is not int:
            raise TypeError('res_seq must be an int')
        elif len(str(value)) <= 4:
            self._res_seq = value
        else:
            raise ValueError('res_seq is out of bounds')

    # i_code
    @property
    def i_code(self):
        return self._i_code

    @i_code.setter
    def i_code(self, value):
        if type(value) is not str:
            raise TypeError('i_code must be a string')
        elif len(value) <= 1:
            self._i_code = value
        else:
            raise ValueError('i_code is out of bounds')

    # loc
    @property
    def loc(self):
        return self._loc

    @loc.setter
    def loc(self, value):
        if type(value) is not tuple:
            raise TypeError('loc must be a tuple')
        elif len(value) is 3:
            if (len(str(value[0])) <= 8 and len(str(value[1])) <=8 and
                                               len(str(value[2])) <=8):
                self._loc = value
            else:
                #raise ValueError('One or more coords in loc is out of bounds')
                # Don't throw error, just round when printing
                self._loc = value
        else:
            raise ValueError('loc must have 3 coordinates')

    # occ
    @property
    def occ(self):
        return self._occ

    @occ.setter
    def occ(self, value):
        if type(value) is not float:
            raise TypeError('occ must be a float')
        elif len(str(value)) <= 6:
            self._occ = value
        else:
            raise ValueError('occ is out of bounds')

    # temp
    @property
    def temp(self):
        return self._temp

    @temp.setter
    def temp(self, value):
        if type(value) is not float:
            raise TypeError('temp must be a float')
        else:
            self._temp = value

    # ele
    @property
    def ele(self):
        return self._ele

    @ele.setter
    def ele(self, value):
        if type(value) is not str:
            raise TypeError('ele must be a string')
        elif len(value) <= 2:
            self._ele = value
        else:
            raise ValueError('ele is out of bounds')

    # charge
    @property
    def charge(self):
        return self._charge

    @charge.setter
    def charge(self, value):
        if type(value) is not str:
            raise TypeError('charge must be a string')
        elif len(value) <= 2:
            self._charge = value
        else:
            raise ValueError('charge is out of bounds')

    def is_backbone(self):
        bb_atoms = ['C', 'CA', 'O', 'N', 'H', 'HA' ]
        if self.name in bb_atoms:
            return True
        else:
            return False

    def distance(self, other):
        """Calculates the cartesian distance to other atom

        @param other: The atom object to which to calculate distance
        """
        return math.sqrt((self.loc[0] - other.loc[0])**2 +
                         (self.loc[1] - other.loc[1])**2 +
                         (self.loc[2] - other.loc[2])**2
                        )

    @staticmethod
    def str_to_atom(atom_line):
        """Parse a string into an Atom object.

        Args:
            atom_line(str): A string from a pdb file. Must begin with ATOM

        Returns:
            Atom object

        """
        # Pad string to ensure sufficient length
        atom_line = atom_line.ljust(80)

        # Ensure line begins with 'ATOM' or 'HETATM'
        if atom_line[:6].strip() != 'ATOM' and atom_line[:6].strip() != 'HETATM':
            raise ValueError('Error: Input string is neither an ATOM nor a HETATM.')

        # The splits here are defined in wwpdb.org/documentation
        return Atom(atom_line[:6].strip(),
                    int(atom_line[6:11].strip()), atom_line[12:16].strip(),
                    atom_line[16].strip(), atom_line[17:20].strip(),
                    atom_line[21].strip(), int(atom_line[22:26].strip()),
                    atom_line[26].strip(), (float(atom_line[30:38].strip()),
                                            float(atom_line[38:46].strip()),
                                            float(atom_line[46:54].strip())),
                    float(atom_line[54:60].strip()),
                    float(atom_line[60:66].strip()), atom_line[76:78].strip(),
                    atom_line[78:80].strip())

    @staticmethod
    def str_to_atoms(atom_string):
        """Parse a string into a list of Atom objects.

        Args:
            atom_line(str): A pdb formatted string

        Returns:
            A list of Atom objects

        """
        atom_list = [Atom.str_to_atom(line) for line in atom_string.splitlines()]
        return atom_list

    @staticmethod
    def sele_to_atoms(sele_name):
        """Generate a list of atom objects from a PyMol selection

        Args:
            sele_name (string): PyMol selection name

        Returns:
            A list of Atom objects
        """
        atom_list = []
        # Get the list of all atoms in the selection
        stored.atoms = []
        try:
            cmd.iterate(sele_name,
                        "stored.atoms.append([type, ID, name, alt, resn, chain,\
                        resi,\"\",(0,0,0), q, b, elem, \"\"])")
            # Get coordinates
            # NOTE: PYMOL MIGHT HAVE ROUNDING ERRORS FOR COORDINATES!
            # The below only works in PyMol 1.7.4 and beyond
            #coords_list = cmd.get_coords(sele_name, 1)
            coords_list = cmd.get_model(sele_name, 1).get_coord_list()
        except:
            raise CmdException

        # Write coordinates, extract icode, make atom lists, reslist
        for e, coord in zip(stored.atoms, coords_list):
            e[8] = tuple(coord)
            m = re.match(r"(-?\d+)(\w?)", e[6])
            try:
                e[6] = int(m.group(1))
                e[7] = m.group(2)
            except:
                print("Error in resi" + e[6])
                raise Exception

            #make residues
            a = Atom(*e)
            atom_list.append(a)

        return atom_list

    @staticmethod
    def pdb_to_atoms(file_path):
        """Parse a pdb file (wwpdb.org) and return a list of atoms contained in
            the pdb file.

            Args:
                file_path(str): The location of the pdb file to read.

            Returns:
                A list of Atom objects

        """
        atom_list = []
        with open(file_path) as fp:
            for line in fp:
                if line[:6].strip() == 'ATOM' or line[:6].strip() == 'HETATM':
                    a = Atom.str_to_atom(line)
                    atom_list.extend(a)
        return atom_list

#EOF

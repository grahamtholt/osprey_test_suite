#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Defines a class and tools to define a pymol-based amino-acid residue

The Residue class is a collection of Atom objects.

ToDo:
    * Implement TER handling
    * Add method to calc dist between side chain and any residue in other
"""

##################################################
# Version
#   0.0.0  -- First version
#   0.1.0 -- Added methods for distances
#   0.2.0 -- Added methods for ca-cb vectors
# Author
#   Graham Holt
##################################################
import re
from numpy import dot
from numpy.linalg import norm
from numpy import subtract
from math import acos, cos
from math import pi
from pymol import cmd, CmdException, stored
import atom

__author__ = "Graham Holt"
__version__ = "0.2.0"
__status__ = "Prototype"

class Residue(object):
    """A residue containing an arbitrary number of Atom objects. Residues have
    the following properties:

    Any change to residue attributes (other than atom_list) will change the
        attributes of any Atoms in atom_list.

    Attributes:
        atom_list(list<Atom>): A list of Atom objects.
        res_name(str): The name of the residue. Restricted to three characters.
        res_seq(int): Residue sequence number. Follows Atom class restrictions.
        chain_id(char): Chain identifier. Follows Atom class restrictions.
        i_code(char): Code for insertion of residues. Follows Atom class
            restrictions.

    """

    def __init__(self, atom_list = []):
        """Instantiate a Residue object with a given list of atoms. Defaults to
        an empty list.

        Args:
            atom_list(list): A list of Atom objects. Defaults to empty list.

        """

        self.atom_list = atom_list
        self.res_name = atom_list[0].res_name
        self.res_seq = atom_list[0].res_seq
        self.chain_id = atom_list[0].chain_id
        self.i_code = atom_list[0].i_code

    # atom_list
    @property
    def atom_list(self):
        return self._atom_list

    @atom_list.setter
    def atom_list(self,value):
        if not (all(e.res_seq == value[0].res_seq for e in value) and
                all(e.res_name == value[0].res_name for e in value) and
                all(e.chain_id == value[0].chain_id for e in value) and
                all(e.i_code == value[0].i_code for e in value)):
            raise ValueError('All Atoms in Residue must have same values for'\
                            'res_seq, res_name, chain_id, and i_code.')
        else:
            self._atom_list = value

    @property
    def res_name(self):
        return self._res_name

    @res_name.setter
    def res_name(self, value):
        for e in self.atom_list:
            e.res_name = value
        self._res_name = value

    @property
    def res_seq(self):
        return self._res_seq

    @res_seq.setter
    def res_seq(self, value):
        for e in self.atom_list:
            e.res_seq = value
        self._res_seq = value

    @property
    def chain_id(self):
        return self._chain_id

    @chain_id.setter
    def chain_id(self, value):
        for e in self.atom_list:
            e.chain_id = value
        self._chain_id = value

    @property
    def i_code(self):
        return self._i_code

    @i_code.setter
    def i_code(self, value):
        for e in self.atom_list:
            e.i_code = value
        self._i_code = value


    def add_atom(self, atom):
        """Add an atom to the atom_list"""
        if not self.can_contain(atom):
            raise ValueError('All Atoms in Residue must have same values for'\
                            'res_seq, res_name, chain_id, and i_code.')
        else:
            self.atom_list.append(atom)

    def can_contain(self, atom):
        """Returns True if an Atom has res_name, res_seq, chain_id, and i_code
        values that match with the Residue object.

        Args:
            atom(Atom): The Atom that we want to check against the Residue

        Returns:
            boolean value

        """
        if (atom.res_name == self.res_name and atom.res_seq == self.res_seq
                and atom.chain_id == self.chain_id
                and atom.i_code == self.i_code):
            return True
        else:
            return False

    def remove_atom(self, atom):
        """Remove an atom from the atom_list"""
        self.atom_list.remove(atom)

    def remove_atom_byname(self, atom_name):
        """Remove an atom from the atom_list by specifying atom.name.

        Args:
            atom_name(str): The name field of the atom to be removed.

        """
        for e in self.atom_list:
            if e.name == atom_name:
                self.atom_list.remove(e)

    def res_to_pdb(self):
        return ("\n".join((str(e) for e in self.atom_list)))+"\n"

    def __str__(self):
        return ''.join((self.res_name, str(self.res_seq)))

    def __eq__(self, other):
        """Overrides the default implementation"""
        if isinstance(other, Residue):
            return ( self.chain_id+str(self.res_seq)+self.i_code ==
                    other.chain_id+str(other.res_seq)+other.i_code )
        return False

    def __ne__(self, other):
        """Overrides the default implementation (Unnecessary in Python 3)"""
        return not self.__eq__(other)


    def min_distance(self, other):
        """Calculate the smallest distance between two residues

        Calculates the smallest distance between any atom in this residue and
        any atom in the other residue
        """
        distance = float("inf")
        for atom in other.atom_list:
            for myatom in self.atom_list:
                new_dist = myatom.distance(atom)
                if new_dist < distance:
                    distance = new_dist

        return distance

    def max_distance(self, other):
        """Calculate the largest distance between two residues

        Calculates the largest distance between any atom in this residue and
        any atom in the other residue
        """
        distance = float("-inf")
        for atom in other.atom_list:
            for myatom in self.atom_list:
                new_dist = myatom.distance(atom)
                if new_dist > distance:
                    distance = new_dist

        return distance

    def min_sc_any_distance(self, other):
        """Calculate the smallest distance between side chain and residue

        Calculates the smallest distance between any sidechain atom in
        this residue and any atom in the other residue
        """
        distance = float("inf")
        for atom in other.atom_list:
            for myatom in self.atom_list:
                if not myatom.is_backbone():
                    new_dist = myatom.distance(atom)
                    if new_dist < distance:
                        distance = new_dist

        return distance

    def max_sc_any_distance(self, other):
        """Calculate the largest distance between side chain and residue

        Calculates the largest distance between any sidechain atom in
        this residue and any atom in the other residue
        """
        distance = float("-inf")
        for atom in other.atom_list:
            for myatom in self.atom_list:
                if not myatom.is_backbone():
                    new_dist = myatom.distance(atom)
                    if new_dist > distance:
                        distance = new_dist

        return distance

    def min_sc_sc_distance(self, other):
        """Calculate the smallest distance between side chains

        Calculates the smallest distance between any sidechain atom in
        this residue and any sidechain atom in the other residue
        """
        distance = float("inf")
        for atom in other.atom_list:
            if not atom.is_backbone():
                for myatom in self.atom_list:
                    if not myatom.is_backbone():
                        new_dist = myatom.distance(atom)
                        if new_dist < distance:
                            distance = new_dist

        return distance

    def max_sc_sc_distance(self, other):
        """Calculate the largest distance between side chains

        Calculates the largest distance between any sidechain atom in
        this residue and any sidechain atom in the other residue
        """
        distance = float("-inf")
        for atom in other.atom_list:
            if not atom.is_backbone():
                for myatom in self.atom_list:
                    if not myatom.is_backbone():
                        new_dist = myatom.distance(atom)
                        if new_dist > distance:
                            distance = new_dist

        return distance
    def get_calpha(self):
        """Return c-alpha atom"""
        calpha = [ a for a in self.atom_list if a.name == 'CA' ]
        if len(calpha) != 1:
            print calpha
        assert len(calpha)==1
        return calpha[0]

    def get_cbeta(self):
        """Return c-beta atom"""
        cbeta = [ a for a in self.atom_list if a.name == 'CB' ]

        assert len(cbeta)==1
        return cbeta[0]

    def angles_diverge(self, other, angle=130):
        """
        """
        # Since glycines don't have c-betas, just prune, since we want to prune
        # them anyway
        if self.res_name == "GLY":
            return True
        # Convert to radians
        #print("Original angle is %f" % angle)
        angle = angle * pi / 180
        #print("Radian angle is %f" % angle)
        # Create vectors:
        # Vector goes from ca of other (mutable) to ca of this (flexible)
        ca_ca = subtract(self.get_calpha().loc, other.get_calpha().loc)
        # Vector goes from ca of other (mutable) to cb of this (flexible)
        ca_cb = subtract(self.get_cbeta().loc, other.get_calpha().loc)
        # Vector goes from ca to cb of other (mutable)
        mut_ca_cb = subtract(other.get_cbeta().loc, other.get_calpha().loc)

        # Normalize all vectors
        ca_ca_unit = ca_ca/norm(ca_ca)
        ca_cb_unit = ca_cb/norm(ca_cb)
        mut_ca_cb_unit = mut_ca_cb/norm(mut_ca_cb)

        side_chain = subtract(self.get_calpha().loc,
                                    self.get_cbeta().loc)

        # Generate dot products
        ca_dot = dot(ca_ca_unit, mut_ca_cb_unit)
        cb_dot = dot(ca_cb_unit, mut_ca_cb_unit)

        # Find angle between vectors
        ca_angle =  acos(ca_dot)
        cb_angle =  acos(cb_dot)

        # For Debugging
        #ca_angle_deg = ca_angle * 180 / pi
        #cb_angle_deg = cb_angle * 180 / pi
        #pruning_angle_deg = angle *180 / pi
        #print(self.res_name+str(self.res_seq), other.res_name+str(other.res_seq))
        #print("Pruning angle: %f" % pruning_angle_deg)
        #print("Angle: (cb_mut, ca_mut, ca_flex): %f" % ca_angle_deg)
        #print("Angle: (cb_mut, ca_mut, cb_flex): %f" % cb_angle_deg)

        if ca_angle > angle and cb_angle > angle:
            return True
        else:
            return False
    @staticmethod
    def sele_to_res(sele_name):
        """Generate a list of residue objects from a PyMol selection

        Args:
            sele_name (string): PyMol selection name

        Returns:
            A list of Residue objects
        """

        res_list = []
        for a in atom.Atom.sele_to_atoms(sele_name):
            if not res_list or not res_list[-1].can_contain(a):
                res_list.append(Residue([a]))
            else:
                res_list[-1].add_atom(a)

        return res_list

    @staticmethod
    def pdb_to_res(file_path):
        """Parse a pdb file (wwpdb.org) and return a list of residues contained in
            the pdb file.

            Args:
                file_path(str): The location of the pdb file to read.

            Returns:
                A list of Residue objects

        """
        res_list = []
        for a in atom.Atom.pdb_to_atoms(file_path):
            if not res_list or not res_list[-1].can_contain(a):
                res_list.append(Residue([a]))
            else:
                res_list[-1].add_atom(a)
        return res_list

    @staticmethod
    def str_to_res(pdb_string):
        """Parse a pdb formatted string (wwpdb.org) and return a list of
        residues contained in the pdb file.

            Args:
                pdb_string (str): The pdb formatted text

            Returns:
                A list of Residue objects

        """
        res_list = []
        for a in atom.Atom.str_to_atoms(pdb_string):
            if not res_list or not res_list[-1].can_contain(a):
                res_list.append(Residue([a]))
            else:
                res_list[-1].add_atom(a)
        return res_list

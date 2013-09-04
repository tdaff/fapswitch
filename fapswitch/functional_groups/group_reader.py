#!/usr/bin/env python

"""
fapswitch.py

Alter a known structure with new functional groups ready for fapping.

"""


import copy
try:
    import configparser
except ImportError:
    import ConfigParser as configparser
import hashlib
import glob
try:
    import cPickle as pickle
except ImportError:
    import pickle
import random
import re
import socket
import sys
import textwrap
import time
from itertools import chain, combinations, product
from logging import debug, info, warn, error, critical
from os import path

import numpy as np
from numpy import array, identity, asarray, dot, cross, outer, sin, cos
from numpy import roll
from numpy.linalg import norm

from fapswitch.functional_groups import functiaonl_groups
from fapswitch.core.components import Structure, Atom
from fapswitch.core.components import vecdist3, subgroup
from fapswitch.config.config import Options
from fapswitch.core.elements import CCDC_BOND_ORDERS


DOT_FAPSWITCH_VERSION = (6, 0)


class FunctionalGroupLibrary(dict):
    """
    Container for all the available functional groups just subclasses
    the standard dict adding some new methods.

    """
    def __init__(self, *args, **kwargs):
        """Initialise a standard dict and update with functional groups."""
        # Make sure that the dict is properly initilaised
        dict.__init__(self, *args, **kwargs)

        # Functional groups can be put in several places
        # Defaults '*.lib' files in library directory always loaded
        library_glob = path.join(path.dirname(path.realpath(__file__)),
                                 'library', '*.flib')
        for library in glob.glob(library_glob):
            self._from_file(path.join(library))

        # Custom user groups stored in .faps supercede defaults
        faps_glob = path.join(path.expanduser('~'), '.faps', '*.flib')
        for library in glob.glob(faps_glob):
            self._from_file(path.join(library))

        # Job specific groups in the working directory
        for library in glob.glob('*.flib'):
            self._from_file(path.join(library))

    def _from_file(self, library_file_name='functional_groups.flib'):
        """Parse groups from the configparser .ini style file."""
        # just a standard configparser conversion to a dict of
        # FunctionalGroup objects
        library_file = configparser.SafeConfigParser()
        debug("Reading groups from %s" % library_file_name)
        library_file.read(library_file_name)
        for group_name in library_file.sections():
            try:
                if group_name in self:
                    debug("Overriding group %s" % group_name)
                self[group_name] = FunctionalGroup(library_file.items(group_name))
            except KeyError:
                error("Group %s is missing data; update library" % group_name)

    @property
    def group_list(self):
        """Compile a soretd list of all the available groups."""
        return sorted(list(self))


class FunctionalGroup(object):
    """
    Substitutable functional group. Bunch of atoms with information on how
    to connect to the framework

    """

    def __init__(self, items):
        """Initialize from a list of tuples as the attributes."""
        # These are defaults, best that they are overwritten
        self.atoms = []
        self.bonds = {}
        self.orientation = [0, 1, 0]

        # pop the items from a dict giving neater code
        items = dict(items)

        self._parse_atoms(items.pop('atoms'))
        self.orientation = normalise(string_to_tuple(items.pop('orientation'), float))
        self.normal = normalise(string_to_tuple(items.pop('normal'), float))
        self.bond_length = float(items.pop('carbon_bond'))
        self._parse_bonds(items.pop('bonds'))
        self.idx = 0
        self.connection_point = 0  # always connect to the first atom
        # Arbitrary attributes can be set
        self.__dict__.update(items)
        self._gen_neighbours()

    def _parse_atoms(self, atom_text):
        """Read atom information from the file."""
        self.atoms = []
        for atom in atom_text.splitlines():
            atom = atom.strip().split()
            if not atom:
                continue
            new_atom = Atom(atom[0], [float(x) for x in atom[2:5]])
            new_atom.uff_type = atom[1]
            new_atom.site = label_atom(new_atom.element)
            self.atoms.append(new_atom)

    def _parse_bonds(self, bond_block):
        """
        Extract the bonds from the text and calculate the distances between the
        atoms as these are needed for the cif file.

        """

        bond_split = subgroup(bond_block.split(), width=3)
        for bond_trio in bond_split:
            bond = (int(bond_trio[0]), int(bond_trio[1]))
            distance = vecdist3(self.atoms[bond[0]].pos,
                                self.atoms[bond[1]].pos)
            self.bonds[bond] = (distance, float(bond_trio[2]))

    def _gen_neighbours(self):
        """Update atoms with neighbouring atoms."""
        # Iterate over all pairs
        # Assume non periodic and geometries are good
        # dummy atom at the tether point
        self.atoms.append(Atom('C', [-self.bond_length*x for x in self.orientation]))
        for atom in self.atoms:
            # distance matrix for all neighbours except self
            neighbours = []
            for ot_idx, other in enumerate(self.atoms):
                if atom is other:
                    continue
                length = vecdist3(atom.pos, other.pos)
                neighbours.append((length, ot_idx))
            # these are expected to be is distance order
            atom.neighbours = sorted(neighbours)
        # remove the dummy atom
        self.atoms.pop()
        # first atom should be bonded to the tether, flag this as '-1'
        #self.atoms[0].bonds[self.atoms[0].bonds.index(self.natoms)] = -1

    def atoms_attached_to(self, point, direction, normal, attach_point, start_index, bond_length=None):
        """Return a list of atoms at the specified position."""
        if bond_length is None:
            bond_length = self.bond_length
        new_atoms = [copy.copy(atom) for atom in self.atoms]
        rotate_matrix = matrix_rotate(self.orientation, direction)
        my_rotated_normal = np.dot(rotate_matrix, self.normal)
        orient_matrix = matrix_rotate(my_rotated_normal, normal)
        for index, atom in enumerate(new_atoms):
            atom.pos = np.dot(orient_matrix, np.dot(rotate_matrix, atom.pos))
            atom.pos = (atom.pos + point + bond_length*np.array(direction))
            atom.idx = start_index + index
        new_bonds = {}
        for bond_pair, bond_info in self.bonds.items():
            new_bond = (bond_pair[0] + start_index, bond_pair[1] + start_index)
            new_bonds[new_bond] = bond_info
        # bond to structure is single...
        bond_to_structure = (attach_point, self.connection_point + start_index)
        new_bonds[bond_to_structure] = (self.bond_length, 1.0)

        return new_atoms, new_bonds

    @property
    def natoms(self):
        """The number of atoms in the functional group."""
        return len(self.atoms)


def matrix_rotate(source, target):
    """Create a rotation matrix that will rotate source on to target."""
    # Normalise so there is no scaling in the array
    source = asarray(source)/norm(source)
    target = asarray(target)/norm(target)
    v = cross(source, target)
    vlen = dot(v, v)
    if vlen == 0.0:
        # already aligned, no rotation needed
        return identity(3)
    c = dot(source, target)
    h = (1 - c)/vlen
    return array([[c + h*v[0]*v[0], h*v[0]*v[1] - v[2], h*v[0]*v[2] + v[1]],
                  [h*v[0]*v[1] + v[2], c + h*v[1]*v[1], h*v[1]*v[2] - v[0]],
                  [h*v[0]*v[2] - v[1], h*v[1]*v[2] + v[0], c + h*v[2]*v[2]]])


def normalise(vector):
    """Return an array with magnitude 1."""
    return asarray(vector)/norm(vector)

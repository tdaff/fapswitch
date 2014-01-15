"""
group_reader.py

Containers to parse and hold for the functional group information.

"""


import copy
try:
    from configparser import ConfigParser
except ImportError:
    from ConfigParser import SafeConfigParser as ConfigParser
import glob
import re
from collections import OrderedDict
from os import path

import numpy as np
from numpy import array, identity, asarray, dot, cross
from numpy.linalg import norm

from fapswitch.config import options
from fapswitch.config import info, debug, error
from fapswitch.core.components import Atom
from fapswitch.core.elements import UFF
from fapswitch.core.util import vecdist3, subgroup


class FunctionalGroupLibrary(OrderedDict):
    """
    Container for all the available functional groups just subclasses
    the standard dict adding some new methods.

    """
    def __init__(self, *args, **kwargs):
        """Initialise a dictionary and update with functional groups."""
        # Make sure that the dict is properly initilaised
        OrderedDict.__init__(self, *args, **kwargs)

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

    def _from_file(self, flib_file_name='functional_groups.flib'):
        """Parse groups from the configparser .ini style file."""
        # just a standard configparser conversion to a dict of
        # FunctionalGroup objects
        mepo_only = options.getbool('mepo_only')
        flib_file = ConfigParser()
        debug("Reading groups from {}".format(flib_file_name))
        flib_file.read(flib_file_name)
        for group_name in flib_file.sections():
            try:
                new_group = FunctionalGroup(group_name,
                                            flib_file.items(group_name))
                if mepo_only and not new_group.mepo_compatible:
                    debug("Skipped non-MEPO {}".format(group_name))
                else:
                    if group_name in self:
                        debug("Overriding group {}".format(group_name))
                    self[group_name] = new_group
            except KeyError:
                error("Group {} is missing data".format(group_name))

    @property
    def group_list(self):
        """Compile a soretd list of all the available groups."""
        return sorted(list(self))


class FunctionalGroup(object):
    """
    Substitutable functional group. Bunch of atoms with information on how
    to connect to the framework

    """

    def __init__(self, group_name, items):
        """Initialize from a list of tuples as the attributes."""
        # These are defaults, best that they are overwritten
        self.ident = group_name
        self.atoms = []
        self.bonds = {}
        self.orientation = [0, 1, 0]
        self.smiles = ""
        self.mepo_compatible = False

        # pop the items from a dict giving neater code
        items = dict(items)

        self._parse_atoms(items.pop('atoms'))
        self.orientation = normalise(string_to_tuple(items.pop('orientation'),
                                                     float))
        self.normal = normalise(string_to_tuple(items.pop('normal'), float))
        self.bond_length = float(items.pop('carbon_bond'))
        self._parse_bonds(items.pop('bonds'))
        self.connection_point = 0  # always connect to the first atom
        if 'mepo_compatible' in items:
            # This parsing is taken from config.py as it doesn't like "False"
            compat = items.pop('mepo_compatible').lower()
            if compat in ["1", "yes", "true", "on"]:
                self.mepo_compatible = True
        # FIXME(tdaff) Add the synthesis stuff
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
            new_atom.vdw_radius = UFF[new_atom.type][0]/2.0
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
        self.atoms.append(Atom('C', [-self.bond_length*x
                                     for x in self.orientation]))
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

    def atoms_attached_to(self, point, direction, normal, attach_point,
                          start_index):
        """Return a list of atoms at the specified position."""
        new_atoms = [copy.copy(atom) for atom in self.atoms]
        rotate_matrix = matrix_rotate(self.orientation, direction)
        my_rotated_normal = np.dot(rotate_matrix, self.normal)
        orient_matrix = matrix_rotate(my_rotated_normal, normal)
        for index, atom in enumerate(new_atoms):
            atom.pos = np.dot(orient_matrix, np.dot(rotate_matrix, atom.pos))
            atom.pos = (atom.pos + point + self.bond_length*np.array(direction))
            atom.idx = start_index + index
        new_bonds = {}
        for bond_pair, bond_info in self.bonds.items():
            new_bond = (bond_pair[0] + start_index, bond_pair[1] + start_index)
            new_bonds[new_bond] = bond_info
        # bond to structure is single...
        bond_to_structure = (attach_point, self.connection_point + start_index)
        new_bonds[bond_to_structure] = (self.bond_length, 1.0)

        return new_atoms, new_bonds

    def log_info(self):
        """Send all the information about the group to the logging functions."""
        info("[{}]".format(self.ident))
        info("name = {}".format(self.name))
        info("smiles = {}".format(self.smiles))
        info("mepo_compatible = {}".format(self.mepo_compatible))

        debug("atoms =")
        for atom in self.atoms:
            debug("    {0:4} {1:5} {2[0]:10.6f} {2[1]:10.6f} {2[2]:10.6f}".format(atom.type, atom.uff_type, atom.pos))
        debug("orientation = {0[0]:.1f} {0[1]:.1f} {0[2]:.1f}".format(self.orientation))
        debug("normal = {0[0]:.1f} {0[1]:.1f} {0[2]:.1f}".format(self.normal))
        debug("carbon_bond = {}".format(self.bond_length))
        debug("bonds =")
        for bond in self.bonds:
            debug("    {0[0]:4} {0[1]:4} {1[1]:5.2f}".format(bond, self.bonds[bond]))
        info("")

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

def string_to_tuple(value, dtype=None):
    """Parse a list of items, ignoring whitespace, brackets and commas."""
    value = [x for x in re.split(r'[\s,\(\)\[\]]*', value) if x]
    if dtype is not None:
        return tuple([dtype(x) for x in value])
    else:
        return tuple(value)

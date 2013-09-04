#!/usr/bin/env python

"""
fapswitch core components

Includes structures and atoms

"""

import os
import re
import shlex
import shutil
from itertools import count
from logging import warning, debug, error, info
from math import ceil
from os import path

import numpy as np
from numpy import pi, cos, sin, sqrt, arccos
from numpy import array, identity, dot, cross
from numpy.linalg import norm

from fapswitch.core.elements import WEIGHT, ATOMIC_NUMBER, UFF
from fapswitch.core.elements import CCDC_BOND_ORDERS, METALS
from fapswitch.core.elements import COVALENT_RADII

# Global constants
DEG2RAD = pi / 180.0

class Structure(object):
    """
    Hold information on atom positions and cell parameters.
    
    Slim version for use with cif and fapswitch only.

    """
    def __init__(self, name):
        """Just instance an empty structure initially."""
        self.name = name
        self.cell = Cell()
        self.atoms = []
        self.space_group = None

    def from_cif(self, filename=None, string=None):
        """Genereate structure from a .cif file."""
        if filename is not None:
            info("Reading positions from cif file: %s" % filename)
            filetemp = open(filename)
            cif_file = filetemp.readlines()
            filetemp.close()
        elif string is not None:
            info("Positions from cif string")
            cif_file = string.splitlines()
        else:
            error("No source for cif file")
        cif_file = strip_blanks(cif_file)
        params = [None, None, None, None, None, None]
        atoms = []
        cif_bonds = {}
        symmetry = []
        loops = []
        idx = 0
        while idx < len(cif_file):
            # line text needs to be manageable; cif guidelines can be
            # permissive
            # Can no longer just check for _underscores in lines
            # as UFF types can have them and mess up parsing
            line = cif_file[idx].lower().strip()
            if '_cell_length_a' in line:
                params[0] = ufloat(line.split()[1])
            elif '_cell_length_b' in line:
                params[1] = ufloat(line.split()[1])
            elif '_cell_length_c' in line:
                params[2] = ufloat(line.split()[1])
            elif '_cell_angle_alpha' in line:
                params[3] = ufloat(line.split()[1])
            elif '_cell_angle_beta' in line:
                params[4] = ufloat(line.split()[1])
            elif '_cell_angle_gamma' in line:
                params[5] = ufloat(line.split()[1])
            elif '_symmetry_space_group_name_h-m' in line:
                self.space_group = line.split()[1]
            elif 'loop_' in line:
                # loops for _atom_site, _symmetry and _geom
                heads = []
                body = []
                while line.startswith('_') or 'loop_' in line:
                    # must keep the loop_ line as this can still
                    # contain headers
                    heads.extend(line.split())
                    idx += 1
                    # don't lower these to keep atomic symbols
                    line = cif_file[idx].strip()
                while idx < len(cif_file) and not line.startswith('_') and not 'loop_' in line:
                    # shlex keeps 'quoted items' as one
                    # Some cifs seem to have primed atom symbols
                    # posix=False should help
                    # using .shlex instead of .split works with '#'
                    # comments too
                    split_line = shlex.shlex(line, posix=False)
                    split_line.whitespace_split = True
                    split_line = list(split_line)
                    body.extend([x.strip("'").strip('"') for x in split_line])
                    idx += 1
                    try:
                        line = cif_file[idx]
                    except IndexError:
                        line = ''
                if 'loop_' in heads:
                    heads.remove('loop_')
                loops.append((heads, body))
                continue
            idx += 1

        # cell first
        if np.all(params):
            self.cell.params = params
        else:
            error("No cell or incomplete cell found in cif file")

        # parse loop contents
        for heads, body in loops:
            if '_atom_site_fract_x' in heads:
                while body:
                    atoms.append(dict(zip(heads, body)))
                    body = body[len(heads):]
            if '_symmetry_equiv_pos_as_xyz' in heads:
                while body:
                    sym_dict = dict(zip(heads, body))
                    symmetry.append(
                        Symmetry(sym_dict['_symmetry_equiv_pos_as_xyz']))
                    body = body[len(heads):]
            if '_ccdc_geom_bond_type' in heads:
                while body:
                    bond_dict = dict(zip(heads, body))
                    # bond is sorted so there are no duplicates
                    # and tuple so it can be hashed
                    bond = (bond_dict['_geom_bond_atom_site_label_1'],
                            bond_dict['_geom_bond_atom_site_label_2'])
                    bond = tuple(sorted(bond))
                    # bond distance and type defualts to None if not specified
                    distance = bond_dict.get('_geom_bond_distance')
                    if distance is not None:
                        distance = float(distance)
                    bond_type = bond_dict.get('_ccdc_geom_bond_type')
                    cif_bonds[bond] = (distance, bond_type)
                    body = body[len(heads):]

        if not symmetry:
            debug('No symmetry found; assuming identity only')
            symmetry = [Symmetry('x,y,z')]

        newatoms = []
        for site_idx, atom in enumerate(atoms):
            for sym_op in symmetry:
                newatom = Atom(parent=self)
                newatom.from_cif(atom, self.cell.cell, sym_op, site_idx)
                newatoms.append(newatom)

        self.atoms = newatoms

        if len(symmetry) > 1:
            # can skip if just identity operation as it's slow for big systems
            # Some of pete's symmetrised mofs need a higher tolerence
            duplicate_tolerance = 0.2  # Angstroms
            self.remove_duplicates(duplicate_tolerance)

        bonds = {}
        # TODO(tdaff): this works for the one tested MOF; 0.1 was not enough
        # only check for bonds that are too long, not too short.
        bond_tolerence = 0.25
        # Assign bonds by index
        for bond, bond_data in cif_bonds.items():
            for first_index, first_atom in enumerate(self.atoms):
                if first_atom.site == bond[0]:
                    for second_index, second_atom in enumerate(self.atoms):
                        if second_atom is first_atom:
                            continue
                        elif second_atom.site == bond[1]:
                            # TODO(tdaff): symmetry implementation for cif bonding
                            distance = min_distance(first_atom, second_atom)
                            bond_dist = bond_data[0]
                            if bond_dist is None:
                                bond_dist = (first_atom.covalent_radius +
                                             second_atom.covalent_radius)
                            if distance < 0.6 * bond_dist:
                                warning("Short contact ignored: "
                                        "%s(%i) and %s(%i) = %.2f A" %
                                        (first_atom.site, first_index,
                                         second_atom.site, second_index,
                                         distance))
                            elif distance < (bond_dist + bond_tolerence):
                                # use the sorted index as bonds between the
                                # same type are doubly specified
                                bond_id = tuple(sorted((first_index, second_index)))
                                bonds[bond_id] = (distance, CCDC_BOND_ORDERS[bond_data[1]])
                                if first_atom.is_metal or second_atom.is_metal:
                                    first_atom.is_fixed = True
                                    second_atom.is_fixed = True

        self.bonds = bonds
        self.symmetry = symmetry

    def remove_duplicates(self, tolerance=0.02):
        """Find overlapping atoms and remove them."""
        uniq_atoms = []
        for atom in self.atoms:
            for uniq_atom in uniq_atoms:
                if atom.type != uniq_atom.type:
                    continue
                elif min_distance(atom, uniq_atom) < tolerance:
                    break
            # else excutes when not found here
            else:
                uniq_atoms.append(atom)
        debug("Found %i unique atoms in %i" % (len(uniq_atoms), self.natoms))
        self.atoms = uniq_atoms

    def check_close_contacts(self, absolute=1.0, covalent=None):
        """
        Check for atoms that are too close. Specify either an absolute distance
        in Angstrom or a scale factor for the sum of covalent radii. If a
        covalent factor is specified it will take priority over an absolute
        distance. Return True if close contacts found, else return False.
        """
        close_contact_found = False
        for atom_idx, atom in enumerate(self.atoms):
            for other_idx, other in enumerate(self.atoms):
                if other_idx >= atom_idx:
                    # short circuit half the calculations
                    # Can we do combinations with idx in 2.7
                    break
                if covalent is not None:
                    tolerance = covalent * (atom.covalent_radius +
                                            other.covalent_radius)
                else:
                    tolerance = absolute
                if min_distance(atom, other) < tolerance:
                    bond_ids = tuple(sorted([atom_idx, other_idx]))
                    if bond_ids not in self.bonds:
                        warning("Close atoms: %s(%i) and %s(%i)" %
                                (atom.site, atom_idx, other.site, other_idx))
                        close_contact_found = True

        return close_contact_found

    def bond_length_check(self, too_long=1.25, too_short=0.7):
        """
        Check if all bonds fall within a sensible range of scale factors
        of the sum of the covalent radii. Return True if bad bonds are found,
        otherwise False.

        """
        bad_bonds = False
        for bond in self.bonds:
            atom = self.atoms[bond[0]]
            other = self.atoms[bond[1]]
            distance = min_distance(atom, other)
            bond_dist = (atom.covalent_radius + other.covalent_radius)
            if distance > bond_dist * too_long:
                warning("Long bond found: %s(%i) and %s(%i) = %.2f A" %
                        (atom.site, bond[0], other.site, bond[1], distance))
                bad_bonds = True
            elif distance < bond_dist * too_short:
                warning("Short bond found: %s(%i) and %s(%i) = %.2f A" %
                        (atom.site, bond[0], other.site, bond[1], distance))
                bad_bonds = True

        return bad_bonds

    def gen_neighbour_list(self, force=False):
        """All atom pair distances."""
        # This can be expensive so skip if already calcualted
        if not force:
            for atom in self.atoms:
                if not hasattr(atom, 'neighbours'):
                    break
            else:
                # finished loop over all atoms
                debug("Neighbour list already calculated")
                return

        debug("Calculating neighbour list.")

        cell = self.cell.cell
        inv_cell = self.cell.inverse
        cpositions = [atom.ipos(cell, inv_cell) for atom in self.atoms]
        fpositions = [atom.ifpos(inv_cell) for atom in self.atoms]
        cell = cell.tolist()

        # loop over all pairs to find minimum periodic distances
        for atom, a_cpos, a_fpos in zip(self.atoms, cpositions, fpositions):
            neighbours = []
            for o_idx, o_cpos, o_fpos in zip(count(), cpositions, fpositions):
                sep = min_dist(a_cpos, a_fpos, o_cpos, o_fpos, cell)
                neighbours.append((sep, o_idx))
            # First one is self == 0
            # save in incresaing distance order
            atom.neighbours = sorted(neighbours)[1:]


    @property
    def types(self):
        """Ordered list of atom types."""
        return [atom.type for atom in self.atoms]

    @property
    def atomic_numbers(self):
        """Ordered list of atomic numbers."""
        return [atom.atomic_number for atom in self.atoms]

    @property
    def weight(self):
        """Unit cell weight."""
        return sum([atom.mass for atom in self.atoms])

    @property
    def volume(self):
        """Unit cell volume."""
        return self.cell.volume

    @property
    def natoms(self):
        """Number of atoms in the unit cell."""
        return len(self.atoms)


class Cell(object):
    """
    Crystollagraphic cell representations and interconversion methods.

    Setter methods can be defined for different file types, however
    .cell and .params will be self-consistent if set directly.

    """

    def __init__(self):
        """Default to a 1A cubic box."""
        self._cell = array([[1.0, 0.0, 0.0],
                            [0.0, 1.0, 0.0],
                            [0.0, 0.0, 1.0]])
        self._params = (1.0, 1.0, 1.0, 90.0, 90.0, 90.0)
        self._inverse = None

    @property
    def crystal_system(self):
        """Return the IUCr designation for the crystal system."""
        #FIXME(tdaff): must be aligned with x to work
        if self.alpha == self.beta == self.gamma == 90:
            if self.a == self.b == self.c:
                return 'cubic'
            elif self.a == self.b or self.a == self.c or self.b == self.c:
                return 'tetragonal'
            else:
                return 'orthorhombic'
        elif self.alpha == self.beta == 90:
            if self.a == self.b and self.gamma == 120:
                return 'hexagonal'
            else:
                return 'monoclinic'
        elif self.alpha == self.gamma == 90:
            if self.a == self.c and self.beta == 120:
                return 'hexagonal'
            else:
                return 'monoclinic'
        elif self.beta == self.gamma == 90:
            if self.b == self.c and self.alpha == 120:
                return 'hexagonal'
            else:
                return 'monoclinic'
        elif self.a == self.b == self.c and self.alpha == self.beta == self.gamma:
            return 'trigonal'
        else:
            return 'triclinic'

    @property
    def volume(self):
        """Calculate cell volume a.bxc."""
        b_cross_c = cross(self.cell[1], self.cell[2])
        return dot(self.cell[0], b_cross_c)

    def get_cell(self):
        """Get the 3x3 vector cell representation."""
        return self._cell

    def set_cell(self, value):
        """Set cell and params from the cell representation."""
        # Class internally expects an array
        self._cell = array(value).reshape((3, 3))
        self.__mkparam()
        self._inverse = np.linalg.inv(self.cell.T)

    # Property so that params are updated when cell is set
    cell = property(get_cell, set_cell)

    def get_params(self):
        """Get the six parameter cell representation as a tuple."""
        return tuple(self._params)

    def set_params(self, value):
        """Set cell and params from the cell parameters."""
        self._params = value
        self.__mkcell()
        self._inverse = np.linalg.inv(self.cell.T)

    # Property so that cell is updated when params are set
    params = property(get_params, set_params)

    @property
    def inverse(self):
        """Inverted cell matrix for converting to fractional coordinates."""
        try:
            if self._inverse is None:
                self._inverse = np.linalg.inv(self.cell.T)
        except AttributeError:
            self._inverse = np.linalg.inv(self.cell.T)
        return self._inverse

    @property
    def a(self):
        """Magnitude of cell a vector."""
        return self.params[0]

    @property
    def b(self):
        """Magnitude of cell b vector."""
        return self.params[1]

    @property
    def c(self):
        """Magnitude of cell c vector."""
        return self.params[2]

    @property
    def alpha(self):
        """Cell angle alpha."""
        return self.params[3]

    @property
    def beta(self):
        """Cell angle beta."""
        return self.params[4]

    @property
    def gamma(self):
        """Cell angle gamma."""
        return self.params[5]

    # Implementation details -- directly access the private _{cell|param}
    # attributes; please don't break.
    def __mkcell(self):
        """Update the cell representation to match the parameters."""
        a_mag, b_mag, c_mag = self.params[:3]
        alpha, beta, gamma = [x * DEG2RAD for x in self.params[3:]]
        a_vec = array([a_mag, 0.0, 0.0])
        b_vec = array([b_mag * cos(gamma), b_mag * sin(gamma), 0.0])
        c_x = c_mag * cos(beta)
        c_y = c_mag * (cos(alpha) - cos(gamma) * cos(beta)) / sin(gamma)
        c_vec = array([c_x, c_y, (c_mag**2 - c_x**2 - c_y**2)**0.5])
        self._cell = array([a_vec, b_vec, c_vec])

    def __mkparam(self):
        """Update the parameters to match the cell."""
        cell_a = sqrt(sum(x**2 for x in self.cell[0]))
        cell_b = sqrt(sum(x**2 for x in self.cell[1]))
        cell_c = sqrt(sum(x**2 for x in self.cell[2]))
        alpha = arccos(sum(self.cell[1, :] * self.cell[2, :]) /
                       (cell_b * cell_c)) * 180 / pi
        beta = arccos(sum(self.cell[0, :] * self.cell[2, :]) /
                      (cell_a * cell_c)) * 180 / pi
        gamma = arccos(sum(self.cell[0, :] * self.cell[1, :]) /
                       (cell_a * cell_b)) * 180 / pi
        self._params = (cell_a, cell_b, cell_c, alpha, beta, gamma)


class Atom(object):
    """Base atom object."""

    def __init__(self, at_type=False, pos=False, parent=None, **kwargs):
        """Accept arbritary kwargs as attributes."""
        if parent is not None:
            self._parent = parent
        self.type = at_type
        self.pos = pos
        self.charge = 0.0
        self.idx = False
        self.site = None
        self.mass = 0.0
        self.molecule = None
        self.uff_type = None
        self.is_fixed = False
        # Sets anything else specified as an attribute
        for key, val in kwargs.items():
            setattr(self, key, val)

    def __str__(self):
        return "%s %f %f %f" % tuple([self.type] + list(self.pos))

    def __repr__(self):
        return "Atom(%r,%r)" % (self.type, self.pos)

    def fpos(self, inv_cell):
        """Fractional position within a given cell."""
        return dot(inv_cell, self.pos).tolist()

    def ifpos(self, inv_cell):
        """In cell fractional position."""
        return [i % 1 for i in self.fpos(inv_cell)]

    def ipos(self, cell, inv_cell):
        """In cell cartesian position."""
        return dot(self.ifpos(inv_cell), cell)

    def from_cif(self, at_dict, cell, symmetry=None, idx=None):
        """Extract an atom description from dictionary of cif items."""
        self.site = at_dict['_atom_site_label']
        # type_symbol takes precedence but need not be specified
        self.type = at_dict.get('_atom_site_type_symbol', self.site)
        self.mass = WEIGHT[self.type]
        frac_pos = [ufloat(at_dict['_atom_site_fract_x']),
                    ufloat(at_dict['_atom_site_fract_y']),
                    ufloat(at_dict['_atom_site_fract_z'])]
        if symmetry is not None:
            frac_pos = symmetry.trans_frac(frac_pos)
        self.pos = dot(frac_pos, cell)
        if re.match('[0-9]', self.site) and idx is not None:
            debug("Site label may not be unique; appending index")
            self.site = "%s%i" % (self.site, idx)
        if '_atom_site_description' in at_dict:
            self.uff_type = at_dict['_atom_site_description']
        elif '_atom_type_description' in at_dict:
            self.uff_type = at_dict['_atom_type_description']
        if '_atom_type_partial_charge' in at_dict:
            self.charge = float(at_dict['_atom_type_partial_charge'])
        elif '_atom_type_parital_charge' in at_dict:
            #TODO(tdaff): remove for 2.0
            self.charge = float(at_dict['_atom_type_parital_charge'])

    def from_pdb(self, line, charges=False):
        """
        Parse the ATOM line from a pdb file.
        Occupancy field may be used to specify the charge as in a '.pqr' file.
        """
        # pdb is defined with fixed width fields rather than splitting
        self.idx = try_int(line[6:11])
        self.site = line[12:16].strip()
        self.molecule = try_int(line[22:26])
        at_pos = float(line[30:38]), float(line[38:46]), float(line[47:54])
        self.pos = at_pos
        self.type = line[76:78].strip()
        self.mass = WEIGHT[self.type]
        if charges:
            self.charge = float(line[54:60])

    def from_vasp(self, line, at_type=None, cell=identity(3)):
        """Set the atom data from vasp input. Only pass cell if fractional."""
        self.pos = dot([float(x) for x in line.split()[:3]], cell)
        if at_type is not None:
            self.type = at_type
            self.mass = WEIGHT[at_type]

    def from_siesta(self, line, cell):
        """Parse line from SIESTA.STRUCT_OUT file."""
        self.pos = dot([float(x) for x in line.split()[2:5]], cell)
        self.atomic_number = int(line.split()[1])
        self.mass = WEIGHT[self.type]

    def from_xyz(self, line):
        """Parse line from generic xyz file."""
        split_line = line.split()
        self.pos = [float(x) for x in split_line[1:4]]
        if len(split_line) > 4:
            try:
                self.charge = float(split_line[4])
            except ValueError:
                # Assume a comment so skip it
                pass
        self.type = line.split()[0]
        self.mass = WEIGHT[self.type]

    def translate(self, vec):
        """Move the atom by the given vector."""
        self.pos = [x + y for x, y in zip(self.pos, vec)]

    def get_atomic_number(self):
        """The atomic number for the element, or closest match."""
        name = self.type
        atomic_number = None
        while name:
            try:
                atomic_number = ATOMIC_NUMBER.index(name)
                break
            except ValueError:
                name = name[:-1]
        return atomic_number

    def set_atomic_number(self, value):
        """Set the atom type based on the atomic number."""
        self.type = ATOMIC_NUMBER[value]

    atomic_number = property(get_atomic_number, set_atomic_number)

    def get_fractional_coordinate(self):
        """Retrieve the fractional coordinates or calculate from the parent."""
        try:
            # Hopefully the attribute is just set
            return self._fractional
        except AttributeError:
            # Not set yet
            try:
                self._fractional = self.fpos(self._parent.cell.inverse)
                return self._fractional
            except AttributeError:
                return None

    def set_fractional_coordinate(self, value):
        """Set the position using the fractional coordinates."""
        fractional = [x % 1.0 for x in value]
        self._fractional = fractional
        self.pos = dot(fractional, self._parent.cell.cell)

    def del_fractional_coordinate(self):
        """Remove the fractional coordinate; run after updating cell."""
        try:
            del self._fractional
        except AttributeError:
            # Attribute does not exist, so can't delete it
            pass

    fractional = property(get_fractional_coordinate, set_fractional_coordinate, del_fractional_coordinate)

    @property
    def element(self):
        """Guess the element from the type, fall back to type."""
        name = self.type
        while name:
            if name in ATOMIC_NUMBER:
                return name
            else:
                name = name[:-1]
        return self.type

    @property
    def vdw_radius(self):
        """Get the vdw radius from the UFF parameters."""
        return UFF[self.type][0]/2.0

    @property
    def covalent_radius(self):
        """Get the covalent radius from the library parameters."""
        if self.type == 'C' and self.uff_type:
            return COVALENT_RADII[self.uff_type]
        return COVALENT_RADII[self.type]

    @property
    def is_metal(self):
        """Return True if element is in a predetermined set of metals."""
        return self.atomic_number in METALS



class Symmetry(object):
    """Apply symmetry operations to atomic coordinates."""
    def __init__(self, blob):
        """Read the operation from the argument."""
        self.sym_ops = []
        self.blob = blob
        self.parse_blob()
        self.cell = None

    def parse_blob(self):
        """Interpret a symmetry line from a cif."""
        # convert integers to floats to avoid integer division
        self.sym_ops = [re.sub(r'([\d]+)', r'\1.0', x.strip())
                        for x in re.split(',', self.blob) if x.strip()]

    def trans_frac(self, pos):
        """Apply symmetry operation to the supplied position."""
        new_pos = [eval(sym_op.replace('x', str(pos[0]))
                        .replace('y', str(pos[1]))
                        .replace('z', str(pos[2]))) for sym_op in self.sym_ops]
        # TODO(tdaff): should we translate into cell?
        new_pos = [x%1.0 for x in new_pos]
        return new_pos




def unique(in_list, key=None):
    """Unique values in list ordered by first occurance"""
    uniq = []
    if key is not None:
        keys = []
        for item in in_list:
            item_key = key(item)
            if item_key not in keys:
                uniq.append(item)
                keys.append(item_key)
    else:
        for item in in_list:
            if item not in uniq:
                uniq.append(item)
    return uniq



# General utility functions
def mkdirs(directory):
    """Create a directory if it does not exist."""
    if not path.exists(directory):
        os.makedirs(directory)



def move_and_overwrite(src, dest):
    """Move src to dest and overwrite if it is an existing file."""
    # As src and dest can be files or directories, do some checks.
    if path.exists(dest):
        if path.isdir(dest):
            dest_full = path.join(dest, path.basename(src))
            if path.exists(dest_full):
                if path.isfile(dest_full):
                    os.remove(dest_full)
                    shutil.move(src, dest)
                else:
                    raise OSError("Directory %s already exists" % dest_full)
            else:
                shutil.move(src, dest)
        elif path.isfile(dest):
            os.remove(dest)
            shutil.move(src, dest)
        else:
            raise OSError("%s is not a folder or file" % dest)
    else:
        shutil.move(src, dest)



def ufloat(text):
    """Convert string to float, ignoring the uncertainty part."""
    return float(re.sub('\(.*\)', '', text))



def try_int(text, default=0):
    """Try to parse an integer but return a default if it fails."""
    try:
        return int(text)
    except ValueError:
        return default


def prod(seq):
    """Calculate the product of all members of a sequence."""
    # numpy.prod will silently overflow 32 bit integer values
    # so we can use python bignums natively
    product = 1
    for item in seq:
        product *= item
    return product


def dot3(vec1, vec2):
    """Calculate dot product for two 3d vectors."""
    return vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2]


def min_distance(first_atom, second_atom, cell=None):
    """Helper to find mimimum image criterion distance."""
    if cell is None:
        cell = first_atom._parent.cell.cell
    return min_dist(first_atom.pos,
                    first_atom.fractional,
                    second_atom.pos,
                    second_atom.fractional,
                    cell)


def min_dist(c_coa, f_coa, c_cob, f_cob_in, box):
    """Calculate the closest distance assuming fractional, in-cell coords."""
    f_cob = f_cob_in[:]
    fdx = f_coa[0] - f_cob[0]
    if fdx < -0.5:
        f_cob[0] -= 1
    elif fdx > 0.5:
        f_cob[0] += 1
    fdy = f_coa[1] - f_cob[1]
    if fdy < -0.5:
        f_cob[1] -= 1
    elif fdy > 0.5:
        f_cob[1] += 1
    fdz = f_coa[2] - f_cob[2]
    if fdz < -0.5:
        f_cob[2] -= 1
    elif fdz > 0.5:
        f_cob[2] += 1
    if f_cob == f_cob_in:
        # if nothing has changed, use initial values
        return vecdist3(c_coa, c_cob)
    else:
        new_b = [f_cob[0]*box[0][0] + f_cob[1]*box[1][0] + f_cob[2]*box[2][0],
                 f_cob[0]*box[0][1] + f_cob[1]*box[1][1] + f_cob[2]*box[2][1],
                 f_cob[0]*box[0][2] + f_cob[1]*box[1][2] + f_cob[2]*box[2][2]]
        return vecdist3(c_coa, new_b)


def vecdist3(coord1, coord2):
    """Calculate vector between two 3d points."""
    #return [i - j for i, j in zip(coord1, coord2)]
    # Twice as fast for fixed 3d vectors
    vec = [coord2[0] - coord1[0],
           coord2[1] - coord1[1],
           coord2[2] - coord1[2]]

    return (vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2])**0.5




def strip_blanks(lines):
    """Strip lines and remove blank lines."""
    return [line.strip() for line in lines if line.strip() != '']


def subgroup(iterable, width, itype=None):
    """Split an iterable into nested sub-itypes of width members."""
    # Return the same type as iterable
    if itype is None:
        if isinstance(iterable, list):
            itype = list
        else:
            itype = tuple
    # Will leave short groups if not enough members
    return itype([itype(iterable[x:x+width])
                  for x in range(0, len(iterable), width)])


def name_from_types(sites, guest):
    """Generate a string that gives the atoms represented in sites."""
    stypes = []
    for site_ident in sites:
        if site_ident is 0:
            stypes.append('COM')
        else:
            stypes.append(guest.atoms[site_ident-1].element)
    stypes = unique(stypes)
    if len(stypes) > 1:
        site_name = "-".join(stypes)
    else:
        site_name = stypes[0]
    return site_name


def other_bond_index(bond, index):
    """Return the atom index for the other atom in a bond."""
    if bond[0] == index:
        return bond[1]
    elif bond[1] == index:
        return bond[0]
    else:
        raise ValueError("Index %s not found in bond %s" % (index, bond))






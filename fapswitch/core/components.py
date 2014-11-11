
from __future__ import division

"""
fapswitch core components

Includes structures and atoms

"""

import os
import re
import shlex
from itertools import count
from logging import warning, debug, error, info
from os import path

import numpy as np
from numpy import pi, cos, sin, sqrt, arccos
from numpy import array, dot, cross
from numpy.linalg import norm

from fapswitch.core.util import min_vect, normalise, strip_blanks, vecdist3
from fapswitch.core.util import arbitrary_normal, min_vector
from fapswitch.core.elements import WEIGHT, ATOMIC_NUMBER, UFF
from fapswitch.core.elements import CCDC_BOND_ORDERS, METALS, OB_BOND_ORDERS
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
        self.bonds = {}
        self.attachments = {}
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
            # Each atom needs to know which other atom it was generated from
            # and by which symmetry operations.
            asymmetric_site = None
            for sym_op in symmetry:
                newatom = Atom(parent=self)
                newatom.from_cif(atom, self.cell.cell, sym_op, site_idx)
                if asymmetric_site is None:
                    asymmetric_site = newatom
                else:
                    newatom.asymmetric_site = asymmetric_site
                newatoms.append(newatom)

        self.atoms = newatoms

        if len(symmetry) > 1:
            # can skip if just identity operation as it's slow for big systems
            # Some of pete's symmetrised mofs need a higher tolerence
            duplicate_tolerance = 0.2  # Angstroms
            self.remove_duplicates(duplicate_tolerance)

        bonds = {}
        # Only check for bonds that are too long, not too short.
        # Tolerence was tested on some difficult MOFs, needed to
        # be > 0.1
        bond_tolerance = 0.25
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
                            elif distance < (bond_dist + bond_tolerance):
                                # use the sorted index as bonds between the
                                # same type are doubly specified
                                bond_id = tuple(sorted((first_index, second_index)))
                                bonds[bond_id] = (distance, CCDC_BOND_ORDERS[bond_data[1]])

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
                    uniq_atom.sym_ops.extend(atom.sym_ops)
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

    def gen_factional_positions(self):
        """
        Precalculate the fractional positions for all the atoms in the
        current cell. These will be incorrect if the cell chagnes!

        """

        cell = self.cell.cell
        inv_cell = self.cell.inverse
        for atom in self.atoms:
            atom.cellpos = atom.ipos(cell, inv_cell)
            atom.cellfpos = atom.ifpos(inv_cell)

    def gen_normals(self):
        """
        Calculate the normal to the vectors between the first two neighbours.

        """
        self.gen_neighbour_list()
        broken_symmetry = False

        cell = self.cell.cell
        inv_cell = self.cell.inverse
        cpositions = [atom.ipos(cell, inv_cell) for atom in self.atoms]
        fpositions = [atom.ifpos(inv_cell) for atom in self.atoms]

        for at_idx, atom in enumerate(self.atoms):
            # If we have some symmetry operations, use them to get the
            # other normals to keep as much symmetry as possible
            if atom.asymmetric_site is not None:
                ref = atom.asymmetric_site
                normals = []
                for sym_op in atom.sym_ops:
                    # Place an atom at the normal, and use the symmetry
                    # operation to find the normal rotation
                    normal = Atom(parent=self)
                    normal.pos = ref.pos + ref.normal
                    normal.fractional = sym_op.trans_frac(normal.fractional)
                    unit_normal = normalise(min_vector(normal, atom))
                    for other_normal in normals:
                        if np.allclose(unit_normal, other_normal):
                            break
                    else:
                        normals.append(unit_normal)
                if len(normals) > 1:
                    print(normals)
                    broken_symmetry = True
                    atom.normals = normals
                # Only use the first operation as it will be the first for
                # all atoms.
                atom.normal = normals[0]
            else:
                # Speed is not so important here but these are derived from the
                # surface ares distance calculations so they are explicitly
                # given the fractional positions.
                conex = [x[1] for x in atom.neighbours[:2]]
                left = min_vect(cpositions[at_idx], fpositions[at_idx],
                                cpositions[conex[0]], fpositions[conex[0]],
                                cell)
                right = min_vect(cpositions[at_idx], fpositions[at_idx],
                                 cpositions[conex[1]], fpositions[conex[1]],
                                 cell)
                normal = cross(left, right)
                if norm(normal) == 0.0:
                    normal = arbitrary_normal(left)
                    # already normalised
                    atom.normal = normal
                else:
                    atom.normal = normalise(normal)

        if broken_symmetry:
            info("Symmetry may be broken by asymmetric functional groups.")
        else:
            info("Symmetry maintained with asymmetric functional groups.")

    def gen_attachment_sites(self):
        """
        Find connected atoms (bonds) and organise by sites. Symmetry
        equivalent atoms have the same site.

        """

        atoms = self.atoms
        cell = self.cell.cell
        inv_cell = self.cell.inverse
        self.gen_neighbour_list()
        self.attachments = {}

        for at_idx, atom in enumerate(self.atoms):
            if atom.type == 'H':
                h_index = self.atoms.index(atom)
                for bond in self.bonds:
                    if h_index in bond:
                        o_index = bond[0] if bond[1] == h_index else bond[1]
                        break
                else:
                    error("Unbonded hydrogen")

                # Generate the direction here as we were getting buggy
                # attachment over periodic boundaries
                direction = min_vect(atoms[o_index].ipos(cell, inv_cell),
                                     atoms[o_index].ifpos(inv_cell),
                                     atom.ipos(cell, inv_cell),
                                     atom.ifpos(inv_cell), cell)
                direction = normalise(direction)

                if atom.site in self.attachments:
                    self.attachments[atom.site].append((h_index, o_index, direction))
                else:
                    self.attachments[atom.site] = [(h_index, o_index, direction)]

    def gen_types_from_bonds(self):
        """
        Pass the bonding information into openbabel to get the atomic types.
        """
        # import these locally so we can run fapswitch without them
        import openbabel as ob
        import pybel

        # Construct the molecule from atoms and bonds

        obmol = ob.OBMol()
        obmol.BeginModify()

        for atom in self.atoms:
            new_atom = obmol.NewAtom()
            new_atom.SetAtomicNum(atom.atomic_number)

        for bond, bond_info in self.bonds.items():
            # Remember openbabel indexes from 1
            obmol.AddBond(bond[0]+1, bond[1]+1, OB_BOND_ORDERS[bond_info[1]])

        obmol.EndModify()

        pybel_mol = pybel.Molecule(obmol)

        # need to tell the typing system to ignore all atoms in the setup
        # or it will silently crash with memory issues
        constraint = ob.OBFFConstraints()
        for at_idx in range(pybel_mol.OBMol.NumAtoms()):
            constraint.AddIgnore(at_idx)
        uff = ob.OBForceField_FindForceField('uff')
        uff.Setup(pybel_mol.OBMol, constraint)
        uff.GetAtomTypes(pybel_mol.OBMol)
        # Dative nitrogen bonds break aromaticity determination from resonant
        # structures, so make anything with an aromatic bond be aromatic
        for at_idx, atom, ob_atom in zip(count(), self.atoms, pybel_mol):
            uff_type = ob_atom.OBAtom.GetData("FFAtomType").GetValue()
            if atom.type in ['C', 'N', 'O', 'S']:
                for bond, bond_info in self.bonds.items():
                    if at_idx in bond and bond_info[1] == 1.5:
                        uff_type = uff_type[0] + '_R'
                        break

            atom.uff_type = uff_type

    def gen_babel_uff_properties(self):
        """
        Process a supercell with babel to calculate UFF atom types and
        bond orders.
        """
        # import these locally so we can run fapswitch without them
        import openbabel as ob
        import pybel
        # Pass as free form fractional
        # GG's periodic should take care of bonds
        warning("Assuming periodic openbabel; not generating supercell")
        cell = self.cell
        as_fffract = ['generated fractionals\n', '%f %f %f %f %f %f\n' % cell.params]
        for atom in self.atoms:
            atom_line = ("%s " % atom.type + "%f %f %f\n" % tuple(atom.cellfpos))
            as_fffract.append(atom_line)
        pybel_string = ''.join(as_fffract)
        pybel_mol = pybel.readstring('fract', pybel_string)
        # need to tell the typing system to ignore all atoms in the setup
        # or it will silently crash with memory issues
        constraint = ob.OBFFConstraints()
        for at_idx in range(pybel_mol.OBMol.NumAtoms()):
            constraint.AddIgnore(at_idx)
        uff = ob.OBForceField_FindForceField('uff')
        uff.Setup(pybel_mol.OBMol, constraint)
        uff.GetAtomTypes(pybel_mol.OBMol)
        for atom, ob_atom in zip(self.atoms, pybel_mol):
            atom.uff_type = ob_atom.OBAtom.GetData("FFAtomType").GetValue()

        bonds = {}
        max_idx = self.natoms
        # look at all the bonds separately from the atoms
        for bond in ob.OBMolBondIter(pybel_mol.OBMol):
            # These rules are translated from ob/forcefielduff.cpp...
            start_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            if start_idx > max_idx and end_idx > max_idx:
                continue
            if end_idx > max_idx:
                end_idx %= max_idx
            if start_idx > max_idx:
                start_idx %= max_idx

            start_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()

            bond_order = bond.GetBondOrder()
            if bond.IsAromatic():
                bond_order = 1.5
            # e.g., in Cp rings, may not be "aromatic" by OB
            # but check for explicit hydrogen counts
            #(e.g., biphenyl inter-ring is not aromatic)
            #FIXME(tdaff): aromatic C from GetType is "Car" is this correct?
            if start_atom.GetType()[-1] == 'R' and end_atom.GetType()[-1] == 'R' and start_atom.ExplicitHydrogenCount() == 1 and end_atom.ExplicitHydrogenCount() == 1:
                bond_order = 1.5
            if bond.IsAmide():
                bond_order = 1.41
            # save the indicies as zero based
            bond_length = bond.GetLength()
            bond_id = tuple(sorted((start_idx-1, end_idx-1)))
            bonds[bond_id] = (bond_length, bond_order)

        self.bonds = bonds

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
        self.uff_type = None
        self.vdw_radius = 0.0
        self.covalent_radius = 0.0
        self.is_metal = False
        self.sym_ops = [Symmetry('x,y,z')]
        self.asymmetric_site = None

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
            self.sym_ops = [symmetry]
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
            self.charge = float(at_dict['_atom_type_parital_charge'])

        # set these as atributes, not properties
        # breaks use with other initialisations
        self.vdw_radius = UFF[self.type][0]/2.0
        if self.type == 'C' and self.uff_type:
            self.covalent_radius = COVALENT_RADII[self.uff_type]
        else:
            self.covalent_radius = COVALENT_RADII[self.type]
        self.is_metal = self.atomic_number in METALS

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

    fractional = property(get_fractional_coordinate,
                          set_fractional_coordinate,
                          del_fractional_coordinate)

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


class Symmetry(object):
    """Apply symmetry operations to atomic coordinates."""
    def __init__(self, text):
        """Read the operation from the argument."""
        self.xyz = "".join(text.split())  # remove all spaces
        rt_matrix = []
        for idx, sub_xyz in enumerate(self.xyz.split(',')):
            row = [0.0, 0.0, 0.0, 0.0]
            # Find multipliers of x, y and z for rotation
            for component in re.finditer(r'([+-]?)(\d*)([x-z]+)', sub_xyz):
                # signed?
                value = -1.0 if component.group(1) == '-' else 1.0
                # Some scaling factor?
                if component.group(2):
                    value *= int(component.group(2))
                if component.group(3) == 'x':
                    row[0] += value
                elif component.group(3) == 'y':
                    row[1] += value
                elif component.group(3) == 'z':
                    row[2] += value

            # Find the translation component any trailing numbers or fractions
            for component in re.finditer(r"([+-])(\d+)/*(\d*)$", sub_xyz):
                # signed?
                value = -1.0 if component.group(1) == '-' else 1.0
                value *= int(component.group(2))
                if component.group(3):
                    value /= int(component.group(3))
                row[3] += value

            rt_matrix.append(row)

        rt_matrix.append([0.0, 0.0, 0.0, 1.0])
        self.rt_matrix = np.array(rt_matrix)

    def trans_frac(self, pos, in_cell=True):
        """Apply symmetry operation to the supplied position."""

        new_pos = dot(self.rt_matrix, [pos[0], pos[1], pos[2], 1.0])[:3]
        #print(translated, also_new)
        if in_cell:
            # translate positions into cell; leave as numpy array
            new_pos %= 1.0

        return new_pos

    def rotate(self, vector):
        """
        Apply only the rotation part of the
        symmetry operation to the vector.

        """
        return dot(self.rt_matrix[:3, :3], vector)

    def __repr__(self):
        return "Symmetry('{}')".format(self.xyz)


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


def ufloat(text):
    """Convert string to float, ignoring the uncertainty part."""
    return float(re.sub(r'\(.*\)', '', text))


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


def other_bond_index(bond, index):
    """Return the atom index for the other atom in a bond."""
    if bond[0] == index:
        return bond[1]
    elif bond[1] == index:
        return bond[0]
    else:
        raise ValueError("Index %s not found in bond %s" % (index, bond))

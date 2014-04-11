#!/usr/bin/env python

"""
create_functional_group.py

Generate a functional group block for a functional_groups.lib file. Takes a
smiles string and uses openbabel to deduce the 3D structure and aligns it
based on a benzene ring.

"""

import argparse

import openbabel as ob
import numpy as np
import pybel
from numpy import asarray, cross, dot, array, identity
from numpy.linalg import norm
from scipy.spatial.distance import cdist


ATOMIC_NUMBER = [
    "ZERO", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg",
    "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn",
    "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb",
    "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In",
    "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm",
    "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta",
    "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At",
    "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk",
    "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt",
    "Ds", "Rg", "Cn", "Uut", "Uuq", "Uup", "Uuh", "Uuo"]

UFF = {
    "H": (2.5711, 0.0440),
    "He": (2.1043, 0.0560),
    "Li": (2.1836, 0.0250),
    "Be": (2.4455, 0.0850),
    "B": (3.6375, 0.1800),
    "C": (3.4309, 0.1050),
    "N": (3.2607, 0.0690),
    "O": (3.1181, 0.0600),
    "F": (2.9970, 0.0500),
    "Ne": (2.8892, 0.0420),
    "Na": (2.6576, 0.0300),
    "Mg": (2.6914, 0.1110),
    "Al": (4.0082, 0.5050),
    "Si": (3.8264, 0.4020),
    "P": (3.6946, 0.3050),
    "S": (3.5948, 0.2740),
    "Cl": (3.5164, 0.2270),
    "Ar": (3.4460, 0.1850),
    "K": (3.3961, 0.0350),
    "Ca": (3.0282, 0.2380),
    "Sc": (2.9355, 0.0190),
    "Ti": (2.8286, 0.0170),
    "V": (2.8010, 0.0160),
    "Cr": (2.6932, 0.0150),
    "Mn": (2.6380, 0.0130),
    "Fe": (2.5943, 0.0130),
    "Co": (2.5587, 0.0140),
    "Ni": (2.5248, 0.0150),
    "Cu": (3.1137, 0.0050),
    "Zn": (2.4616, 0.1240),
    "Ga": (3.9048, 0.4150),
    "Ge": (3.8130, 0.3790),
    "As": (3.7685, 0.3090),
    "Se": (3.7462, 0.2910),
    "Br": (3.7320, 0.2510),
    "Kr": (3.6892, 0.2200),
    "Rb": (3.6652, 0.0400),
    "Sr": (3.2438, 0.2350),
    "Y": (2.9801, 0.0720),
    "Zr": (2.7832, 0.0690),
    "Nb": (2.8197, 0.0590),
    "Mo": (2.7190, 0.0560),
    "Tc": (2.6709, 0.0480),
    "Ru": (2.6397, 0.0560),
    "Rh": (2.6094, 0.0530),
    "Pd": (2.5827, 0.0480),
    "Ag": (2.8045, 0.0360),
    "Cd": (2.5373, 0.2280),
    "In": (3.9761, 0.5990),
    "Sn": (3.9128, 0.5670),
    "Sb": (3.9378, 0.4490),
    "Te": (3.9823, 0.3980),
    "I": (4.0090, 0.3390),
    "Xe": (3.9235, 0.3320),
    "Cs": (4.0242, 0.0450),
    "Ba": (3.2990, 0.3640),
    "La": (3.1377, 0.0170),
    "Ce": (3.1680, 0.0130),
    "Pr": (3.2126, 0.0100),
    "Nd": (3.1850, 0.0100),
    "Pm": (3.1600, 0.0090),
    "Sm": (3.1360, 0.0080),
    "Eu": (3.1119, 0.0080),
    "Gd": (3.0005, 0.0090),
    "Tb": (3.0745, 0.0070),
    "Dy": (3.0540, 0.0070),
    "Ho": (3.0371, 0.0070),
    "Er": (3.0210, 0.0070),
    "Tm": (3.0059, 0.0060),
    "Yb": (2.9890, 0.2280),
    "Lu": (3.2429, 0.0410),
    "Hf": (2.7983, 0.0720),
    "Ta": (2.8241, 0.0810),
    "W": (2.7342, 0.0670),
    "Re": (2.6317, 0.0660),
    "Os": (2.7796, 0.0370),
    "Ir": (2.5302, 0.0730),
    "Pt": (2.4535, 0.0800),
    "Au": (2.9337, 0.0390),
    "Hg": (2.4099, 0.3850),
    "Tl": (3.8727, 0.6800),
    "Pb": (3.8282, 0.6630),
    "Bi": (3.8932, 0.5180),
    "Po": (4.1952, 0.3250),
    "At": (4.2318, 0.2840),
    "Rn": (4.2451, 0.2480),
    "Fr": (4.3654, 0.0500),
    "Ra": (3.2758, 0.4040),
    "Ac": (3.0985, 0.0330),
    "Th": (3.0255, 0.0260),
    "Pa": (3.0504, 0.0220),
    "U": (3.0246, 0.0220),
    "Np": (3.0504, 0.0190),
    "Pu": (3.0504, 0.0160),
    "Am": (3.0121, 0.0140),
    "Cm": (2.9631, 0.0130),
    "Bk": (2.9747, 0.0130),
    "Cf": (2.9515, 0.0130),
    "Es": (2.9391, 0.0120),
    "Fm": (2.9275, 0.0120),
    "Md": (2.9168, 0.0110),
    "No": (2.8936, 0.0110),
    "Lr": (2.8829, 0.0110)
}


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
    h = (1 - c)/(np.dot(v, v))
    return array([[c + h*v[0]*v[0], h*v[0]*v[1] - v[2], h*v[0]*v[2] + v[1]],
                  [h*v[0]*v[1] + v[2], c + h*v[1]*v[1], h*v[1]*v[2] - v[0]],
                  [h*v[0]*v[2] - v[1], h*v[1]*v[2] + v[0], c + h*v[2]*v[2]]])


def realign(coordinates, origin_index, up_index, plane_index):
    """
    Take a set of coordinates and align them so that up_index lies along
    the [0, 1, 0] direction and plane index lies in the plane,
    returning the aligned coordinates as an array.
    """
    # move everything
    array_coords = asarray(coordinates)
    array_coords -= array_coords[origin_index]
    # align with axis
    current_axis = array_coords[up_index] - array_coords[origin_index]
    align_matrix = matrix_rotate(-current_axis, [0, 1, 0])
    array_coords = asarray([dot(align_matrix, x) for x in array_coords])
    # align with plane
    in_plane_vector = array_coords[plane_index]
    plane_rotate = matrix_rotate([in_plane_vector[0], 0, in_plane_vector[2]],
                                 [1, 0, 0])
    array_coords = asarray([dot(plane_rotate, x) for x in array_coords])
    return array_coords


def main():
    """
    Create an aligned functional group based on command line arguments.
    """

    parser = argparse.ArgumentParser(
        description='Create a functional group from a smiles pattern',
        epilog='Example usage: %(prog)s -s OPr -n PropylEther -c "Alkyl Ether" -m OCCC')

    parser.add_argument('smi_string', help="Smiles string to generate group")
    parser.add_argument('-s', '--short-name',
                        help='Short name (defaults to smiles string)')
    parser.add_argument('-n', '--name', required=True,
                        help='Descriptive name (e.g. PropylEther)')
    parser.add_argument('-m', '--mepo-compatible', action='store_true',
                        help='Record group as compatible with MEPO-QEq')
    parser.add_argument('-c', '--classification',
                        help='General classification (e.g. "Alkyl Halide")')
    parser.add_argument('-t', '--terminal', action='store_true',
                        help='Output to terminal as well as files')

    args = parser.parse_args()

    fgroup = args.smi_string
    if '%99' in fgroup:
        print('Do not use ring closure 99')
        raise SystemExit
    if not args.short_name:
        args.short_name = fgroup

    # Use an explicitly defined benzene as a base
    # Do rings closure at 99 in case functional group has other closures
    attached = '[cH]%99[cH][cH][cH][cH]c%99'

    # make3D by default gives an optimised structure, great!
    pybel_mol = pybel.readstring('smi', attached + fgroup)
    pybel_mol.title = "[{}] {}".format(args.short_name, args.name)
    pybel_mol.make3D(forcefield='UFF')

    uff = ob.OBForceField_FindForceField('uff')
    uff.Setup(pybel_mol.OBMol)
    uff.GetAtomTypes(pybel_mol.OBMol)

    coordinates = []

    for ob_atom in pybel_mol:
        coordinates.append(ob_atom.coords)

    rotated_coordinates = realign(coordinates, 11, 10, 8)

    bonds = {}

    # look at all the bonds separately from the atoms
    for bond in ob.OBMolBondIter(pybel_mol.OBMol):
        # These rules are translated from ob/forcefielduff.cpp...
        start_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()

        start_atom = bond.GetBeginAtom()
        end_atom = bond.GetEndAtom()

        bond_order = bond.GetBondOrder()
        if bond.IsAromatic():
            bond_order = 1.5
        # e.g., in Cp rings, may not be "aromatic" by OB
        # but check for explicit hydrogen counts
        #(e.g., biphenyl inter-ring is not aromatic)
        #FIXME(tdaff): aromatic C from GetType is "Car" is this correct?
        if (start_atom.GetType()[-1] == 'R' and
            end_atom.GetType()[-1] == 'R' and
            start_atom.ExplicitHydrogenCount() == 1
            and end_atom.ExplicitHydrogenCount() == 1):
            bond_order = 1.5
        if bond.IsAmide():
            bond_order = 1.41
        # Zero the indicies for the connecting atom so that
        # negative indexes are benzene atoms
        bond_length = bond.GetLength()
        bond_id = tuple(sorted((start_idx-12, end_idx-12)))
        bonds[bond_id] = (bond_length, bond_order)

    # We can start building our output now!

    output_text = [
        "[{}]\n".format(args.short_name),
        "name = {}\n".format(args.name),
        "smiles = {}\n".format(fgroup),
        "mepo_compatible = {}\n".format(args.mepo_compatible)]

    if args.classification:
        output_text.append("class = {}\n".format(args.classification))

    # functional group fingerprint
    nbins = 10
    max_distance = 10.0
    bin_width = max_distance/nbins
    fingerprint = [0.0]*(nbins*3)

    atom_block = []
    base_atom = pybel_mol.atoms[10].OBAtom
    for ob_atom, coord in zip(pybel_mol, rotated_coordinates):
        atom_idx = ob_atom.OBAtom.GetIndex()
        if atom_idx > 10:
            atomicnum = ob_atom.atomicnum
            element = ATOMIC_NUMBER[atomicnum]
            ff_type = ob_atom.OBAtom.GetData("FFAtomType").GetValue()
            atom_block.append("    {0:4} {1:5} {2[0]:10.6f} {2[1]:10.6f} {2[2]:10.6f}\n".format(element, ff_type, coord))

            # Generate fingerprint data
            distance = ob_atom.OBAtom.GetDistance(base_atom)
            if distance > max_distance:
                continue
            # Put in distance bin
            fingerprint[int(distance/bin_width)] += 1
            # Put in electronegativity bin
            electronegativity = ob.etab.GetElectroNeg(atomicnum)
            fingerprint[nbins + int(distance/bin_width)] += electronegativity
            # Put in vdw radii
            vdw_radius = ob.etab.GetVdwRad(atomicnum)
            fingerprint[2*nbins + int(distance/bin_width)] += vdw_radius

    fingerprint = ",".join("{:.2f}".format(i) for i in fingerprint)

    #
    # 3D fingerprint
    #

    xmin, xmax = -5.658385, 6.758497
    ymin, ymax = -2.506779, 7.580274
    zmin, zmax = -2.469688, 4.024162
    spacing = 1.0

    # make gridpoints have the cartesian coordinates of all the
    # points of interest on the grid
    x_range = np.arange(xmin-2.0*spacing, xmax+3.0*spacing, spacing)
    y_range = np.arange(ymin-2.0*spacing, ymax+3.0*spacing, spacing)
    z_range = np.arange(zmin-2.0*spacing, zmax+3.0*spacing, spacing)
    gridpoints = [(x, y, z) for x in x_range for y in y_range for z in z_range]
    grid_shape = (len(x_range), len(y_range), len(z_range))

    # Calculate all the atom-point distances
    distance_matrix = cdist(rotated_coordinates, gridpoints)

    # Find charges for all the atoms manually.
    # Automatically would do gasteiger, but fails for 
    # some elements and we use qeq anyway
    qeq = ob.OBChargeModel_FindType('qeq')
    qeq.ComputeCharges(pybel_mol.OBMol)

    # coulomb = q1q2/4pie0r no units yet...
    coulomb_matrix = np.zeros(len(gridpoints))
    for ob_atom, distances in zip(pybel_mol, distance_matrix):
        coulomb_matrix += ob_atom.partialcharge/distances

    # LJ potential based off UFF also no units yet...
    vdw_matrix = np.zeros(len(gridpoints))
    for ob_atom, distances in zip(pybel_mol, distance_matrix):
        # Lorentz-Berthelot mixing rules
        probe = (3.4309, 0.1050)  # Carbon
        source = UFF[ATOMIC_NUMBER[ob_atom.atomicnum]]
        sigma = (source[0] + probe[0]) / 2.0
        epsilon = (source[1] * probe[1])**0.5
        vdw_matrix += 4*epsilon*((sigma/distances)**12 - (sigma/distances)**6)

    # Make into 3D gridded data
    coulomb_matrix = np.reshape(coulomb_matrix, grid_shape)
    vdw_matrix = np.reshape(vdw_matrix, grid_shape)

    # Can clip the maximums here or elsewhere
    coulomb_matrix = np.clip(coulomb_matrix, -0.1, 0.1)
    vdw_matrix = np.clip(vdw_matrix, -10, 0)

    # 3D plotting for visualisation
    #from mayavi import mlab
    #s = mlab.contour3d(coulomb_matrix)
    #s = mlab.contour3d(vdw_matrix)
    #mlab.show()

    #
    # Output
    #

    output_text.append('atoms =\n')
    output_text.extend(atom_block)
    output_text.append('orientation = 0.0 1.0 0.0\n')
    output_text.append('normal = 0.0 0.0 1.0\n')
    output_text.append('carbon_bond = {:.3f}\n'.format(bonds[(-1, 0)][0]))
    output_text.append('fingerprint = {}\n'.format(fingerprint))

    bonds_block = []
    # no bonds < idx 11
    for bond in sorted(bonds):
        if not bond[0] < 0 and not bond[1] < 0:
            bonds_block.append("    {0[0]:4} {0[1]:4} {1[1]:5.2f}\n".format(bond, bonds[bond]))

    output_text.append('bonds =\n')
    output_text.extend(bonds_block[:])

    # Make some pictures; do this now so the ascii can go in the file
    # But first get rid of the benzene
    for _idx in range(10):
        pybel_mol.OBMol.DeleteAtom(pybel_mol.atoms[0].OBAtom)
    pybel_mol.atoms[0].OBAtom.SetType('R')

    if not 'ascii' in pybel.outformats:
        print("Ascii art not available, please upgrade openbabel")
    else:
        ascii_mol = pybel_mol.write(format='ascii', opt={'a': 2, 'w': 40})
        ascii_mol = ['# {}\n'.format(x) for x in ascii_mol.splitlines() if x.strip()]
        output_text[2:2] = ['#\n'] + ascii_mol + ['#\n']

    basename = args.short_name

    pybel_mol.write(format='mol', filename='{}.mol'.format(basename))

    # Always output to a library
    with open('{}.flib'.format(basename), 'w') as out_lib:
        out_lib.writelines(output_text)

    # Make the image with R groups and implicit hydrogen
    unopt_mol = pybel.readstring('smi', "[*:1]" + fgroup)
    unopt_mol.write(format='svg', filename='{}.svg'.format(basename), opt={'C': None})

    # Make a table row in html
    with open('{}.html'.format(basename), 'w') as out_html:
        out_html.write("""\
                <td>{args.short_name}</td>
                <td><p>name: {args.name}</p>
                    <p>smiles: {args.smi_string}</p>
                    <p>MEPO-QEq compatible: {args.mepo_compatible}</td>
                <td><a href="img/{args.short_name}.svg">
                    <img src="img/{args.short_name}.svg"
                         alt="Group: {args.short_name}"
                         title="[{args.short_name}]
                                {args.name}
                                {args.smi_string})"
                         style="height: 75px"/></a>
                </td>
""".format(args=args))

    if args.terminal:
        print("".join(output_text))


if __name__ == '__main__':
    main()

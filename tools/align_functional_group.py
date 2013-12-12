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
    parser.add_argument('-m', '--mepo-compliant', action='store_true',
                        help='Record group as compliant with MEPO-QEq')
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
        "mepo_compliant = {}\n".format(args.mepo_compliant)]

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
                    <p>smiles: {args.smi_string}</p></td>
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

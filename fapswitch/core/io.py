"""

Conversion to and from file formats, mainly cifs of
structures and pickled structures.


"""

try:
    import cPickle as pickle
except ImportError:
    import pickle
import re
import time
from collections import namedtuple
from os import path

from fapswitch.config import options
from fapswitch.core.components import Structure
from fapswitch.core.elements import CCDC_BOND_ORDERS, OB_BOND_ORDERS
from fapswitch.config import info, debug, warning, error
from fapswitch.extensions import sa_score

DOT_FAPSWITCH_VERSION = (7, 0)

Ligand = namedtuple('Ligand', ['smiles', 'inchi', 'inchikey', 'sa_score'])

def load_structure(name):
    """
    Load a structure from a pickle or generate a new one as required.
    Returns an initialised structure. Caches loaded structure on disk.

    """

    pickle_file = "__{}.fapswitch".format(name)
    loaded = False
    if path.exists(pickle_file):
        info("Existing structure found: {}; loading...".format(pickle_file))
        #TODO(tdaff): deal with errors
        try:
            with open(pickle_file, 'rb') as p_structure:
                structure = pickle.load(p_structure)
            # Negative versions ensure that very old caches will be removed
            if not hasattr(structure, 'fapswitch_version'):
                structure.fapswitch_version = (-1, -1)
            # Need to make sure it is still valid
            if structure.fapswitch_version[0] < DOT_FAPSWITCH_VERSION[0]:
                error("Old dot-fapswitch detected, re-initialising")
                loaded = False
            elif structure.fapswitch_version[1] < DOT_FAPSWITCH_VERSION[1]:
                warning("Cached file {} may be out of date".format(pickle_file))
                loaded = True
            else:
                debug("Finished loading")
                loaded = True
        except EOFError:
            warning("Corrupt pickle; re-initialising")
            loaded = False

    if not loaded:
        info("Initialising a new structure. This may take some time.")
        structure = Structure(name)
        structure.from_cif('{}.cif'.format(name))

        # Ensure that atoms in the structure are properly typed
        structure.gen_factional_positions()
        bonding_src = options.get('connectivity')
        if bonding_src == 'file':
            # Rudimentary checks for poor structures
            if not hasattr(structure, 'bonds'):
                error("No bonding in input structure, will probably fail")
            elif len(structure.bonds) == 0:
                error("Zero bonds found, will fail")
            elif not hasattr(structure.atoms[0], 'uff_type'):
                warning("Atoms not properly typed, expect errors")
            else:
                info("Bonding from input file used")
        elif bonding_src == 'openbabel':
            info("Generating topology with Open Babel")
            structure.gen_babel_uff_properties()

        # A couple of structure checks
        structure.check_close_contacts()
        structure.bond_length_check()

        # Initialise the sites after bonds are perceived
        structure.gen_attachment_sites()
        structure.gen_normals()

        # Cache the results
        info("Dumping cache of structure to {}".format(pickle_file))
        debug("dot-fapswitch version {}.{}".format(*DOT_FAPSWITCH_VERSION))
        structure.fapswitch_version = DOT_FAPSWITCH_VERSION
        with open(pickle_file, 'wb') as p_structure:
            pickle.dump(structure, p_structure, protocol=-1)

    return structure


def atoms_to_cif(atoms, cell, bonds, name, identifiers=None):
    """Return a CIF file with bonding and atom types."""

    inv_cell = cell.inverse

    type_count = {}

    atom_part = []
    for atom in atoms:
        if atom is None:
            # blanks are left in here
            continue
        if hasattr(atom, 'uff_type') and atom.uff_type is not None:
            uff_type = atom.uff_type
        else:
            uff_type = '?'
        if atom.element in type_count:
            type_count[atom.element] += 1
        else:
            type_count[atom.element] = 1
        atom.site = "%s%i" % (atom.element, type_count[atom.element])
        atom_part.append("%-5s %-5s %-5s " % (atom.site, atom.element,
                                              uff_type))
        atom_part.append("%f %f %f " % tuple(atom.ifpos(inv_cell)))
        atom_part.append("%f\n" % atom.charge)

    bond_part = []
    for bond, bond_info in bonds.items():
        try:
            bond_length = bond_info[0]
            bond_order = CCDC_BOND_ORDERS[bond_info[1]]
            bond_part.append("%-5s %-5s %f %-5s\n" %
                             (atoms[bond[0]].site, atoms[bond[1]].site,
                              bond_length, bond_order))
        except AttributeError:
            # one of the atoms is None so skip
            #debug("cif NoneType atom")
            pass

    identifiers_part = []
    if identifiers is not None:
        identifiers_part.extend(["loop_\n",
                                 "_chemical_identifier_smiles\n",
                                 "_chemical_identifier_inchi\n",
                                 "_chemical_identifier_inchi_key\n",
                                 "_chemical_synthetic_accessibility\n"])
        for identifiers in identifiers:
            identifiers_part.append("{}  {}  {}  {}\n".format(*identifiers))


    cif_file = [
        "data_%s\n" % name.replace(' ', '_'),
        "%-33s %s\n" % ("_audit_creation_date",
                        time.strftime('%Y-%m-%dT%H:%M:%S%z')),
        "%-33s %s\n" % ("_audit_creation_method",
                        "fapswitch_%i.%i" % DOT_FAPSWITCH_VERSION),
        "%-33s %s\n" % ("_symmetry_space_group_name_H-M", "P1"),
        "%-33s %s\n" % ("_symmetry_Int_Tables_number", "1"),
        "%-33s %s\n" % ("_space_group_crystal_system", cell.crystal_system),
        "%-33s %-.10s\n" % ("_cell_length_a", cell.a),
        "%-33s %-.10s\n" % ("_cell_length_b", cell.b),
        "%-33s %-.10s\n" % ("_cell_length_c", cell.c),
        "%-33s %-.10s\n" % ("_cell_angle_alpha", cell.alpha),
        "%-33s %-.10s\n" % ("_cell_angle_beta", cell.beta),
        "%-33s %-.10s\n" % ("_cell_angle_gamma", cell.gamma),
        "%-33s %s\n" % ("_cell_volume", cell.volume),
        # start of atom loops
        "\nloop_\n",
        "_atom_site_label\n",
        "_atom_site_type_symbol\n",
        "_atom_site_description\n",
        "_atom_site_fract_x\n",
        "_atom_site_fract_y\n",
        "_atom_site_fract_z\n",
        "_atom_type_partial_charge\n"] + atom_part + [
        # bonding loop
        "\nloop_\n",
        "_geom_bond_atom_site_label_1\n",
        "_geom_bond_atom_site_label_2\n",
        "_geom_bond_distance\n",
        "_ccdc_geom_bond_type\n"] + bond_part + ["\n"] + identifiers_part

    return cif_file


def atoms_to_identifiers(atoms, bonds):
    """Derive the smiles for all the organic ligands."""
    try:
        import openbabel as ob
        import pybel
    except ImportError:
        # Don't bother if no openbabel'
        return

    obmol = ob.OBMol()
    obmol.BeginModify()

    # Translation table for indexes
    seen_atoms = {}

    babel_idx = 1

    for idx, atom in enumerate(atoms):
        if atom is None or atom.is_metal: # or atom.atomic_number == 1:
            # If we ignore them it should split the
            # ligands into fragments
            continue
        else:
            new_atom = obmol.NewAtom()
            new_atom.SetAtomicNum(atom.atomic_number)
            # so we correlate the bond index
            # to the index for the babel_mol
            seen_atoms[idx] = babel_idx
            babel_idx += 1

    for bond, bond_info in bonds.items():
        if bond[0] in seen_atoms and bond[1] in seen_atoms:

            obmol.AddBond(seen_atoms[bond[0]],
                          seen_atoms[bond[1]],
                          OB_BOND_ORDERS[bond_info[1]])

    obmol.EndModify()

    pybelmol = pybel.Molecule(obmol)

    # Strip out stereochemistry
    full_molecule = pybelmol.write('can', opt={'i': None}).strip()
    # Fix for delocalised carboxylate detached from metals
    full_molecule = re.sub(r'C\(O\)O([)$.])', r'C(=O)O\1', full_molecule)

    # remove any lone atoms
    unique_smiles = (set(full_molecule.split(".")) -
                     {'O', 'H', 'N'})

    identifiers = []
    for smile in unique_smiles:
        pybelmol = pybel.readstring('smi', smile)
        can_smiles = pybelmol.write('can', opt={'i': None}).strip()
        smol = Ligand(can_smiles,
                      pybelmol.write('inchi', opt={'w': None}).strip(),
                      pybelmol.write('inchikey', opt={'w': None}).strip(),
                      sa_score(can_smiles))
        identifiers.append(smol)

    return identifiers

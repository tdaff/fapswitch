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
from os import path

import numpy as np
from numpy import array, identity, asarray, dot, cross, outer, sin, cos
from numpy import roll
from numpy.linalg import norm

from fapswitch.functional_groups import functional_groups
from fapswitch.core.components import Structure, Atom
from fapswitch.core.components import vecdist3, subgroup
from fapswitch.core.util import min_vect
from fapswitch.core.io import atoms_to_cif
from fapswitch.config import options
from fapswitch.config import debug, info, warning, error, critical
from fapswitch.core.elements import CCDC_BOND_ORDERS


DOT_FAPSWITCH_VERSION = (6, 0)





def rotation_about_angle(axis_in, angle):
    """
    Create a rotation matrix corresponding to the rotation around a general
    axis by a specified angle.

    """
    axis = axis_in/norm(axis_in)

    aprod = outer(axis, axis)
    skew = array([[0, axis[2], -axis[1]], [-axis[2], 0, axis[0]],
        [axis[1], -axis[0], 0]])

    # R = dd^T + cos(a) (I - dd^T) + sin(a) skew(d)
    return aprod+cos(angle)*(identity(3)-aprod)+sin(angle)*skew



def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))




def make_collision_tester(options=None, test_method=None, test_scale=None):
    """
    Create a function that will test atom overlap based on separation
    and atomic radii.

    """

    if test_method is None:
        test_method = options.get('fapswitch_collision_method').lower()
    if test_scale is None:
        test_scale = options.getfloat('fapswitch_collision_scale')

    if test_method == 'covalent':
        info('Covalent radii collision test, scale factor: %f' % test_scale)
        def collision(test_atom, atoms, cell, ignore=()):
            """Covalent radii collision test."""
            # TODO(tdaff): Can't use cleaner .fractional tests as the structure
            # is not passed around to make _parent; can this be fixed
            pos = test_atom.ipos(cell.cell, cell.inverse)
            ipos = test_atom.ifpos(cell.inverse)
            for idx, atom in enumerate(atoms):
                # Skip non atoms
                if atom is None or idx in ignore:
                    continue
                dist = min_vect(pos, ipos, atom.ipos(cell.cell, cell.inverse),
                                atom.ifpos(cell.inverse), cell.cell)
                dist = (dot(dist, dist))**0.5
                min_dist = test_scale*(test_atom.covalent_radius + atom.covalent_radius)
                if dist < min_dist:
                    return False
            return True

    elif test_method == 'vdw':
        info('VdW radii collision test, scale factor: %f' % test_scale)
        def collision(test_atom, atoms, cell, ignore=()):
            """Covalent radii collision test."""
            pos = test_atom.ipos(cell.cell, cell.inverse)
            ipos = test_atom.ifpos(cell.inverse)
            for idx, atom in enumerate(atoms):
                # Skip non atoms
                if atom is None or idx in ignore:
                    continue
                dist = min_vect(pos, ipos, atom.ipos(cell.cell, cell.inverse),
                                atom.ifpos(cell.inverse), cell.cell)
                dist = (dot(dist, dist))**0.5
                min_dist = test_scale*(test_atom.vdw_radius + atom.vdw_radius)
                if dist < min_dist:
                    return False
            return True

    elif test_method == 'cvdw':
        from ccollision import wdist
        info('CVdW radii collision test, scale factor: %f' % test_scale)
        def collision(test_atom, atoms, cell, ignore=()):
            """Covalent radii collision test."""
            pos = test_atom.ipos(cell.cell, cell.inverse)
            ipos = test_atom.ifpos(cell.inverse)
            for idx, atom in enumerate(atoms):
                # Skip non atoms
                if atom is None or idx in ignore:
                    continue
                dist = wdist(pos, ipos, atom.ipos(cell.cell, cell.inverse),
                             atom.ifpos(cell.inverse), cell.cell)
                min_dist = test_scale*(test_atom.vdw_radius + atom.vdw_radius)
                if dist < min_dist:
                    return False
            return True
    else:
        info('Collison test with absolute distance: %f' % test_scale)
        def collision(test_atom, atoms, cell, ignore=()):
            """Covalent radii collision test."""
            pos = test_atom.ipos(cell.cell, cell.inverse)
            ipos = test_atom.ifpos(cell.inverse)
            for idx, atom in enumerate(atoms):
                # Skip non atoms
                if atom is None or idx in ignore:
                    continue
                dist = min_vect(pos, ipos, atom.ipos(cell.cell, cell.inverse),
                                atom.ifpos(cell.inverse), cell.cell)
                dist = (dot(dist, dist))**0.5
                if dist < test_scale:
                    return False
            return True

    # Make a closure for the tester function
    return collision

def arbitrary_normal(vector):
    """Create a normalised normal to an input, does not use random values."""
    if vector == [0, 0, 0]:
        # everything is norma to a zero length vector
        return asarray([1, 0 , 0])
    elif vector[0] == vector[1] == vector[2]:
        # just need cross product with something with non-equal elements
        return normalise(cross(vector, [0.0, 1.0, 0.0]))
    else:
        #different elements so we can swap some to create a different vector
        return normalise(cross(vector, roll(vector, 1)))


def all_combinations_replace(structure, groups, rotations=12, replace_only=None, groups_only=None, max_different=None, backends=()):
    """
    Replace every functional point with every combination of functional groups.

    """

    if replace_only is not None:
        local_attachments = [att_id for att_id in structure.attachments if att_id in replace_only]
        debug("Replacing only: %s" % list(local_attachments))
    else:
        local_attachments = structure.attachments
        debug("Replacing all sites: %s" % list(local_attachments))
    sites = powerset(sorted(local_attachments))

    if groups_only is not None:
        local_groups = [x for x in groups if x in groups_only]
        debug("Using only: %s" % local_groups)
    else:
        local_groups = list(groups)
        debug("Using all groups: %s" % local_groups)

    if max_different is None or max_different <= 0:
        max_different = len(local_groups)

    for site_set in sites:
        for group_set in product(local_groups, repeat=len(site_set)):
            #TODO(tdaff): make this more efficient
            if len(set(group_set)) > max_different:
                continue
            replace_list = zip(group_set, site_set)
            site_replace(structure, groups, replace_list, backends=backends)


def random_combination_replace(structure, groups, rotations=12, replace_only=None, groups_only=None, max_different=0, prob_unfunc=-1.0, backends=()):
    """
    Make a random structure in the site symmetry constrained sample space.

    """

    if replace_only is not None:
        local_attachments = [att_id for att_id in structure.attachments if att_id in replace_only]
        debug("Replacing only: %s" % list(local_attachments))
    else:
        local_attachments = structure.attachments
        debug("Replacing all sites: %s" % list(local_attachments))

    if groups_only is not None:
        local_groups = [x for x in groups if x in groups_only]
        debug("Using only: %s" % local_groups)
    else:
        local_groups = list(groups)
        debug("Using all groups: %s" % local_groups)

    # Limit mixing chemistry with many groups
    if len(local_groups) > max_different > 0:
        local_groups = random.sample(local_groups, max_different)
        debug("Restricted to: %s" % local_groups)

    replace_list = []

    if prob_unfunc < 0:
        # Negative probability means try a random proportion of
        # functionalisation
        prob_unfunc = random.random()
        debug("Random proportion of functionalisation (%f)" % prob_unfunc)

    for site in sorted(local_attachments):
        if random.random() < prob_unfunc:
            # no functional group here
            continue
        else:
            replace_list.append((random.choice(local_groups), site))
    # Do the replacement
    return site_replace(structure, groups, replace_list, backends=backends)


def site_replace(structure, groups, replace_list, rotations=12, backends=()):
    """
    Use replace list to modify the structure with (group, site) pairs.

    Will dump a cif to the backends on success, and return 1 for failed attempt.

    """

    rotation_angle = 2*np.pi/rotations
    new_mof_name = []
    new_mof_friendly_name = []
    # copy the atoms and bonds so we don't alter the original structure
    new_mof = list(structure.atoms)
    new_mof_bonds = dict(structure.bonds)
    for this_group, this_site in replace_list:
        attachment = groups[this_group]
        new_mof_name.append("%s@%s" % (this_group, this_site))
        new_mof_friendly_name.append("%s@%s" % (attachment.name, this_site))
        for this_point in structure.attachments[this_site]:
            attach_id = this_point[0]
            attach_to = this_point[1]
            attach_at = structure.atoms[attach_to].pos
            attach_towards = this_point[2]
            attach_normal = structure.atoms[attach_to].normal
            new_mof[attach_id:attach_id+1] = [None]
            start_idx = len(new_mof)
            for _trial_rotation in range(rotations):
                incoming_group, incoming_bonds = attachment.atoms_attached_to(attach_at, attach_towards, attach_normal, attach_to, start_idx)
                for atom in incoming_group:
                    if not test_collision(atom, new_mof, structure.cell, ignore=[attach_to]):
                        debug("Rotating group")
                        attach_normal = dot(rotation_about_angle(attach_towards, rotation_angle), attach_normal)
                        # Don't need to test more atoms
                        break
                else:
                    # Fits, so add and move on
                    new_mof.extend(incoming_group)
                    new_mof_bonds.update(incoming_bonds)
                    break
            else:
                # Did not attach
                fail_name = ".".join(["@".join(x) for x in replace_list])
                error("Failed: %s@%s from %s" % (this_group, this_site, fail_name))
                return False

    new_mof_name = ".".join(new_mof_name)
    new_mof_friendly_name = ".".join(new_mof_friendly_name)
    info("Generated (%i): [%s]" % (count(), new_mof_friendly_name))

    full_mof_name = "%s_func_%s" % (structure.name, new_mof_name)
    cif_file = atoms_to_cif(new_mof, structure.cell, new_mof_bonds, full_mof_name)

    for backend in backends:
        backend.add_symmetry_structure(structure.name, replace_list, cif_file)

    # successful
    return True


def freeform_replace(structure, groups, replace_only=None, groups_only=None, num_groups=None, custom=None, rotations=36, max_different=0, prob_unfunc=0.5, backends=()):
    """
    Replace sites with no symmetry constraint and with random rotations
    for successive insertion trials (i.e. there will be variation for the
    same structure)

    """
    # Assume that the replace only is passed as a list or iterable
    # default to everything
    if replace_only is None:
        replace_only = list(structure.attachments)
    # Valid list is True where allowed to attach
    valid_list = []
    for attachment, points in structure.attachments.items():
        if attachment in replace_only:
            valid_list.extend([True]*len(points))
        else:
            valid_list.extend([False]*len(points))

    debug("Attachment mask %s" % valid_list)

    nsites = sum(valid_list)

    if groups_only is not None:
        local_groups = [x for x in groups if x in groups_only]
        debug("Using only: %s" % local_groups)
    else:
        local_groups = list(groups)
        debug("Using all groups: %s" % local_groups)

    if len(local_groups) > max_different > 0:
        local_groups = random.sample(local_groups, max_different)
        debug("Restricted to: %s" % local_groups)

    if custom is not None:
        # Specific functionalisation requested
        debug("Processing custom string: %s" % custom)
        func_repr = custom.strip('{}').split(".")
        if len(func_repr) > nsites:
            error("Expected %s sites; got %s" % (nsites, len(func_repr)))
            func_repr = func_repr[:nsites]
            warning("Truncated to {%s}" % ".".join(func_repr))
        elif len(func_repr) < nsites:
            error("Expected %s sites; got %s" % (nsites, len(func_repr)))
            func_repr = func_repr + ['']*(nsites - len(func_repr))
            warning("Padded to {%s}" % ".".join(func_repr))
        for unmasked, site in zip(valid_list, func_repr):
            if unmasked is False and site != '':
                warning("Replacing masked site")
    else:
        # Randomise the selection
        if num_groups is None:
            num_groups = random.randint(1, nsites)
            debug("Randomly replacing %i sites" % num_groups)
        elif num_groups > nsites:
            warning("Too many sites requested; changing all %i" % nsites)
            num_groups = nsites
        #TODO(tdaff): selected groups only
        func_repr = [random.choice(local_groups) for _counter in range(num_groups)]
        # Pad to the correct length
        func_repr.extend([""]*(nsites - num_groups))
        # Randomise
        random.shuffle(func_repr)
        # These need to be put in the unmasked slots
        masked_func_repr = []
        for unmasked in valid_list:
            if unmasked:
                masked_func_repr.append(func_repr.pop(0))
            else:
                masked_func_repr.append('')
        func_repr = masked_func_repr
    # Unique-ish
    unique_name = hashlib.md5(str(func_repr).encode('utf-8')).hexdigest()
    new_mof_name = []
    new_mof = list(structure.atoms)
    new_mof_bonds = dict(structure.bonds)
    for this_point, this_group in zip(chain(*[structure.attachments[x] for x in sorted(structure.attachments)]), func_repr):
        if this_group == "":
            new_mof_name.append("")
            continue
        else:
            new_mof_name.append(this_group)
        attachment = groups[this_group]
        attach_id = this_point[0]
        attach_to = this_point[1]
        attach_at = structure.atoms[attach_to].pos
        attach_towards = this_point[2]
        attach_normal = structure.atoms[attach_to].normal
        #extracted_atoms = new_mof[attach_id:attach_id+1]
        new_mof[attach_id:attach_id+1] = [None]
        start_idx = len(new_mof)
        for trial_rotation in range(rotations):
            incoming_group, incoming_bonds = attachment.atoms_attached_to(attach_at, attach_towards, attach_normal, attach_to, start_idx)
            for atom in incoming_group:
                if not test_collision(atom, new_mof, structure.cell, ignore=[attach_to]):
                    debug("Randomly rotating group")
                    attach_normal = dot(rotation_about_angle(attach_towards, random.random()*np.pi*2), attach_normal)
                    break
            else:
                # Fits, so add and move on
                new_mof.extend(incoming_group)
                new_mof_bonds.update(incoming_bonds)
                break
        else:
            # this_point not valid
            error("Failed to generate: %s" % ".".join([x or "" for x in func_repr]))
            warning("Stopped after: %s" % ".".join(new_mof_name))
            return False

    new_mof_name = ".".join(new_mof_name)
    info("Generated (%i): {%s}" %  (count(), new_mof_name))
    info("With unique name: %s" % unique_name)

    full_mof_name = "%s_free_%s" % (structure.name, new_mof_name)
    cif_file = atoms_to_cif(new_mof, structure.cell, new_mof_bonds, full_mof_name)

    for backend in backends:
        backend.add_freeform_structure(structure.name, func_repr, cif_file)

    # completed sucessfully
    return True


def label_atom(element=None, site=None):
    """Produce unique atom labels."""
    if not hasattr(label_atom, 'seen'):
        label_atom.seen = set()
    if not hasattr(label_atom, 'index'):
        label_atom.index = 1
    if element is not None:
        label = "%s%i" % (element, label_atom.index)
        while label in label_atom.seen:
            label_atom.index += 1
            label = "%s%i" % (element, label_atom.index)
        label_atom.seen.add(label)
        return label
    elif site is not None:
        label_atom.seen.add(site)


def count(reset=False):
    """Return the next itneger from a global state."""
    if not hasattr(count, 'idx') or reset is True:
        count.idx = 1
    elif reset:
        count.idx = reset
    else:
        count.idx += 1
    return count.idx


def fapswitch_deamon(options, structure, f_groups, backends):
    """
    Use sockets to listen and receive structures.

    """
    timeout = 7200
    # set this to zero for random available port
    port = options.getint("fapswitch_port")
    listener = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    # no host as this is running locally
    listener.bind(('', port))
    port = listener.getsockname()[1]
    if options.getint('verbosity') < 0:
        critical("Listening on port %i ..." % port)
    else:
        info("Listening on port %i ..." % port)
    listener.settimeout(timeout)
    listener.listen(1)
    # We only wait for a little bit so that process doesn't zombie forever
    try:
        conn, addr = listener.accept()
    except socket.timeout:
        error("No connection within %i seconds; exiting" % timeout)
        listener.close()
        return False
    info('Connected by %s' % str(addr))

    # Server will continue until an empty input is sent or something times out
    conn.settimeout(timeout)
    while 1:
        try:
            line = conn.recv(1024)
        except socket.timeout:
            error("Timed out after %i sconds waiting for input" % timeout)
            return False

        # Empty input closes server
        if not line:
            break

        # collect information to send what has been done back
        processed = []

        # freeform strings are in braces {}, no spaces
        free_strings = re.findall('\{(.*?)\}', line)
        debug("Freeform strings: %s" % str(free_strings))
        for free_string in free_strings:
            complete = freeform_replace(structure, f_groups, custom=free_string, backends=backends)
            processed.append('{%s}' % free_string)
            processed.append('%s' % complete)

        # site replacements in square brackets [], no spaces
        site_strings = re.findall('\[(.*?)\]', line)
        debug("Site replacement strings: %s" % str(site_strings))
        for site_string in site_strings:
            site_list = [x.split('@') for x in site_string.split('.') if x]
            debug(str(site_list))
            complete = site_replace(structure, f_groups, site_list, backends=backends)
            processed.append('[%s]' % site_string)
            processed.append('%s' % complete)

        try:
            conn.sendall(':'.join(processed))
        except conn.timeout:
            error("Timed out sending status after %i seconds" % timeout)
            return False

    conn.close()
    return True


def main():
    """
    Run the substitution for an input structure.

    """

    job_options = options
    job_name = job_options.get('job_name')

    # Load an existing pickled structure or generate a new one
    pickle_file = "__%s.fapswitch" % job_name
    loaded = False
    if path.exists(pickle_file):
        info("Existing structure found: %s; loading..." % pickle_file)
        with open(pickle_file, 'rb') as load_structure:
            input_structure = pickle.load(load_structure)
        # Negative versions ensure that very old caches will be removed
        if not hasattr(input_structure, 'fapswitch_version'):
            input_structure.fapswitch_version = (-1, -1)
        # Need to make sure it is still valid
        if input_structure.fapswitch_version[0] < DOT_FAPSWITCH_VERSION[0]:
            error("Old dot-fapswitch detected, re-initialising")
            loaded = False
        elif input_structure.fapswitch_version[1] < DOT_FAPSWITCH_VERSION[1]:
            warning("Cached file %s may be out of date" % pickle_file)
            loaded = True
        else:
            debug("Finished loading")
            loaded = True

    if not loaded:
        info("Initialising a new structure. This may take some time.")
        input_structure = Structure(job_name)
        input_structure.from_cif('{}.cif'.format(job_name))

        # Ensure that atoms in the structure are properly typed
        input_structure.gen_factional_positions()
        bonding_src = job_options.get('fapswitch_connectivity')
        if bonding_src == 'file':
            # Rudimentary checks for poor structures
            if not hasattr(input_structure, 'bonds'):
                error("No bonding in input structure, will probably fail")
            elif len(input_structure.bonds) == 0:
                error("Zero bonds found, will fail")
            elif not hasattr(input_structure.atoms[0], 'uff_type'):
                warning("Atoms not properly typed, expect errors")
            else:
                info("Bonding from input file used")
        elif bonding_src == 'openbabel':
            info("Generating topology with Open Babel")
            input_structure.gen_babel_uff_properties()

        # A couple of structure checks
        input_structure.check_close_contacts()
        input_structure.bond_length_check()

        # Initialise the sites after bonds are perceived
        input_structure.gen_attachment_sites()
        input_structure.gen_normals()

        # Cache the results
        info("Dumping cache of structure connectivity to %s" % pickle_file)
        debug("dot-fapswitch version %i.%i" % DOT_FAPSWITCH_VERSION)
        input_structure.fapswitch_version = DOT_FAPSWITCH_VERSION
        with open(pickle_file, 'wb') as save_structure:
            pickle.dump(input_structure, save_structure)

    # Structure is ready!
    # Begin processing
    info("Structure attachment sites: %s" % list(input_structure.attachments))

    # Will use selected sites if specified, otherwise use all
    replace_only = job_options.gettuple('fapswitch_replace_only')
    if replace_only == ():
        replace_only = None

    # label_atom has a global state that
    for atom in input_structure.atoms:
        label_atom(site=atom.site)

    # Functional group library is self initialising
    f_groups = functional_groups
    info("Groups in library: %s" % str(f_groups.group_list))

    global test_collision
    test_collision = make_collision_tester(job_options)

    #Define some backends for where to send the structures
    backends = []
    backend_options = job_options.gettuple('fapswitch_backends')

    if 'sqlite' in backend_options:
        # Initialise and add the database writer
        debug("Initialising the sqlite backend")
        try:
            from fapswitch.backend.sql import AlchemyBackend
            backend = AlchemyBackend(job_name)
            backend.populate_groups(f_groups)
            backends.append(backend)
        except ImportError:
            error("SQLAlchemy not installed; sql backend unavailable")
        # done

    if 'file' in backend_options:
        # Just dumps to a named file
        debug("Initialising cif file writer backend")
        from fapswitch.backend.cif_file import CifFileBackend
        backends.append(CifFileBackend())


    # Decide if we should run the server mode
    if job_options.getbool('daemon'):
        # Make the program die if the daemon is called unsuccessfully
        success = fapswitch_deamon(job_options, input_structure, f_groups, backends=backends)
        if success is False:
            raise SystemExit

    # User defined, single-shot functionalisations
    custom_strings = job_options.get('fapswitch_custom_strings')
    # Pattern matching same as in the daemon
    # freeform strings are in braces {}, no spaces
    freeform_strings = re.findall('\{(.*?)\}', custom_strings)
    debug("Freeform option strings: %s" % str(freeform_strings))
    for freeform_string in freeform_strings:
        freeform_replace(input_structure, f_groups, custom=freeform_string, backends=backends)
    # site replacements in square brackets [], no spaces
    site_strings = re.findall('\[(.*?)\]', custom_strings)
    debug("Site replacement options strings: %s" % str(site_strings))
    for site_string in site_strings:
        # These should be functional_group1@site1.functional_group2@site2
        site_list = [x.split('@') for x in site_string.split('.') if x]
        debug(str(site_list))
        site_replace(input_structure, f_groups, site_list, backends=backends)

    # Full systematic replacement of everything start here

    # Only use these functional groups for replacements
    replace_groups = job_options.gettuple('fapswitch_replace_groups')
    if replace_groups == ():
        replace_groups = None

    max_different = job_options.getint('fapswitch_max_different')

    prob_unfunc = job_options.getfloat('fapswitch_unfunctionalised_probability')

    if job_options.getbool('fapswitch_replace_all_sites'):
        all_combinations_replace(input_structure, f_groups, replace_only=replace_only, groups_only=replace_groups, max_different=max_different, backends=backends)

    # group@site randomisations
    random_count = job_options.getint('fapswitch_site_random_count')
    successful_randoms = 0
    while successful_randoms < random_count:
        #function returns true if structure is generated
        if random_combination_replace(input_structure, f_groups, replace_only=replace_only, groups_only=replace_groups, max_different=max_different, prob_unfunc=prob_unfunc, backends=backends):
            successful_randoms += 1
            info("Generated %i of %i site random structures" % (successful_randoms, random_count))

    # fully freeform randomisations
    random_count = job_options.getint('fapswitch_full_random_count')
    successful_randoms = 0
    while successful_randoms < random_count:
        #function returns true if structure is generated
        if freeform_replace(input_structure, f_groups, replace_only=replace_only, groups_only=replace_groups, max_different=max_different, prob_unfunc=prob_unfunc, backends=backends):
            successful_randoms += 1
            info("Generated %i of %i fully random structures" % (successful_randoms, random_count))


if __name__ == '__main__':

    main()

else:
    test_collision = make_collision_tester(test_method='vdw', test_scale=0.5)

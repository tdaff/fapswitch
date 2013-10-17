"""

Core methods that do functional group replacement.

"""

import hashlib
import random
from itertools import product, chain

import numpy as np
from numpy import dot

from fapswitch.config import debug, info, warning, error, critical
from fapswitch.core.util import powerset
from fapswitch.core.collision import make_collision_tester
from fapswitch.core.io import atoms_to_cif, atoms_to_smiles
from fapswitch.core.util import rotation_about_angle
from fapswitch.functional_groups import functional_groups

test_collision = make_collision_tester(test_method='cvdw', test_scale=1.122462048309373)



def count(reset=False):
    """Return the next itneger from a global state."""
    if not hasattr(count, 'idx') or reset is True:
        count.idx = 1
    elif reset:
        count.idx = reset
    else:
        count.idx += 1
    return count.idx


def all_combinations_replace(structure, rotations=12, replace_only=None, groups_only=None, max_different=None, backends=()):
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
        local_groups = [x for x in functional_groups if x in groups_only]
        debug("Using only: %s" % local_groups)
    else:
        local_groups = list(functional_groups)
        debug("Using all groups: %s" % local_groups)

    if max_different is None or max_different <= 0:
        max_different = len(local_groups)

    for site_set in sites:
        for group_set in product(local_groups, repeat=len(site_set)):
            #TODO(tdaff): make this more efficient
            if len(set(group_set)) > max_different:
                continue
            replace_list = zip(group_set, site_set)
            site_replace(structure, replace_list, backends=backends)


def random_combination_replace(structure, rotations=12, replace_only=None, groups_only=None, max_different=0, prob_unfunc=-1.0, backends=()):
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
        local_groups = [x for x in functional_groups if x in groups_only]
        debug("Using only: %s" % local_groups)
    else:
        local_groups = list(functional_groups)
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
        debug("Random functionalisation proportion ({})".format(prob_unfunc))

    for site in sorted(local_attachments):
        if random.random() < prob_unfunc:
            # no functional group here
            continue
        else:
            replace_list.append((random.choice(local_groups), site))
    # Do the replacement
    return site_replace(structure, replace_list, backends=backends)


def site_replace(structure, replace_list, rotations=12, backends=()):
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
        attachment = functional_groups[this_group]
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
    atoms_to_smiles(new_mof, new_mof_bonds)

    for backend in backends:
        backend.add_symmetry_structure(structure.name, replace_list, cif_file)

    # successful
    return True


def freeform_replace(structure, replace_only=None, groups_only=None, num_groups=None, custom=None, rotations=36, max_different=0, prob_unfunc=0.5, backends=()):
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
        local_groups = [x for x in functional_groups if x in groups_only]
        debug("Using only: %s" % local_groups)
    else:
        local_groups = list(functional_groups)
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
        attachment = functional_groups[this_group]
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

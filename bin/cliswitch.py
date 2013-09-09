#!/usr/bin/env python

"""
fapswitch.py

Alter a known structure with new functional groups ready for fapping.

"""


try:
    import cPickle as pickle
except ImportError:
    import pickle
import re
import socket
from os import path


from fapswitch.functional_groups import functional_groups
from fapswitch.core.components import Structure
from fapswitch.config import options
from fapswitch.config import debug, info, warning, error, critical
from fapswitch.core.methods import site_replace, all_combinations_replace
from fapswitch.core.methods import freeform_replace, random_combination_replace


DOT_FAPSWITCH_VERSION = (6, 0)


def fapswitch_deamon(structure, f_groups, backends):
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
        free_strings = re.findall('{(.*?)}', line)
        debug("Freeform strings: %s" % str(free_strings))
        for free_string in free_strings:
            complete = freeform_replace(structure, f_groups,
                                        custom=free_string, backends=backends)
            processed.append('{%s}' % free_string)
            processed.append('%s' % complete)

        # site replacements in square brackets [], no spaces
        site_strings = re.findall(r'\[(.*?)\]', line)
        debug("Site replacement strings: %s" % str(site_strings))
        for site_string in site_strings:
            site_list = [x.split('@') for x in site_string.split('.') if x]
            debug(str(site_list))
            complete = site_replace(structure, f_groups,
                                    site_list, backends=backends)
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

    # Functional group library is self initialising
    f_groups = functional_groups
    info("Groups in library: %s" % str(f_groups.group_list))

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
        success = fapswitch_deamon(input_structure, f_groups,
                                   backends=backends)
        if success is False:
            raise SystemExit

    # User defined, single-shot functionalisations
    custom_strings = job_options.get('fapswitch_custom_strings')
    # Pattern matching same as in the daemon
    # freeform strings are in braces {}, no spaces
    freeform_strings = re.findall('{(.*?)}', custom_strings)
    debug("Freeform option strings: %s" % str(freeform_strings))
    for freeform_string in freeform_strings:
        freeform_replace(input_structure, f_groups, custom=freeform_string,
                         backends=backends)
    # site replacements in square brackets [], no spaces
    site_strings = re.findall(r'\[(.*?)\]', custom_strings)
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
        all_combinations_replace(input_structure, f_groups,
                                 replace_only=replace_only,
                                 groups_only=replace_groups,
                                 max_different=max_different,
                                 backends=backends)

    # group@site randomisations
    random_count = job_options.getint('fapswitch_site_random_count')
    successful_randoms = 0
    while successful_randoms < random_count:
        #function returns true if structure is generated
        if random_combination_replace(input_structure, f_groups,
                                      replace_only=replace_only,
                                      groups_only=replace_groups,
                                      max_different=max_different,
                                      prob_unfunc=prob_unfunc,
                                      backends=backends):
            successful_randoms += 1
            info("Generated %i of %i site random structures" %
                 (successful_randoms, random_count))

    # fully freeform randomisations
    random_count = job_options.getint('fapswitch_full_random_count')
    successful_randoms = 0
    while successful_randoms < random_count:
        #function returns true if structure is generated
        if freeform_replace(input_structure, f_groups,
                            replace_only=replace_only,
                            groups_only=replace_groups,
                            max_different=max_different,
                            prob_unfunc=prob_unfunc,
                            backends=backends):
            successful_randoms += 1
            info("Generated %i of %i fully random structures" %
                 (successful_randoms, random_count))


if __name__ == '__main__':

    main()

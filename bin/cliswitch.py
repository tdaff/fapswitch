#!/usr/bin/env python

"""
cliswitch.py

Command line interface to the fapswitch library functions
which alter a known CIF structure with new functional groups
ready for fapping.

"""


import re
import sys
from os.path import dirname, realpath
# Put the parent fapswitch first in the path
sys.path.insert(1, dirname(dirname(realpath(__file__))))

from fapswitch.config import options
from fapswitch.config import debug, info, error
from fapswitch.core.io import load_structure
from fapswitch.core.methods import site_replace, all_combinations_replace
from fapswitch.core.methods import freeform_replace, random_combination_replace
from fapswitch.functional_groups import functional_groups


def main():
    """
    Simple application that will read options and run the substitution
    for an input structure.

    """

    # Name for a the single structure
    job_name = options.get('job_name')

    # Load it
    input_structure = load_structure(job_name)
    # Structure is ready!

    # Begin processing
    info("Structure attachment sites: "
         "{}".format(list(input_structure.attachments)))

    # Will use selected sites if specified, otherwise use all
    replace_only = options.gettuple('replace_only')
    if replace_only == ():
        replace_only = None

    # Functional group library is self initialising
    info("Groups in library: {}".format(functional_groups.group_list))

    #Define some backends for where to send the structures
    backends = []
    backend_options = options.gettuple('backends')

    if 'sqlite' in backend_options:
        # Initialise and add the database writer
        debug("Initialising the sqlite backend")
        try:
            from fapswitch.backend.sql import AlchemyBackend
            backend = AlchemyBackend(job_name)
            backend.populate_groups(functional_groups)
            backends.append(backend)
        except ImportError:
            error("SQLAlchemy not installed; sql backend unavailable")
        # done

    if 'file' in backend_options:
        # Just dumps to a named file
        debug("Initialising cif file writer backend")
        from fapswitch.backend.cif_file import CifFileBackend
        backends.append(CifFileBackend())


    ##
    # User defined, single-shot functionalisations
    ##

    custom_strings = options.get('custom_strings')
    # Pattern matching same as in the daemon
    # freeform strings are in braces {}, no spaces
    freeform_strings = re.findall('{(.*?)}', custom_strings)
    debug("Freeform option strings: {}".format(freeform_strings))
    for freeform_string in freeform_strings:
        freeform_replace(input_structure, custom=freeform_string,
                         backends=backends)

    # site replacements in square brackets [], no spaces
    site_strings = re.findall(r'\[(.*?)\]', custom_strings)
    debug("Site replacement options strings: {}".format(site_strings))
    for site_string in site_strings:
        # These should be functional_group1@site1.functional_group2@site2
        site_list = [x.split('@') for x in site_string.split('.') if x]
        debug(str(site_list))
        site_replace(input_structure, site_list, backends=backends)

    ##
    # Full systematic replacement of everything start here
    ##

    # Only use these functional groups for replacements
    replace_groups = options.gettuple('replace_groups')
    if replace_groups == ():
        replace_groups = None

    max_different = options.getint('max_different')

    prob_unfunc = options.getfloat('unfunctionalised_probability')

    # Do absolutely every combination (might take a while)
    if options.getbool('replace_all_sites'):
        all_combinations_replace(input_structure,
                                 replace_only=replace_only,
                                 groups_only=replace_groups,
                                 max_different=max_different,
                                 backends=backends)

    # group@site randomisations
    random_count = options.getint('site_random_count')
    successful_randoms = 0
    while successful_randoms < random_count:
        #function returns true if structure is generated
        if random_combination_replace(input_structure,
                                      replace_only=replace_only,
                                      groups_only=replace_groups,
                                      max_different=max_different,
                                      prob_unfunc=prob_unfunc,
                                      backends=backends):
            successful_randoms += 1
            info("Generated %i of %i site random structures" %
                 (successful_randoms, random_count))

    # fully freeform randomisations
    random_count = options.getint('full_random_count')
    successful_randoms = 0
    while successful_randoms < random_count:
        #function returns true if structure is generated
        if freeform_replace(input_structure,
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

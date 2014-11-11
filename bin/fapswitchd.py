#!/usr/bin/env python

"""
fapswitchd.py

A socket server for functionalisation of MOF structures.
Will start up a server that accepts incoming socket connections
and responds to

"""


import re
import socket
import sys
from os.path import dirname, realpath
sys.path.insert(1, dirname(dirname(realpath(__file__))))

import fapswitch
from fapswitch.config import options
from fapswitch.config import debug, info, error, critical
from fapswitch.core.io import load_structure
from fapswitch.core.methods import site_replace
from fapswitch.core.methods import freeform_replace
from fapswitch.functional_groups import functional_groups


def fapswitch_deamon(structure, backends, rotations=12):
    """
    Use sockets to listen and receive structures.

    """
    timeout = options.getint('timeout')
    # set this to zero for random available port
    port = options.getint('port')
    listener = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    # no host as this is running locally
    listener.bind(('', port))
    port = listener.getsockname()[1]

    # Ensure that we always know the port
    critical("Listening on port {} ...".format(port))

    listener.settimeout(timeout)
    listener.listen(1)
    # We only wait for a little bit so that process doesn't zombie forever
    try:
        conn, addr = listener.accept()
    except socket.timeout:
        error("No connection within {} seconds; exiting".format(timeout))
        listener.close()
        return False
    info('Connected by {}'.format(addr))

    # Server will continue until an empty input is sent or something times out
    conn.settimeout(timeout)
    while 1:
        try:
            line = conn.recv(1024).decode('utf-8')
        except socket.timeout:
            error("Timed out after {} sconds waiting for input".format(timeout))
            return False

        # Empty input closes server
        if not line:
            break

        # collect information to send what has been done back
        processed = []

        # freeform strings are in braces {}, no spaces
        free_strings = re.findall('{(.*?)}', line)
        debug("Freeform strings: {}".format(free_strings))
        for free_string in free_strings:
            complete = freeform_replace(structure, custom=free_string,
                                        backends=backends, rotations=12)
            processed.append('{{{}}}'.format(free_string))
            processed.append('{}'.format(complete))

        # site replacements in square brackets [], no spaces
        site_strings = re.findall(r'\[(.*?)\]', line)
        debug("Site replacement strings: {}".format(site_strings))
        for site_string in site_strings:
            site_list = [x.split('@') for x in site_string.split('.') if x]
            debug("{}".format(site_list))
            complete = site_replace(structure, site_list, backends=backends,
                                    rotations=12)
            processed.append('[{}]'.format(site_string))
            processed.append('{}'.format(complete))

        try:
            conn.sendall((':'.join(processed)).encode('utf-8'))
        except conn.timeout:
            error("Timed out sending status after {} seconds".format(timeout))
            return False

    conn.close()
    return True


def main():
    """
    Initialise everything needed to get the daemon running and start it up.

    """

    info("Welcome to fapswitchd; the daemon interface to fapswitch")
    info("Using fapswitch version {}".format(fapswitch.__version__))

    # Name for a the single structure
    job_name = options.get('job_name')

    # Load it
    input_structure = load_structure(job_name)
    # Structure is ready!

    # Begin processing
    info("Structure attachment sites: "
         "{}".format(list(input_structure.attachments)))
    info("Structure attachment multiplicities :"
         "{}".format(dict((key, len(val)) for key, val in
                          input_structure.attachments.items())))

    # Functional group library is self initialising
    info("Groups in library: {}".format(functional_groups.group_list))

    #Define some backends for where to send the structures
    backends = []
    backend_options = options.gettuple('backends')

    rotations = options.getint('rotations')
    info("Will rotate each group a maximum of {} times.".format(rotations))

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

    # Make the program die if the daemon is called unsuccessfully
    if fapswitch_deamon(input_structure, backends=backends,
                        rotations=rotations):
        info("Daemon completed successfully")
    else:
        error("Daemon did not complete successfully; check output")


if __name__ == '__main__':

    main()

#!/usr/bin/env python

"""
webswitch.py

A web server frontend that answers requests with newly functionalised
MOFs.

"""


import glob
import os
import random
import re
import sys
from os import path
from os.path import dirname, realpath
# Put the parent fapswitch first in the path
sys.path.insert(1, dirname(dirname(realpath(__file__))))

import tornado.web
import tornado.ioloop
import tornado.template

import fapswitch
from fapswitch.bibliography import references
from fapswitch.config import options
from fapswitch.config import debug, info
from fapswitch.core.io import load_structure
from fapswitch.functional_groups import functional_groups
from fapswitch.core.methods import random_combination_replace
from fapswitch.backend.web_store import WebStoreBackend

TORNADO_PORT = 8888

web_dir = path.join(path.dirname(fapswitch.__file__), 'web')
datastore = os.path.join(web_dir, 'datastore')
ligand_dir = os.path.join(web_dir, 'generated')

templates = tornado.template.Loader(path.join(web_dir, 'templates'))


# Send the index file
class IndexHandler(tornado.web.RequestHandler):
    """
    Respond to the root path with an intro page.
    """

    def get(self, url='/'):
        """Render index."""
        self.index()

    def post(self, url='/'):
        """Render index."""
        self.index()

    def index(self):
        """Make the index page and send it out."""
        my_refs = [references['Daff2013']]
        page = templates.load('index.html').generate(references=my_refs)
        self.write(page)


# Random webpage generator
class RandomHandler(tornado.web.RequestHandler):
    """
    Respond to a request for a random structure with a web view.
    """

    #both GET and POST requests have the same responses
    def get(self, url='/'):
        """GET request"""
        debug("GET request")
        self.handle_request()

    def post(self, url='/'):
        """POST request"""
        debug('POST request')
        self.handle_request()

    def handle_request(self):
        """Generate a random structure and return the rendered page."""

        debug("Arguments: {}".format(self.request.arguments))
        groups = self.get_arguments('groups')
        if not groups:
            groups = None
        info("Groups selected: {}".format(groups))

        max_trials = 10
        # Possible options:
        # replace_only: tuple of sites to replace
        # groups_only: only use specific groups
        # max_different: resrict simultaneous types of groups

        # This is new every time and keeps all the information
        # we need specific to the web version
        backends = [WebStoreBackend()]

        for _trial in range(max_trials):

            debug("Trial {}".format(_trial))
            base_structure = random.choice(initialised_structures)
            status = random_combination_replace(groups_only=groups,
                                                structure=base_structure,
                                                backends=backends)
            if status:
                cif_info = backends[0].cifs[0]
                # MEPO compatibility if all groups are okay
                if all(functional_groups[function[0]].mepo_compatible for
                       function in cif_info['functions']):
                    mepo_compatible = "Yes"
                else:
                    mepo_compatible = "No"

                collision_tester = options.get('collision_method')
                collision_cutoff = options.getfloat('collision_scale')

                if cif_info['ligands'] is None:
                    ligands = []
                else:
                    ligands = cif_info['ligands']

                if ligands:
                    sa_score = max(ligand.sa_score for ligand in ligands)
                else:
                    sa_score = 0.0

                processed_ligands = make_ligands(ligands)

                extra_info = """<h4>Hypothetical functionalised MOF</h4>
                <p>Functional groups have been added using the crystal
                symmetry. A collision detection routine with a {} radius
                at {:.2f} was used to carry out the functionalisation.
                Note that although atoms may appear close, the bonding
                connectivity defined in the cif file will be correct.</p>
                """.format(collision_tester, collision_cutoff)

                # These references are always required
                local_references = [
                    references['Kadantsev2013'],
                    references['Ertl2009']]

                # Find all the references and add them too
                for reference in re.findall(r'\[(.*?)\]', extra_info):
                    local_references.append(references[reference])

                # Raw HTML anchors. Ugly.
                extra_info = re.sub(
                    r'\[(.*?)\]',  # non-greedy(?) find in square brackets
                    r'[<a href="#\1">\1</a>]',  # replace raw html
                    extra_info)

                page = templates.load('random.html').generate(
                    mepo_compatible=mepo_compatible,
                    references=local_references,
                    functional_groups=functional_groups,
                    extra_info=extra_info,
                    sa_score=sa_score,
                    processed_ligands=processed_ligands,
                    **cif_info)
                self.write(page)

                break
        else:
            page = templates.load('failed.html').generate()
            self.write(page)


# adds event handlers for commands and file requests
application = tornado.web.Application([
    (r"/(random.*)", RandomHandler),
    (r"/", IndexHandler),
    (r"/(index.*)", IndexHandler),
    (r"/(.*\.cif)", tornado.web.StaticFileHandler,
     {"path": path.join(web_dir, 'generated')}),
    (r"/(.*\.png)", tornado.web.StaticFileHandler,
     {"path": path.join(web_dir, 'generated')}),
    (r"/(.*\.jpg)", tornado.web.StaticFileHandler,
     {"path": path.join(web_dir, 'static')}),
    (r"/(.*\.js)", tornado.web.StaticFileHandler,
     {"path": path.join(web_dir, 'js')}),
    (r"/(.*\.css)", tornado.web.StaticFileHandler,
     {"path": path.join(web_dir, 'css')}),
])


def make_ligands(ligands):
    """
    Process lignds to find the maximum SA Score and make png files.
    Takes a list of Ligand namedtuples and returns a list of:
    (smiles, png_filename, sa_score).
    """

    out_ligands = []

    try:
        import pybel
        import openbabel as ob
        # Test if we have png2 or _png2 as it changes in the
        # 2.3.1 and 2.3.2 versions
        obconversion = ob.OBConversion()
        formatok = obconversion.SetOutFormat("png2")
        if formatok:
            png_fmt = 'png2'
        else:
            png_fmt = '_png2'

    except ImportError:
        return out_ligands

    for ligand in ligands:
        png_file = "{}.png".format(ligand.inchikey)
        png_path = path.join(ligand_dir, png_file)

        # Make the png if it doesn't exist yet
        if not os.path.exists(png_path):
            smi = pybel.readstring('smi', ligand.smiles)
            smi.write(png_fmt, png_path)

        # Add everything to the list
        out_ligands.append((ligand.smiles, png_file, ligand.sa_score))

    return out_ligands


if __name__ == "__main__":

    #
    # Fun fapswitch as a web server
    #
    # Everything lives in the global scope

    # Directory to store pickles

    initial_directory = os.getcwd()

    os.chdir(datastore)

    available_structures = glob.glob('*.cif')

    initialised_structures = []

    for available_structure in available_structures:
        sname = available_structure[:-4]
        initialised_structures.append(load_structure(sname))

    os.chdir(initial_directory)

    # Functional group library is self initialising
    info("Groups in library: %s" % str(functional_groups.group_list))

    #start tornado
    application.listen(TORNADO_PORT)
    info("Starting server on port number {}...".format(TORNADO_PORT))
    info("Open at http://127.0.0.1:{}/index.html".format(TORNADO_PORT))
    tornado.ioloop.IOLoop.instance().start()

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
from collections import OrderedDict
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
from fapswitch.core.methods import random_combination_replace, site_replace
from fapswitch.backend.web_store import WebStoreBackend

TORNADO_PORT = 8888
STRUCTURE_CACHE = 20

web_dir = path.join(path.dirname(fapswitch.__file__), 'web')
datastore = path.join(web_dir, 'datastore')
ligand_dir = path.join(web_dir, 'generated')

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

        max_trials = 20
        top_50_groups = ["Me", "Ph", "Cl", "OMe", "OH", "Et", "OEt", "F", "Br",
                         "NO2", "NH2", "CN", "COOEt", "COMe", "COOH", "Bnz",
                         "COOMe", "iPr", "pTol", "4ClPh", "tBu", "4OMePh",
                         "CF3", "COPh", "Pr", "NMe2", "Bu", "OBnz", "4NO2Ph",
                         "OAc", "4FPh", "I", "4BrPh", "2ClPh", "All", "COH",
                         "SMe", "CONH2", "NPh", "24DClPh", "CHex", "Morph",
                         "HCO", "3ClPh", "oTol", "2Fur", "iBu", "NCOMe"]
        small_groups = ["F", "Cl", "Me", "NH2", "OH", "CN"]
        # Possible options:
        # replace_only: tuple of sites to replace
        # groups_only: only use specific groups
        # max_different: restrict simultaneous types of groups

        # This is new every time and keeps all the information
        # we need specific to the web version
        backends = [WebStoreBackend()]
        failed = ""

        if 'mof-choice' in self.request.arguments:
            # Selected a specific MOF
            chosen_structure = self.get_argument('mof-choice')
            base_structure = get_structure(chosen_structure)
            replace_list = []
            for site in available_structures[chosen_structure]:
                group = self.get_argument(site, None)
                if group is None or 'None' in group:
                    continue
                elif 'Random' in group:
                    replace_list.append([random.choice(top_50_groups), site])
                else:
                    replace_list.append([group, site])

            # Now make the MOF
            status = site_replace(base_structure, replace_list=replace_list,
                                  backends=backends)
            if not status:
                # couldn't make it so just use clean structure
                failed = ".".join("{}@{}".format(x[0], x[1])
                                  for x in replace_list)
                site_replace(base_structure, replace_list=[],
                             backends=backends)
        else:
            # Completely random
            chosen_structure = random.choice(list(available_structures))
            # Make sure we have functionalisation sites
            while len(available_structures[chosen_structure]) == 0:
                chosen_structure = random.choice(list(available_structures))
            # Here's the actual structure
            base_structure = get_structure(chosen_structure)

            # Use several combinations to try to get something functionalised
            trial_number = 0
            while trial_number < max_trials:
                if trial_number < max_trials/4.0:
                    debug("Trial all groups: {}".format(trial_number))
                    status = random_combination_replace(
                        structure=base_structure, backends=backends,
                        max_different=2)
                elif trial_number < 2.0*max_trials/4.0:
                    debug("Trial max one group: {}".format(trial_number))
                    status = random_combination_replace(
                        structure=base_structure, backends=backends,
                        max_different=1)
                elif trial_number < 3.0*max_trials/4.0:
                    debug("Trial top 50: {}".format(trial_number))
                    status = random_combination_replace(
                        structure=base_structure, backends=backends,
                        groups_only=top_50_groups, max_different=2)
                else:
                    debug("Trial small groups: {}".format(trial_number))
                    status = random_combination_replace(
                        structure=base_structure, backends=backends,
                        groups_only=small_groups, max_different=1)

                # If functionalisation attempted
                if status:
                    if backends[0].cifs[-1]['functions']:
                        # it was successful; done here
                        break
                else:
                    # only increment if we actually tried to add groups
                    trial_number += 1
            else:
                site_replace(base_structure, replace_list=[],
                             backends=backends)
                failed = "{} random combinations".format(max_trials)


        # Should always have a structure, even if it is clean; but failed will
        # be True for that
        cif_info = backends[0].cifs[-1]
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
            available_structures=available_structures,
            failed=failed,
            **cif_info)
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


def get_structure(name):
    """
    Retrieve the structure from the cache or file. Return the structure object.

    Uses globals, so be careful
    """

    # Pick the structure
    # Update the cache and make room if it is too big
    if name in initialised_structures:
        base_structure = initialised_structures[name]
    else:
        base_structure = load_structure(name)
        initialised_structures[name] = base_structure
        if len(initialised_structures) > STRUCTURE_CACHE:
                    initialised_structures.popitem(last=True)

    return base_structure



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

    available_structures = OrderedDict()
    available_structure_files = sorted(glob.glob('*.cif'))

    # Pulling out sites only works with well formatted cifs with site label
    # as the first item on the line and ' H ' within the line for attachment
    # sites.
    # Since you are running the web component, you should know what
    # you're doing...
    for available_structure in available_structure_files:
        sname = available_structure[:-4]
        available_structures[sname] = []
        for line in open(available_structure):
            if ' H ' in line:
                available_structures[sname].append(line.split()[0])

    # Populate this later
    initialised_structures = OrderedDict()

    # Functional group library is self initialising
    info("Groups in library: %s" % str(functional_groups.group_list))

    #start tornado
    application.listen(TORNADO_PORT)
    info("Starting server on port number {}...".format(TORNADO_PORT))
    info("Open at http://127.0.0.1:{}/index.html".format(TORNADO_PORT))
    tornado.ioloop.IOLoop.instance().start()

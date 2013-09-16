#! /usr/bin/env python
import glob
import os
import json
import logging
import pickle
import random
import re
from os import path

import tornado.web
import tornado.ioloop
import tornado.template

import fapswitch
from fapswitch.bibliography import references
from fapswitch.config import debug, info
from fapswitch.core.io import load_structure
from fapswitch.functional_groups import functional_groups
from fapswitch.core.methods import random_combination_replace
from fapswitch.backend.web_store import WebStoreBackend

TORNADO_PORT = 8888

web_dir = path.join(path.dirname(fapswitch.__file__), 'web')

templates = tornado.template.Loader(path.join(web_dir, 'templates'))

# Send the index file
class IndexHandler(tornado.web.RequestHandler):
    def get(self, url = '/'):
        self.render('index.html')
    def post(self, url ='/'):
        self.render('index.html')


# Random generator webpage
class RandomHandler(tornado.web.RequestHandler):
    """
    Respond to a request for a random structure with a web view.
    """
    #both GET and POST requests have the same responses
    def get(self, url='/'):
        """GET request"""
        print("get")
        self.handle_request()

    def post(self, url='/'):
        """POST request"""
        print('post')
        self.handle_request()

    def handle_request(self):
        """Generate a random structure and return the rendered page."""

        info("{}".format(self))
        info("{}".format(self.request.arguments))
        groups = self.get_arguments('groups')
        debug(groups)
        if not groups:
            groups = None
        info("GROUPS !! {}".format(groups))

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
                # MEPO compliance if all groups are okay
                if all(functional_groups[function[0]].mepo_compliant for
                       function in cif_info['functions']):
                    mepo_compliant = "Yes"
                else:
                    mepo_compliant = "No"

                synthesis = """<h4>Step one:</h4>
                <p>Add some reagents and heat for a long time[Blank2017]</p>
                <h4>Step two:</h4>
                <p>More reagents and put in ice[Blank2015]<p>
                """

                local_references = []
                # This is always there for MEPO
                local_references.append(references['Kadantsev2013'])

                # Find all the references and add them too
                for reference in re.findall(r'\[(.*?)\]', synthesis):
                    local_references.append(references[reference])

                # Raw HTML anchors. Ugly.
                synthesis = re.sub(
                    r'\[(.*?)\]',  # non-greedy(?) find in square brackets
                    r'[<a href="#\1">\1</a>]',  # replace raw html
                    synthesis)


                page = templates.load('random.html').generate(
                    mepo_compliant=mepo_compliant,
                    references=local_references,
                    functional_groups=functional_groups,
                    synthesis=synthesis,
                    **cif_info)
                self.write(page)

                break
        else:
            page = templates.load('failed.html').generate()
            self.write(page)



# adds event handlers for commands and file requests
application = tornado.web.Application([
    (r"/(random.*)", RandomHandler ),
    (r"/", IndexHandler),
    (r"/(index\.html)", tornado.web.StaticFileHandler,{"path": web_dir}),
    (r"/(.*\.cif)", tornado.web.StaticFileHandler,{"path": path.join(web_dir, 'generated')}),
    (r"/(.*\.png)", tornado.web.StaticFileHandler,{"path": path.join(web_dir, 'static')}),
    (r"/(.*\.jpg)", tornado.web.StaticFileHandler,{"path": path.join(web_dir, 'static')}),
    (r"/(.*\.js)", tornado.web.StaticFileHandler,{"path": path.join(web_dir, 'js')}),
    (r"/(.*\.css)", tornado.web.StaticFileHandler,{"path": path.join(web_dir, 'css') }),
])


if __name__ == "__main__":

    #
    # Fun fapswitch as a web server
    #
    # Everything lives in the global scope

    # Directory to store pickles

    initial_directory = os.getcwd()

    datastore = os.path.join(web_dir, 'datastore')

    os.chdir(datastore)

    available_structures = glob.glob('*.cif')

    initialised_structures = []

    for available_structure in available_structures:
        sname = available_structure[:-4]
        initialised_structures.append(load_structure(sname))

    os.chdir(initial_directory)

    # Functional group library is self initialising
    info("Groups in library: %s" % str(functional_groups.group_list))

    #TODO: set collision tester

    #tell tornado to run checkSerial every 10ms
    #serial_loop = tornado.ioloop.PeriodicCallback(checkSerial, 10)
    #serial_loop.start()

    #start tornado
    application.listen(TORNADO_PORT)
    info("Starting server on port number {}...".format(TORNADO_PORT))
    info("Open at http://127.0.0.1:{}/index.html".format(TORNADO_PORT))
    tornado.ioloop.IOLoop.instance().start()

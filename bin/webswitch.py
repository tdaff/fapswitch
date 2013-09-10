#! /usr/bin/env python
import glob
import os
import json
import logging
import pickle
import random
from os import path

import tornado.web
import tornado.ioloop
import tornado.template

import fapswitch
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

        groups = self.get_argument('groups', None)
        debug(groups)
        if groups is not None:
            groups = groups.split(',')
        debug("{}".format(groups))

        max_trials = 10
        # Possible options:
        # replace_only: tuple of sites to replace
        # groups_only: only use specific groups
        # max_different: resrict simultaneous types of groups

        backends = [WebStoreBackend()]
        for _trial in range(max_trials):
            base_structure = random.choice(initialised_structures)
            status = random_combination_replace(groups_only=groups,
                                                structure=base_structure,
                                                backends=backends)
            if status:
                cif_info = backends[0].cifs[0]
                page = templates.load('random.html').generate(**cif_info)
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
    print("Starting server on port number %i..." % TORNADO_PORT)
    print("Open at http://127.0.0.1:%i/index.html" % TORNADO_PORT)
    tornado.ioloop.IOLoop.instance().start()

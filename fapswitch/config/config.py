#!/usr/bin/env python

"""
configuration for fapswitch

Provides the Options class that will transparently handle the different option
sources through the .get() method. Pulls in defaults, and job options plus
command line customisation. Instantiating Options will set up the logging for
the particular job.

"""

__all__ = ['Options']

import argparse
# Python 3 fix
try:
    import configparser
except ImportError:
    import ConfigParser as configparser
import copy
import logging
import os
import re
import sys

from io import StringIO
from os import path
from logging import debug, error


class DecreaseAction(argparse.Action):
    """Decrease the destination by 1 for each call."""
    def __call__(self, parser, namespace, values, option_string=None):
        """Decrease destination by one."""
        previous_value = getattr(namespace, self.dest)
        setattr(namespace, self.dest, previous_value-1)


class Options(object):
    """
    Transparent options handling.

    A single unified way of dealing with input files and command line options
    delivering sensible defaults for unspecified values. Access options with
    the .get() method, or the method that specifies the expected type.

    """
    def __init__(self, job_name=None):
        """Initialize options from all .fap files and the commandline."""
        # use .get{type}() to read attributes, only access args directly
        self.job_dir = ''
        self.script_dir = ''
        self.dot_faps = path.join(path.expanduser('~'), '.faps')
        self.job_name = job_name
        self.options = {}
        self.cmdopts = {}
        self.optfiles = configparser.SafeConfigParser()
        # populate options
        self._init_paths()
        self.commandline()
        self._init_logging()

        # Load in the defualts at the bottom
        self.load_without_sections(path.join(self.script_dir, 'defaults.fap'))

        # Specify any number of job options files
        for job_type in getattr(self.options, 'job_type', []):
            job_fap = '{}.fap'.format(job_type)
            self.load_without_sections(os.path.join(self.dot_faps, job_fap))

        # Job specific faps file
        self.load_without_sections(path.join(self.job_dir,
                                             '{}.fap'.format(self.job_name)))


    def get(self, item):
        """Map values from different sources based on priorities."""
        if item in self.__dict__:
            # Instance attributes, such as job_name and job_dir
            debug("an attribute: %s" % item)
            return object.__getattribute__(self, item)
        elif self.options.__dict__.get(item) is not None:
            # Commandline options from optparse where option is set
            debug("an option: %s" % item)
            return self.options.__dict__[item]
        elif item in self.cmdopts:
            # Commandline -o custom key=value options
            debug("a custom -o option: %s" % item)
            return self.cmdopts[item]
        elif self.optfiles.has_option('options', item):
            debug("a file option: %s" % item)
            return self.optfiles.get('options', item)
        else:
            # Most things have a default, but not always. Error properly.
            debug("unspecified option: %s" % item)
            raise AttributeError(item)

    def getbool(self, item):
        """
        Parse option and if the value of item is not already a bool return
        True for "1", "yes", "true" and "on" and False for "0", "no", "false"
        and "off". Case-insensitive.

        """
        value = self.get(item)
        if isinstance(value, bool):
            return value
        # Can't use isinstance with basestring to be 2.x and 3.x compatible
        # fudge it by assuming strings can be lowered
        elif hasattr(value, 'lower'):
            if value.lower() in ["1", "yes", "true", "on"]:
                return True
            elif value.lower() in ["0", "no", "false", "off"]:
                return False
            else:
                # Not a valid bool
                raise ValueError(value)
        else:
            return bool(item)

    def getint(self, item):
        """Return item's value as an integer."""
        value = self.get(item)
        return int(value)

    def getfloat(self, item):
        """Return item's value as a float."""
        value = self.get(item)
        return float(value)

    def gettuple(self, item, dtype=None):
        """Return item's value interpreted as a tuple of 'dtype' [strings]."""
        value = self.get(item)
        # Regex strips bracketing so can't nest, but safer than eval
        value = [x for x in re.split('[\s,\(\)\[\]]*', value) if x]
        if dtype is not None:
            return tuple([dtype(x) for x in value])
        else:
            return tuple(value)

    def _init_paths(self):
        """Find the script directory and set up working directory"""
        # Where the script is has the config defaults.
        if __name__ != '__main__':
            self.script_dir = os.path.dirname(__file__)
        else:
            self.script_dir = os.path.abspath(sys.path[0])
        # Where we run the job.
        self.job_dir = os.getcwd()

    def _init_logging(self):
        """
        Setup the logging to terminal and .flog file, with levels as required.
        Must run before any logging calls so we need to access attributes
        rather than using self.get()!

        """
        root_logger = logging.getLogger()

        # Have the logger itself set with the lowest possible level
        root_logger.setLevel(logging.DEBUG)
        # Reset any handlers that might have been set accidentally
        root_logger.handlers = []

        # Always at least INFO in .flog
        file_level = logging.INFO
        flog_filename = '{}.flog'.format(self.job_name)

        verbosity = self.options.verbosity

        if verbosity <= -2:
            # -qq
            stdout_level = logging.CRITICAL
        elif verbosity <= -1:
            # -q
            stdout_level = logging.ERROR
        elif verbosity >= 1:
            # -v
            stdout_level = logging.DEBUG
            file_level = logging.DEBUG
        else:
            stdout_level = logging.INFO

        # Easier to do simple file configuration then add the stdout
        file_handler = logging.FileHandler(flog_filename)
        file_handler.setLevel(file_level)
        formatter = logging.Formatter('[%(asctime)s] %(levelname)s %(message)s',
                                      datefmt='%Y%m%d %H:%M:%S')
        file_handler.setFormatter(formatter)
        root_logger.addHandler(file_handler)

        # Make these uniform widths
        logging.addLevelName(10, '--')
        logging.addLevelName(20, '>>')
        logging.addLevelName(30, '**')
        logging.addLevelName(40, '!!')
        logging.addLevelName(50, 'XX')

        # Use nice coloured console output
        console = ColouredConsoleHandler(stream=sys.stdout)
        console.setLevel(stdout_level)
        formatter = logging.Formatter('%(levelname)s %(message)s')
        console.setFormatter(formatter)
        # add the handler to the root logger
        root_logger.addHandler(console)


    def commandline(self):
        """Specified options, highest priority."""
        DESCRIPTION = """Fapswitch function switching program"""
        # use description for the script, not for this module
        parser = argparse.ArgumentParser(description=DESCRIPTION)
        parser.add_argument("-v", "--verbose", action="count", default=0,
                            dest="verbosity", help="Increase verbosity. "
                            "Specify more times for more debugging "
                            "information. Cancels '--quiet'.")
        parser.add_argument("-q", "--quiet", action=DecreaseAction, nargs=0,
                            dest="verbosity", help="Decrease verbosity. "
                            "Specify more times for less output. Cancels "
                            "'--verbose'.")
        parser.add_argument("-o", "--option", action="append", dest="cmdopts",
                            default=[], help="Set program options as "
                            "section.key=value pairs. Use \"quotation marks\" "
                            "if options contain spaces.")
        parser.add_argument("-j", "--job-type", dest="job_type",
                            action="append", default=[], help="Read "
                            "preconfigured job settings from job-type.fap in "
                            "the user ~/.faps/ directory")
        # Always have the job name at the end
        parser.add_argument('job_name', help="Name for job", nargs='?',
                            default='default')
        local_args = parser.parse_args()

        if self.job_name is None:
            self.job_name = local_args.job_name

        # key value options from the command line
        for pair in local_args.cmdopts:
            if '=' in pair:
                pair = pair.split('=')
                self.cmdopts[pair[0]] = pair[1]
            else:
                self.cmdopts[pair] = True

        self.options = local_args

    def load_without_sections(self, filename):
        """Load a configuration with no header into the options."""
        # ConfigParser requires header sections so we add them to a StringIO
        # of the file if they are missing.
        try:
            with open(filename, 'r') as filetemp:
                file_contents_fp = StringIO('[options]\n' + filetemp.read())
            self.optfiles.readfp(file_contents_fp)
            debug('Incorporated options from: {}'.format(filename))
        except IOError:
            # file does not exist so we just use a blank string
            debug('Options file not found: {}'.format(filename))


class ColouredConsoleHandler(logging.StreamHandler):
    """Makes colourised output for the console."""
    def emit(self, record):
        """Colourise leve id and emit a record."""
        # Need to make a actual copy of the record
        # to prevent altering the message for other loggers
        myrecord = copy.copy(record)
        levelno = myrecord.levelno
        if levelno >= 50:  # CRITICAL / FATAL
            front = '\033[30;41m'  # black/red
        elif levelno >= 40:  # ERROR
            front = '\033[30;41m'  # black/red
        elif levelno >= 30:  # WARNING
            front = '\033[30;43m'  # black/yellow
        elif levelno >= 20:  # INFO
            front = '\033[30;42m'  # black/green
        elif levelno >= 10:  # DEBUG
            front = '\033[30;46m'  # black/cyan
        else:  # NOTSET and anything else
            front = '\033[0m'  # normal

        myrecord.levelname = '{}{}\033[0m'.format(front, myrecord.levelname)
        logging.StreamHandler.emit(self, myrecord)


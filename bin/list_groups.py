#!/usr/bin/env python

"""
list_groups.py

Output a list of all the functional groups in the libraries
supplied to fapswitch in the current context.

"""


import sys
from os.path import dirname, realpath
sys.path.insert(1, dirname(dirname(realpath(__file__))))

from fapswitch.functional_groups import functional_groups


def main():
    """
    Just tell each group to log itself to the outputs.

    """

    for group_name in functional_groups:
        functional_groups[group_name].log_info()


if __name__ == '__main__':

    main()

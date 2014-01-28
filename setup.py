#!/usr/bin/env python

from distutils.core import setup
from glob import glob

scripts = glob('bin/*')

setup(name='fapswitch',
      version='7.1.2',
      description='Functional group modification on CIF files',
      long_description=open('README.rst').read(),
      author='Tom Daff',
      author_email='tdaff@uottawa.ca',
      license='BSD',
      url='http://titan.chem.uoattawa.ca/fapswitch/',
      packages=['fapswitch', 'fapswitch.core', 'fapswitch.extensions',
                'fapswitch.bibliography', 'fapswitch.config',
                'fapswitch.functional_groups', 'fapswitch/backend'],
      package_data={'fapswitch':
                        ['config/*.fap', 'functional_groups/library/*.flib',
                         'web/templates/*html', 'web/css/*.css',
                         'web/datastore/*cif']},
      scripts=scripts,
      classifiers=["Programming Language :: Python",
                   "Programming Language :: Python :: 3",
                   "Development Status :: 5 - Production/Stable",
                   "Intended Audience :: Science/Research",
                   "License :: OSI Approved :: BSD License",
                   "Operating System :: OS Independent",
                   "Topic :: Scientific/Engineering :: Chemistry"])

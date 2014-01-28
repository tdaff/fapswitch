#!/usr/bin/env python

from distutils.core import setup
from glob import glob

try:
    import Cython.Build
    from Cython.Distutils import build_ext
except ImportError:
    use_cython = False
else:
    use_cython = True

cmdclass = {}
ext_modules = []

if use_cython:
    ext_modules = Cython.Build.cythonize(['fapswitch/*.py', 'fapswitch/*/*.py', 'fapswitch/*/*/*.py'], exclude=['fapswitch/extensions/fragments.py'])
    cmdclass.update({ 'build_ext': build_ext })

scripts = glob('bin/*')

setup(name='fapswitch',
      version='7.1.1',
      description='Functional group modification on cifs',
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
      cmdclass=cmdclass,
      ext_modules=ext_modules,
      classifiers=["Programming Language :: Python",
                   "Programming Language :: Python :: 3",
                   "Development Status :: 5 - Production/Stable",
                   "Intended Audience :: Science/Research",
                   "License :: OSI Approved :: BSD License",
                   "Operating System :: OS Independent",
                   "Topic :: Scientific/Engineering :: Chemistry",
                   ])

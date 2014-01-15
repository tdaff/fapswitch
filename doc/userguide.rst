User Guide
==========

Installation
------------

Requirements
############

Fapswitch has the following prequisites, this is the minimum required to
run the code:

* Python

  * Version 2.7, 3.3 or greater
  * http://www.python.org/


* NumPy

  * Any recent version
  * http://www.numpy.org/


To build the documentation you will need:

* Sphinx

  * Any version
  * http://sphinx-doc.org/


Added functionality is accessible with these optional libraries but the
code will function without them.

* Open Babel

  * Required for ligand structure perception
  * http://openbabel.org


* Periodic Open Babel

  * Required for bond perception in initial structure
  * http://software-lisc.fbk.eu/obgmx/


* RDKit

  * Required for synthetic accessibility scores
  * http://www.rdkit.org/


* SQLAlchemy

  * Store cifs in the database backend in addition to files
  * Useful if you generate lots of structures
  *  http://www.sqlalchemy.org/


* SciPy

  * A `scipy.weave` optimised collision tester, `cvdw`, can speed up the code
  * Currently does not work with Python 3.3 or where no C++ compiler is found
  * http://www.scipy.org/


* Tornado

  * Used to provide a web interface to functionalisation
  * http://www.tornadoweb.org/


Source
######

* Obtain a copy of the source code from https://bitbucket.org/tdaff/fapswitch.
  If you do not have access, contact a member of the Woo lab.
* Unpack the source somewhere.
* Add the `bin/` to your path to have access to the scripts for running
  the code.
* Add the root directory to your PYTHONPATH to make the modules importable
  from your own code.


Running on Wooki
----------------

* Use the provided module for fapswitch.
* Run `module avail` to see which versions are available
* Use `module load fapswitch/7.0.0` where 7.0.0 can be any later version
* The module will try and load the modules for RDKit

Running on Sharcnet
-------------------

* The current version is maintained in `tdaff`'s home
* Add the following to your `~/.bashrc`:

  .. code-block:: sh

    # for fapswitch
    export PATH=$PATH:/home/tdaff/Compile/fapswitch-7.0/bin
    export PYTHONPATH=/home/tdaff/Compile/fapswitch-7.0:$PYTHONPATH

    # for openbabel/pybel
    export PYTHONPATH=$PYTHONPATH:/home/tdaff/work/Software/Python/lib/python2.7/site-packages
    export LD_LIBRARY_PATH=/home/tdaff/work/Software/lib:$LD_LIBRARY_PATH
    export PATH=/home/tdaff/work/Software/bin:$PATH

    # for RDKit
    export LD_LIBRARY_PATH=/home/tdaff/work/Software/RDKit_2013_09_1/lib:$LD_LIBRARY_PATH
    export PYTHONPATH=/home/tdaff/work/Software/RDKit_2013_09_1:$PYTHONPATH
    export LD_PRELOAD=/usr/lib64/libstdc++.so.6


Preparing Input Files
---------------------

Fapswitch requires your initial structure as a cif file. To make sure you
have the best functionalisation experience, the following preparations
can help you

* Clean any solvent from the structure, otherwise it will prevent
  functional groups fitting in.
* Remove disorder from the structure. During the functionalisation,
  partially occupied sites all become fully occupied.
* Reduce the atoms in the cif to the asymmetric unit and make sure that
  all the symmetry operations are provided.
* (Optional) Add bonding information to the cif:

  * CIFs exported from Materials Studio with the correct bonding will have
    the bonding in the correct format for fapswitch.
  * To manually add the bonding, it is contained within a `loop_` that defines
    atoms in the bond by their atom labels, optionally their distance, and
    the bond type (S=single, D=double, T=Triple, A=Aromatic)::

       loop_
       _geom_bond_atom_site_label_1
       _geom_bond_atom_site_label_2
       _geom_bond_distance
       _ccdc_geom_bond_type
       C1 C2 1.1021 A

* (Optionally) Add atom typing information to the cif:

  * Put the atom type in the `_atom_site_description` field.
  * Assign atom types based on the UFF to work with the atom types on the
    functional groups.
  * These get passed through to the output file.
  * Must be added manually::

      loop_
      _atom_site_label
      _atom_site_type_symbol
      _atom_site_description
      _atom_site_fract_x
      _atom_site_fract_y
      _atom_site_fract_z
      C1    C     C_R   0.943990 0.413360 0.732350
      C2    C     C_R   0.681010 0.211640 0.267650
      O1    O     O_3   1.000000 0.324853 0.517757
      H1    H     H_    0.881320 0.412480 0.761750

Options
-------

All of the fapswitch codes use a unified options and command scheme. Options
are set as simple "option = value". In some cases the value can be a list of
values, for which punctuation is ignored, so that any whitespace, comma or
brackets can be used to separate items. Boolean options accept, "true", "yes",
"on", "1" for positive responses (case insensitive). Options
specified are used in ascending priority as listed below:

* (lowest priority) **defaults** are set for all options. These are shown below
  or can be found in the file `fapswitch/config/defaults.fap`
* **job types** are specified through the commandline. When running the
  command use the `-j` option to specify job types and the code will search
  the user's ~/.faps folder for the corresponding job.fap files, for example:

    cliswitch.py -j basic -j just_halogen job_name

  will read the options from `~/.faps/basic.fap`, followed by
  `~/.faps/just_halogen.fap`
* **job specific** options are read from `job_name.fap` in the working
  directory.
* (highest priority) **commandline** options are specified with
  `-o option=value` and take priority over all other options.

All options
###########

.. envvar:: backends

  Default: file

  Backends to store the output structures. [str, list] {file, sqlite}

.. envvar:: collision_method

  Default: vdw

  Method to use to test for collisions [str] {absolute, covalent, vdw}

.. envvar:: collision_scale

  Default: 1.122462048309373

  Absolute value in Angstrom or scale factor for atomic radii for minimum
  distance in collision test for insertions [float]

.. envvar:: connectivity

  Default: openbabel

  Where to get the connectivity information from [str] {openbabel, file}

.. envvar:: custom_strings

  Default:

  Make functionalisations with the set of {.freeform.srings.} and
  [symm@try.strings]. [str, list]

.. envvar:: full_random_count

  Default: 0

  Number of completely randomised structures to make. [int]

.. envvar:: max_different

  Default: 0

  Maximum number of groups that will be used simultaneously. [int]

.. envvar:: mepo_only

  Default: False

  Only load MEPO-QEq compatible groups. [bool]

.. envvar:: port

  Default: 0

  Socket port to run the server mode on; leave as zero to pick random
  available port as two instances cannot share a port. [int]

.. envvar:: replace_all_sites

  Default: False

  Should fapswitch produce all group@site combinations? [bool]

.. envvar:: replace_groups

  Default:

  Only use the specified groups in systematic functionalisations. [list]

.. envvar:: replace_only

  Default:

  Only replace the listed sites in systematic functionalisations. [list]

.. envvar:: site_random_count

  Default: 0

  Number of symmetry based randomised structures to make. [int]

.. envvar:: timeout

  Default: 7200

  Number of seconds of inactivity after which the fapswitchd.py
  daemon will close[int]

.. envvar:: unfunctionalised_probability

  Default: 0.5

  Probability that a site will have no functionalisation in random switching
  scheme. [float]


Running the code
----------------

For a typical job, you would set the required options in the configuration
files and run `cliswitch.py job_name` to generate functionalised structures.

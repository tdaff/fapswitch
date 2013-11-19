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




Preparing Input Files
---------------------

Fapswitch requires your initial structure as a cif file. To make sure you
have the best functionalisation experience, the following preparations
can help you

* Clean any solvent from te structure, otherwise it will prevent
  functional groups fitting in.
* Remove disorder from the structure. During the functionalisation,
  partially occupied sites all become fully occupied.
* Reduce the atoms in the cif to the asymmetric unit and make sure that
  all the symmetry operations are provided.
* (Optional) Add bonding information to the cif

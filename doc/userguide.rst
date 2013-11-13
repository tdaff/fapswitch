User Guide
==========

Installation
------------

Fapswitch has the following prequisites, this is the minimum required to
run the code:

* Python

  * Version 2.7, 3.3 or greater
  * http://www.python.org/


* NumPy

  * Any recent version
  * http://www.numpy.org/


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






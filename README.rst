Fapswitch
---------

Functional groups switching procedure targeted at hypothetical MOF
generation from high symmetry cif files.

This code is the property of the Woo Lab and must not be used without
contacting a member of the lab first.

Documentation
=============

This code is currently undocumented. See the '--help' and the input
files for more information.

Installation
============

This code is currently undocumented.

Running
=======

This code is currently undocumented.

Changes
=======

7.0.0
------
Version 7 is a complete reorganisation of the code with many
compatibility breaking changes

  * ``fapswitch_`` prefix dropped from all the options.
  * separate scripts for different use cases:

    * ``cliswitch.py`` -- for all combination replacement or interactive use
    * ``fapswitchd.py`` -- start a socket server that accepts functional strings
    * ``webswitch.py`` -- tornado based webapp to demonstrate functionalisation

  * ``functional_groups`` and ``options`` are now in globally accessible modules
  * ``-q`` and ``-v`` cancel each other.
  * ``--silent`` option has been removed, ``-qq`` gives the same effect
  * Collision detection now properly uses 1/2 sum of vdW radii

  * User defined socket timeout


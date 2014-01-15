Changelog
=========

7.1.0
-----

  * Functional groups have been cleaned and assigned IUPAC names
  * MEPO-QEq compatibility has been computed for all functional groups
  * ``mepo_only`` option will only load MEPO compatible groups
  * ``list_groups.py`` lists all the groups that have been loaded in the 
    current context
  * In the API mepo_compliance is renamed to mepo_compatibility
  * Documentation includes the functional groups


7.0.0
------
Version 7 is a complete reorganisation of the code with many compatibility
breaking changes. Future development will be based on this codebase.

  * ``fapswitch_`` prefix dropped from all the options.
  * separate scripts for different use cases:

    * ``cliswitch.py`` -- for all combination replacement or interactive use
    * ``fapswitchd.py`` -- start a socket server that accepts functional strings
    * ``webswitch.py`` -- tornado based webapp to demonstrate functionalisation

  * ``functional_groups`` and ``options`` are now in globally accessible modules
  * ``-q`` and ``-v`` cancel each other.
  * ``--silent`` option has been removed, ``-qq`` gives the same effect
  * output text is no longer automatically wrapped
  * Collision detection now properly uses 1/2 sum of vdW radii

  * User defined socket timeout


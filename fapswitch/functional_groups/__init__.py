"""

fapswitch.functional_groups

Provides the functional_groups library object that is initialised
to contain all the groups available.

"""

from .group_reader import FunctionalGroupLibrary

# instance the library here; the module should make
# the same groups available wherever it is imported
functional_groups = FunctionalGroupLibrary()

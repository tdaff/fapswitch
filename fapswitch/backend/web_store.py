"""
Store backend for webswitch that keeps files in the cif_store directory and
tracks their filenames in object.cifs.

"""

import errno
import hashlib
import os
import os.path

_generated_path = os.path.join('web', 'generated')

class WebStoreBackend(object):
    """Abstraction for keeping cif file in memory to be retrieved."""

    def __init__(self, directory=_generated_path):
        """
        Make a storage containing backend and setup a directory.
        """
        self.cifs = []
        self.directory = directory
        try:
            os.makedirs(self.directory)
        except OSError as exc:
            if exc.errno == errno.EEXIST and os.path.isdir(self.directory):
                pass
            else:
                raise

    def add_symmetry_structure(self, base_structure, functions, cif_file):
        """
        Write out the cif file with a name derived from the base structure
        and the functionalisations and track the name.

        """
        new_mof_name = ".".join(["@".join(x) for x in functions])
        cif_filename = '%s_func_%s.cif' % (base_structure, new_mof_name)
        cif_save_path = os.path.join(self.directory, cif_filename)
        with open(cif_save_path, 'w') as output_file:
            output_file.writelines(cif_file)
        self.cifs.append(
            {'mof_name': new_mof_name,
             'base_structure': base_structure,
             'functions': functions,
             'cif_filename': cif_filename})

    def add_freeform_structure(self, base_structure, functions, cif_file):
        """
        Write a cif file with an md5 fixed-length name based on the
        functionalisation.

        """
        unique_name = hashlib.md5(str(functions)).hexdigest()
        cif_filename = os.path.join(
            self.directory,
            '%s_free_%s.cif' % (base_structure, unique_name))
        with open(cif_filename, 'w') as output_file:
            output_file.writelines(cif_file)
        self.cifs.append(
            {'mof_name': unique_name,
             'base_structure': base_structure,
             'functions': functions,
             'cif_filename': cif_filename})

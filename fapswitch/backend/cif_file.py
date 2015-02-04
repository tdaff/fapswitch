"""
Cif file writing backend for fapswitch.

"""

import hashlib

from fapswitch.config import options
from fapswitch.config import error, debug


class CifFileBackend(object):
    """Abstraction for writing cif files in a pluggable manner."""

    def __init__(self):
        self.hash_filenames = options.get('hash_filenames').lower()

    def add_symmetry_structure(self, base_structure, functions, cif_file,
                               manual_angles=None, **kwargs):
        """
        Write out the cif file with a name derived from the base structure
        and the functionalisations.

        """
        if manual_angles is not None and any(manual_angles):
            new_mof_components = []
            for site, manual_angle in zip(functions, manual_angles):
                if manual_angle is not None:
                    angle_str = "%{}".format(manual_angle)
                else:
                    angle_str = ""
                new_mof_components.append("{}@{}{}".format(site[0], site[1],
                                                           angle_str))
            new_mof_name = ".".join(new_mof_components)
        else:
            new_mof_name = ".".join(["@".join(x) for x in functions])

        hashed_name = hashlib.md5(new_mof_name).hexdigest()

        if self.hash_filenames == 'always':
            cif_filename = '%s_func_%s.cif' % (base_structure, hashed_name)
            with open(cif_filename, 'w') as output_file:
                output_file.writelines(cif_file)
        else:
            cif_filename = '%s_func_%s.cif' % (base_structure, new_mof_name)
            try:
                with open(cif_filename, 'w') as output_file:
                    output_file.writelines(cif_file)
            except IOError:
                if self.hash_filenames == 'never':
                    error('Unable to write file {}'.format(cif_filename))
                    return 1
                # 'onerror', or any other has option leads here; try and use
                # a shorter name
                cif_filename = '%s_func_%s.cif' % (base_structure, hashed_name)
                debug('Automatically shortened file name to '
                      '{}'.format(cif_filename))
                with open(cif_filename, 'w') as output_file:
                    output_file.writelines(cif_file)

    def add_freeform_structure(self, base_structure, functions, cif_file,
                               **kwargs):
        """
        Write a cif file with an md5 fixed-length name based on the
        functionalisation.

        """
        unique_name = hashlib.md5(str(functions).encode('utf-8')).hexdigest()
        cif_filename = '%s_free_%s.cif' % (base_structure, unique_name)
        with open(cif_filename, 'w') as output_file:
            output_file.writelines(cif_file)

"""
Cif file writing backend for fapswitch.

"""

import hashlib


class CifFileBackend(object):
    """Abstraction for writing cif files in a pluggable manner."""

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
        cif_filename = '%s_func_%s.cif' % (base_structure, new_mof_name)
        with open(cif_filename, 'w') as output_file:
            output_file.writelines(cif_file)

    def add_freeform_structure(self, base_structure, functions, cif_file, **kwargs):
        """
        Write a cif file with an md5 fixed-length name based on the
        functionalisation.

        """
        unique_name = hashlib.md5(str(functions).encode('utf-8')).hexdigest()
        cif_filename = '%s_free_%s.cif' % (base_structure, unique_name)
        with open(cif_filename, 'w') as output_file:
            output_file.writelines(cif_file)

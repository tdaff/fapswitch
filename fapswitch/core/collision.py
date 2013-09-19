#!/usr/bin/env python

"""
ccollision.py

SciPy weave compiled collision testers.

"""


#from scipy import weave
#from scipy.weave import converters

from numpy import dot

from fapswitch.core.util import min_vect
from fapswitch.config import debug, info, warning, error, critical




def make_collision_tester(options=None, test_method=None, test_scale=None):
    """
    Create a function that will test atom overlap based on separation
    and atomic radii.

    """

    if test_method is None:
        test_method = options.get('collision_method').lower()
    if test_scale is None:
        test_scale = options.getfloat('collision_scale')

    if test_method == 'covalent':
        info('Covalent radii collision test, scale factor: %f' % test_scale)
        def collision(test_atom, atoms, cell, ignore=()):
            """Covalent radii collision test."""
            # TODO(tdaff): Can't use cleaner .fractional tests as the structure
            # is not passed around to make _parent; can this be fixed
            pos = test_atom.ipos(cell.cell, cell.inverse)
            ipos = test_atom.ifpos(cell.inverse)
            for idx, atom in enumerate(atoms):
                # Skip non atoms
                if atom is None or idx in ignore:
                    continue
                dist = min_vect(pos, ipos, atom.ipos(cell.cell, cell.inverse),
                                atom.ifpos(cell.inverse), cell.cell)
                dist = (dot(dist, dist))**0.5
                min_dist = test_scale*(test_atom.covalent_radius + atom.covalent_radius)
                if dist < min_dist:
                    return False
            return True

    elif test_method == 'vdw':
        info('VdW radii collision test, scale factor: %f' % test_scale)
        half_scale = test_scale*0.5
        def collision(test_atom, atoms, cell, ignore=()):
            """Covalent radii collision test."""
            pos = test_atom.ipos(cell.cell, cell.inverse)
            ipos = test_atom.ifpos(cell.inverse)
            for idx, atom in enumerate(atoms):
                # Skip non atoms
                if atom is None or idx in ignore:
                    continue
                dist = min_vect(pos, ipos, atom.ipos(cell.cell, cell.inverse),
                                atom.ifpos(cell.inverse), cell.cell)
                dist = (dot(dist, dist))**0.5
                min_dist = half_scale*(test_atom.vdw_radius + atom.vdw_radius)
                if dist < min_dist:
                    return False
            return True

    elif test_method == 'cvdw':
        from ccollision import wdist
        info('CVdW radii collision test, scale factor: %f' % test_scale)
        half_scale = test_scale*0.5
        def collision(test_atom, atoms, cell, ignore=()):
            """Covalent radii collision test."""
            pos = test_atom.ipos(cell.cell, cell.inverse)
            ipos = test_atom.ifpos(cell.inverse)
            for idx, atom in enumerate(atoms):
                # Skip non atoms
                if atom is None or idx in ignore:
                    continue
                dist = wdist(pos, ipos, atom.ipos(cell.cell, cell.inverse),
                             atom.ifpos(cell.inverse), cell.cell)
                min_dist = half_scale*(test_atom.vdw_radius + atom.vdw_radius)
                if dist < min_dist:
                    return False
            return True
    else:
        info('Collison test with absolute distance: %f' % test_scale)
        def collision(test_atom, atoms, cell, ignore=()):
            """Covalent radii collision test."""
            pos = test_atom.ipos(cell.cell, cell.inverse)
            ipos = test_atom.ifpos(cell.inverse)
            for idx, atom in enumerate(atoms):
                # Skip non atoms
                if atom is None or idx in ignore:
                    continue
                dist = min_vect(pos, ipos, atom.ipos(cell.cell, cell.inverse),
                                atom.ifpos(cell.inverse), cell.cell)
                dist = (dot(dist, dist))**0.5
                if dist < test_scale:
                    return False
            return True

    # Make a closure for the tester function
    return collision



def wdist(c_coa, f_coa, c_cob, f_cob, box):
    """Rewritten min_vect for weave speedup."""
    code = """
    int i, j, k;
    double p, q, r;
    double distance;
    double diff[3];
    double new_f_cob[3];
    int no_min_image;

    no_min_image = 0;

    for(i=0;i<3;i++)
    {
        p = f_coa[i];
        q = f_cob[i];
        diff[i] = p - q;
        if (diff[i] < -0.5)
        {
            new_f_cob[i] = q - 1.0;
        }
        else
        {
            if (diff[i] > 0.5)
            {
                new_f_cob[i] = q + 1.0;
            }
            else
            {
                new_f_cob[i] = q;
                no_min_image++;
            };
        };
    };

//    if (diff[0] != 0.0 || diff[1] != 0.0 || diff[2] != 0.0)
    if (no_min_image < 3)
    {
//        std::cout << diff[0] << diff[1] << diff[2] << std::endl;

        distance = 0;

        for(i=0;i<3;i++)
        {
            p = new_f_cob[0]*box(0, i) + new_f_cob[1]*box(1, i) + new_f_cob[2]*box(2, i);
            q = c_coa(i);
            distance = distance + pow(p - q, 2);
        };

    }
    else
    {
//        std::cout << diff[0] << std::endl;

        distance = 0;

        for(i=0;i<3;i++)
        {
            p = c_cob(i);
            q = c_coa(i);
            distance = distance + pow(p - q, 2);
        };
    };

    return_val = sqrt(distance);
    """
    return weave.inline(code, ['c_coa', 'f_coa', 'c_cob', 'f_cob', 'box'],
                        type_converters=converters.blitz,
                        support_code='#include <math.h>')

"""

Small utility functions for function switching
e.g. with basic maths and geometry small functions.

"""

from itertools import chain, combinations

from numpy import array, asarray, outer, cos, identity, sin
from numpy import cross, roll
from numpy.linalg import norm


def rotation_about_angle(axis_in, angle):
    """
    Create a rotation matrix corresponding to the rotation around a general
    axis by a specified angle.

    """
    axis = axis_in/norm(axis_in)

    aprod = outer(axis, axis)
    skew = array([[0, axis[2], -axis[1]], [-axis[2], 0, axis[0]],
        [axis[1], -axis[0], 0]])

    # R = dd^T + cos(a) (I - dd^T) + sin(a) skew(d)
    return aprod+cos(angle)*(identity(3)-aprod)+sin(angle)*skew


def arbitrary_normal(vector):
    """Create a normalised normal to an input, does not use random values."""
    if vector == [0, 0, 0]:
        # everything is norma to a zero length vector
        return asarray([1, 0 , 0])
    elif vector[0] == vector[1] == vector[2]:
        # just need cross product with something with non-equal elements
        return normalise(cross(vector, [0.0, 1.0, 0.0]))
    else:
        #different elements so we can swap some to create a different vector
        return normalise(cross(vector, roll(vector, 1)))


def min_vect(c_coa, f_coa, c_cob, f_cob_in, box):
    """Calculate the closest distance assuming fractional, in-cell coords."""
    f_cob = f_cob_in[:]
    fdx = f_coa[0] - f_cob[0]
    if fdx < -0.5:
        f_cob[0] -= 1
    elif fdx > 0.5:
        f_cob[0] += 1
    fdy = f_coa[1] - f_cob[1]
    if fdy < -0.5:
        f_cob[1] -= 1
    elif fdy > 0.5:
        f_cob[1] += 1
    fdz = f_coa[2] - f_cob[2]
    if fdz < -0.5:
        f_cob[2] -= 1
    elif fdz > 0.5:
        f_cob[2] += 1
    if f_cob == f_cob_in:
        # if nothing has changed, use initial values
        return direction3d(c_coa, c_cob)
    else:
        new_b = [f_cob[0]*box[0][0] + f_cob[1]*box[1][0] + f_cob[2]*box[2][0],
                 f_cob[0]*box[0][1] + f_cob[1]*box[1][1] + f_cob[2]*box[2][1],
                 f_cob[0]*box[0][2] + f_cob[1]*box[1][2] + f_cob[2]*box[2][2]]
        return direction3d(c_coa, new_b)

def direction3d(source, target):
    """Return the vector connecting two 3d points."""
    return [target[0] - source[0],
            target[1] - source[1],
            target[2] - source[2]]


def normalise(vector):
    """Return an array with magnitude 1."""
    return asarray(vector)/norm(vector)


def vecdist3(coord1, coord2):
    """Calculate vector between two 3d points."""
    #return [i - j for i, j in zip(coord1, coord2)]
    # Twice as fast for fixed 3d vectors
    vec = [coord2[0] - coord1[0],
           coord2[1] - coord1[1],
           coord2[2] - coord1[2]]

    return (vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2])**0.5


# Itertools derived functions

def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))


# List functions

def strip_blanks(lines):
    """Strip lines and remove blank lines."""
    return [line.strip() for line in lines if line.strip() != '']


def subgroup(iterable, width, itype=None):
    """Split an iterable into nested sub-itypes of width members."""
    # Return the same type as iterable
    if itype is None:
        if isinstance(iterable, list):
            itype = list
        else:
            itype = tuple
    # Will leave short groups if not enough members
    return itype([itype(iterable[x:x+width])
                  for x in range(0, len(iterable), width)])

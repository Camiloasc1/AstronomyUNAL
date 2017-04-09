import numpy as np


def cross(a, b):
    return np.cross(a, b) * a.unit * b.unit
    # x = a[1] * b[2] - a[2] * b[1]
    # y = a[2] * b[0] - a[0] * b[2]
    # z = a[0] * b[1] - a[1] * b[0]
    # return np.array([x, y, z], dtype=object)


def norm(a):
    return np.linalg.norm(a) * a.unit


def angle_orbit(a, e, r, r_, r1_):
    from CelestialMechanics.orbits import ellipse, parable, hyperbola

    if e < 1:  # Ellipse
        angles = ellipse.angles(a, e, r)
    elif e > 1:  # Hyperbola
        angles = hyperbola.angles(a, e, r)
    elif e == 1:  # Parable
        angles = parable.angles(a, r)

    return angles[0 if r1(r_, r1_) > 0 else 1]


def angle_cuadrant(angles, r_, r1_):
    return angles[0 if r1(r_, r1_) > 0 else 1]


def r1(r_, r1_):
    """
    r. = dot(r_, r1_) / norm(r_)

    :param r1_: r1 vector
    :type r1_: (float, float, float)
    :param r_: r vector
    :type r_: (float, float, float)
    :return: r.
    :rtype: float
    """
    return np.dot(r_, r1_) / norm(r_)

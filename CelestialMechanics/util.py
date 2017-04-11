from typing import List, Tuple

import numpy as np


def cross(a: List[float], b: List[float]) -> List[float]:
    """
    np.cross(a, b) keeping the units

    :param a: vector a
    :type a: list
    :param b: vector b
    :type b: list
    :return: np.cross(a, b)
    :rtype: list
    """
    return np.cross(a, b) * a.unit * b.unit


def dot(a: List[float], b: List[float]) -> float:
    """
    np.dot(a, b) keeping the units

    :param a: vector a
    :type a: list
    :param b: vector b
    :type b: list
    :return: np.dot(a, b)
    :rtype: float
    """
    return np.dot(a, b) * a.unit * b.unit


def norm(a: List[float]) -> float:
    """
    np.linalg.norm(a) keeping the units

    :param a: vector a
    :type a: list
    :return: np.linalg.norm(a)
    :rtype: float
    """
    return np.linalg.norm(a) * a.unit


def angle_orbit(a: float, e: float, r: float, r_: List[float], r1_: List[float]) -> float:
    """
    the correct angle given the orbital elements and velocity

    :param a: semi-major axis
    :type a: float
    :param e: eccentricity
    :type e: float
    :param r: radius vector
    :type r: float
    :param r_: radius vector
    :type r_: (float, float, float)
    :param r1_: velocity vector
    :type r1_: (float, float, float)
    :return:
    :rtype:
    """
    from CelestialMechanics.orbits import ellipse, parable, hyperbola

    if e < 1:  # Ellipse
        angles = ellipse.angles(a, e, r)
    elif e > 1:  # Hyperbola
        angles = hyperbola.angles(a, e, r)
    elif e == 1:  # Parable
        angles = parable.angles(a, r)

    return angles[0 if r1(r_, r1_) > 0 else 1]


def angle_cuadrant(angles: Tuple[float, float], r_: List[float], r1_: List[float]) -> float:
    """
    the correct angle

    :param angles: the two possible angles
    :type angles: (float, float)
    :param r_: radius vector
    :type r_: (float, float, float)
    :param r1_: velocity vector
    :type r1_: (float, float, float)
    :return: correct theta angle
    :rtype: float
    """
    return angles[0 if r1(r_, r1_) > 0 else 1]


def r1(r_: List[float], r1_: List[float]) -> float:
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

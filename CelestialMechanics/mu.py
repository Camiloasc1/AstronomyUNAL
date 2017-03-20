import numpy as np

from CelestialMechanics.kepler.constants import K


def mu_sun(m2_over_m1: float) -> float:
    """
    mu = k * sqrt(1 + m2/m1)

    :param m2_over_m1:
    :type m2_over_m1:
    :return: mu
    :rtype: float
    """

    mu = K * np.sqrt(1. + m2_over_m1)
    return mu * mu


def mu_na(n: float, a: float) -> float:
    """
    mu = n^2 / a^3

    :param n: mean motion in degrees
    :type n: float
    :param a: semi-major axis
    :type a: float
    :return: mu
    :rtype: float
    """
    return n * n * a * a * a


def mu_gm1m2(m1: float, m2: float) -> float:
    """
    mu = G (m1 + m2)

    :param m1: mass 1
    :type m1: float
    :param m2: mass 2
    :type m2: float
    :return: mu
    :rtype: float
    """
    from astropy.constants import G

    return G * (m1 + m2)

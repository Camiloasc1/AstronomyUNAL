import numpy as np

from astropy import units as u

from CelestialMechanics.mu import mu_sun, mu_gm1m2


def T(a: float, m1: float, m2: float) -> float:
    """
    T = sqrt((4 * pi^2) / G(m1 + m2) * a^3)

    :param a: semi-major axis
    :type a: float
    :param m1: mass 1
    :type m1: float
    :param m2: mass 2
    :type m2: float
    :return: translation period
    :rtype: float
    """
    T = (4 * np.pi * np.pi) / (mu_gm1m2(m1, m2)) * a * a * a
    T = np.sqrt(T)

    T = float(T)
    return T


def T_sun(a: float, m2_over_m1: float) -> float:
    """
    T = 2 * pi * sqrt(a^3) / (K * sqrt(1 + m2/m1))

    :param a: semi-major axis
    :type a: float
    :param m2_over_m1: m2 / m1
    :type m2_over_m1: float
    :return: translation period
    :rtype: float
    """

    return 2 * np.pi * np.sqrt(a * a * a / mu_sun(m2_over_m1))


def a(T: float, m1: float, m2: float) -> float:
    """
    a = (T^2 * K * (m1 + m2) / (4 * pi^2))^1/3

    :param T: translation period
    :type T: float
    :param m1: mass 1
    :type m1: float
    :param m2: mass 2
    :type m2: float
    :return: semi-major axis
    :rtype: float
    """
    a = T * T * mu_gm1m2(m1, m2) / (4 * np.pi * np.pi)
    a = np.power(a, 1. / 3.)
    return a


def a_sun(T: float, m2_over_m1: float) -> float:
    """
    a = (T * K * (1 + m2/m1) / (4 * pi^2))^2/3

    :param T: translation period
    :type T: float
    :param m2_over_m1: m2 / m1
    :type m2_over_m1: float
    :return: semi-major axis
    :rtype: float
    """
    a = T * T * mu_sun(m2_over_m1) / (4 * np.pi * np.pi)
    a = np.power(a, 1. / 3.)
    return a


def n(T: float) -> float:
    """
    T = 360° / T

    :param T: translation period
    :type T: float
    :return: mean motion in degrees
    :rtype: float
    """
    return 360 / T


def n_sun(a: float, m2_over_m1: float) -> float:
    """
    n = sqrt(mu / a^3)

    :param a: semi-major axis
    :type a: float
    :param m2_over_m1: m2 / m1
    :type m2_over_m1: float
    :return: mean motion in degrees
    :rtype: float
    """
    n = mu_sun(m2_over_m1) / a / a / a
    n = np.sqrt(n) * u.rad
    # n = np.rad2deg(n * u.d) / u.d
    n = n.to(u.deg / u.d)
    return n

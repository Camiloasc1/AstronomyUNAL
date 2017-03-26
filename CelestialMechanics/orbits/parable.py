from typing import Tuple

import numpy as np
from astropy import units as u


def r(q: float, angle: float) -> float:
    """
    r = 2 * q / (1 + np.cos(angle))

    :param q: pericentric distance
    :type q: float
    :param angle: theta angle
    :type angle: float
    :return: radius vector
    :rtype: float
    """
    angle = np.deg2rad(angle)
    return 2 * q / (1 + np.cos(angle))


def angle(q: float, r: float) -> Tuple[float, float]:
    """
    cos(angle) = 2 * q / r - 1
    angle = arccos(2 * q / r - 1)

    :param q: pericentric distance
    :type q: float
    :param r: radius vector
    :type r: float
    :return: angle
    :rtype: float
    """
    angle = np.rad2deg(np.arccos(2 * q / r - 1))
    return angle, (360. * u.deg) - angle


def e() -> float:
    """
    e = 1

    :return: eccentricity
    :rtype: float
    """
    return 1.


def c(a: float) -> float:
    """
    c = -

    :param a: semi-major axis
    :type a: float
    :return: linear eccentricity
    :rtype: float
    """
    return None


def l(a: float) -> float:
    """
    l = 2 * a

    :param a: semi-major axis
    :type a: float
    :return: semi-latus rectum
    :rtype: float
    """
    return 2. * a


def p(a: float) -> float:
    """
    p = 2 * a

    :param a: semi-major axis
    :type a: float
    :return: focal parameter
    :rtype: float
    """
    return 2. * a


def E() -> float:
    """
    E = 0

    :return: energy of the 2-body system
    :rtype: float
    """

    return 0


def v(r: float, m1: float, m2: float) -> float:
    """
    v = sqrt(2 * G * (m1 + m2) / r)

    :param r: radius vector
    :type r: float
    :param m1: mass 1
    :type m1: float
    :param m2: mass 2
    :type m2: float
    :return: velocity
    :rtype: float
    """
    from CelestialMechanics.mu import mu_gm1m2

    v = 2 * mu_gm1m2(m1, m2) / r
    v = np.sqrt(v)

    return v


def solve_C(q: float, mu: float, t_r: float, t: float) -> float:
    """
    C = sqrt(mu / 2 / a ** 3) * (t - t_r)

    :param q: pericentric distance
    :type q: float
    :param mu: G * (m1 + m2)
    :type mu: float
    :param t_r: reference time
    :type t_r: float
    :param t: time
    :type t: float
    :return: C
    :rtype: float
    """
    return np.sqrt(mu / 2 / q ** 3) * (t - t_r)


def solve_S(C: float) -> float:
    """
    S = arctan(2 / 3 / C)

    :param C: C
    :type C: float
    :return: S
    :rtype: float
    """
    return np.arctan(2. / 3. / C)


def solve_FI(S: float) -> float:
    """
    FI = arctan(tan(S / 2) ** (1. / 3.))

    :param S: S
    :type S: float
    :return: FI
    :rtype: float
    """
    S = S.to(u.deg)
    print(S)
    angle = np.tan(S / 2).value
    print(type(angle))
    print(angle)
    print(angle ** (1. / 3.))  ### Why!!!!
    return np.arctan(angle ** (1. / 3.))


def angle_FI(FI: float) -> float:
    """
    theta = 2 * np.arctan(2 / np.tan(2 * FI))

    :param FI: FI
    :type FI: float
    :return: theta angle
    :rtype: float
    """
    return 2 * np.arctan(2 / np.tan(2 * FI))


def r1(q: float, angle: float, mu: float) -> float:
    """
    r. = sqrt(mu / 2 / q) * e * sin(angle)

    :param q: pericentric distance
    :type q: float
    :param angle: theta angle
    :type angle: float
    :param mu: G * (m1 + m2)
    :type mu: float
    :return: r.
    :rtype: float
    """
    r1 = mu / 2 / q
    r1 = np.sqrt(r1)
    r1 *= np.sin(angle)
    return r1


def r_angle1(q: float, r: float, mu: float) -> float:
    """
    rtheta.= sqrt(2 * q * mu) / r

    :param q: pericentric distance
    :type q: float
    :param r: radius vector
    :type r: float
    :param mu: G * (m1 + m2)
    :type mu: float
    :return: r.
    :rtype: float
    """
    r_angle1 = 2 * q * mu
    r_angle1 = np.sqrt(r_angle1)
    r_angle1 /= r
    return r_angle1

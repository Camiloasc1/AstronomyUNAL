import numpy as np
from typing import Tuple

from astropy import units as u


def r(a: float, e: float, angle: float) -> float:
    """
    r = (a * (1 - e * e)) / (1 + e * cos(angle))

    :param a: semi-major axis
    :type a: float
    :param e: eccentricity
    :type e: float
    :param angle: theta angle
    :type angle: float
    :return: radius vector
    :rtype: float
    """
    angle = np.deg2rad(angle)

    r = (a * (1. - e * e)) / (1. + e * np.cos(angle))

    r = float(r)
    return r


def angles(a: float, e: float, r: float) -> Tuple[float, float]:
    """
    theta = arccos((a * (1 - e * e) - r) / (e * r))

    :param a: semi-major axis
    :type a: float
    :param e: eccentricity
    :type e: float
    :param r: radius vector
    :type r: float
    :return: the two theta angles
    :rtype: (float, float)
    """
    angle = a * (1 - e * e) - r
    angle /= e * r
    angle = np.arccos(angle)
    angle = angle.to(u.deg)
    return angle, (360. * u.deg) - angle


def ae(q: float, Q: float) -> Tuple[float, float]:
    """
    a = (Q + q) / 2
    e = (Q - q) / (2 * a)

    :param q: pericentric distance
    :type q: float
    :param Q: apocentric distance
    :type Q: float
    :return: semi-major axis, eccentricity
    :rtype: (float,float)
    """
    a = (Q + q) / 2
    return a, (Q - q) / (2 * a)


def ae_v_sun(vq: float, vQ: float, m2_over_m1: float) -> Tuple[float, float]:
    """
    a = (Q + q) / 2
    e = (Q - q) / (2 * a)

    :param vq: pericentric velocity
    :type vq: float
    :param vQ: apocentric velocity
    :type vQ: float
    :param m2_over_m1: m2 / m1
    :type m2_over_m1: float
    :return: semi-major axis, eccentricity
    :rtype: (float,float)
    """
    from CelestialMechanics.mu import mu_sun

    v = vq / vQ
    e = (v - 1) / (v + 1)
    a = mu_sun(m2_over_m1) / vq ** 2 * (1 + e) / (1 - e)
    a = mu_sun(m2_over_m1) / vQ ** 2 * (1 - e) / (1 + e)
    return a, e


def b(a: float, e: float) -> float:
    """
    b = a * a * (1 - e * e)

    :param a: semi-major axis
    :type a: float
    :param e: eccentricity
    :return: semi-minor axis
    :rtype: float
    """
    return a * a * (1 - e * e)


def qQ(a: float, e: float) -> Tuple[float, float]:
    """
    q = a * (1 - e)
    Q = a * (1 + e)

    :param a: semi-major axis
    :type a: float
    :param e: eccentricity
    :type e: float
    :return: pericentric distance, apocentric distance
    :rtype: (float, float)
    """
    return q(a, e), Q(a, e)


def q(a: float, e: float) -> float:
    """
    q = a * (1 - e)

    :param a: semi-major axis
    :type a: float
    :param e: eccentricity
    :type e: float
    :return: pericentric distance
    :rtype: float
    """
    return a * (1 - e)


def Q(a: float, e: float) -> float:
    """
    Q = a * (1 + e)

    :param a: semi-major axis
    :type a: float
    :param e: eccentricity
    :type e: float
    :return: apocentric distance
    :rtype: float
    """
    return a * (1 + e)


def e(a: float, b: float) -> float:
    """
    e = sqrt(1 - (b * b) / (a * a))

    :param a: semi-major axis
    :type a: float
    :param b: semi-minor axis
    :type b: float
    :return: eccentricity
    :rtype: float
    """
    return np.sqrt(1 - (b * b) / (a * a))


def c(a: float, b: float) -> float:
    """
    c = sqrt(a * a - b * b)

    :param a: semi-major axis
    :type a: float
    :param b: semi-minor axis
    :type b: float
    :return: linear eccentricity
    :rtype: float
    """
    return np.sqrt(a * a - b * b)


def l(a: float, b: float) -> float:
    """
    l = b * b / a

    :param a: semi-major axis
    :type a: float
    :param b: semi-minor axis
    :type b: float
    :return: semi-latus rectum
    :rtype: float
    """
    return b * b / a


def p(a: float, b: float) -> float:
    """
    p = b * b / c(a, b)

    :param a: semi-major axis
    :type a: float
    :param b: semi-minor axis
    :type b: float
    :return: focal parameter
    :rtype: float
    """
    return b * b / c(a, b)


def E(a: float, m1: float, m2: float) -> float:
    """
    E = -G * m1 * m2 / (2 * a)

    :param a: semi-major axis
    :type a: float
    :param m1: mass 1
    :type m1: float
    :param m2: mass 2
    :type m2: float
    :return: energy of the 2-body system
    :rtype: float
    """
    from astropy.constants import G

    return -G * m1 * m2 / (2 * a)


def v(r: float, a: float, m1: float, m2: float) -> float:
    """
    v = sqrt(mu  * (2 / r - 1 / a))

    :param r: radius vector
    :type r: float
    :param a: semi-major axis
    :type a: float
    :param m1: mass 1
    :type m1: float
    :param m2: mass 2
    :type m2: float
    :return: velocity
    :rtype: float
    """
    from CelestialMechanics.mu import mu_gm1m2

    v = mu_gm1m2(m1, m2) * (2 / r - 1 / a)
    v = np.sqrt(v)

    return v


def v_sun(r: float, a: float, m2_over_m1: float) -> float:
    """
    v = sqrt(mu * (2 / r - 1 / a))

    :param r: radius vector
    :type r: float
    :param a: semi-major axis
    :type a: float
    :param m2_over_m1: m2 / m1
    :type m2_over_m1: float
    :return: velocity
    :rtype: float
    """
    from CelestialMechanics.mu import mu_sun

    v = mu_sun(m2_over_m1) * (2 / r - 1 / a)
    v = np.sqrt(v)

    return v


def vqQ(a: float, e: float, m1: float, m2: float) -> Tuple[float, float]:
    """

    :param a: semi-major axis
    :type a: float
    :param e: eccentricity
    :type e: float
    :param m1: mass 1
    :type m1: float
    :param m2: mass 2
    :type m2: float
    :return: velocity in q and Q
    :rtype: (float,float)
    """
    from CelestialMechanics.mu import mu_gm1m2

    vq = mu_gm1m2(m1, m2) / a
    vQ = mu_gm1m2(m1, m2) / a
    vq *= (1 + e) / (1 - e)
    vQ *= (1 - e) / (1 + e)
    vq = np.sqrt(vq)
    vQ = np.sqrt(vQ)
    return vq, vQ


def vqQ_sun(a: float, e: float, m2_over_m1: float) -> Tuple[float, float]:
    """

    :param a: semi-major axis
    :type a: float
    :param e: eccentricity
    :type e: float
    :param m2_over_m1: m2 / m1
    :type m2_over_m1: float
    :return: velocity in q and Q
    :rtype: (float,float)
    """
    from CelestialMechanics.mu import mu_sun

    vq = mu_sun(m2_over_m1) / a
    vQ = mu_sun(m2_over_m1) / a
    vq *= (1 + e) / (1 - e)
    vQ *= (1 - e) / (1 + e)
    vq = np.sqrt(vq)
    vQ = np.sqrt(vQ)
    return vq, vQ


def n(a: float, mu: float) -> float:
    """
    n = np.sqrt(mu / a ** 3)

    :param a: semi-major axis
    :type a: float
    :param mu: G * (m1 + m2)
    :type mu: float
    :return: mean movement
    :rtype: float
    """
    return np.sqrt(mu / a ** 3) * u.rad


def angle_Mr(Lr: float, w_: float) -> float:
    """
    Mr = Lr - w_

    :param Lr: mean longitude at t_r
    :type Lr: float
    :param w_: longitude of the preriapsis
    :type w_: float
    :return: mean anomaly at t_r
    :rtype: float
    """
    return Lr - w_


def angle_M(M_r: float, n: float, t_r: float, t: float) -> float:
    """
    M = Mr + n * (t - tr)

    :param M_r: mean anomaly at tr
    :type M_r: float
    :param n: mean movement
    :type n: float
    :param t_r: reference time
    :type t_r: float
    :param t: time
    :type t: float
    :return: mean anomaly at t
    :rtype: float
    """

    return M_r + n * (t - t_r)


def angle_M_eE(e: float, E: float) -> float:
    """
    M = E - e * np.sin(E)

    :param e: eccentricity
    :type e: float
    :param E: eccentric anomaly
    :type E: float
    :return: mean anomaly
    :rtype: float
    """

    # e = np.rad2deg(e * u.rad)
    e = e * u.rad
    return E - e * np.sin(E)


def solve_E(M: float, e: float, ROUNDS: float) -> float:
    """
    Solve the kepler's equation with:
    E_n = M + e * np.sin(E_n-1)

    :param M: mean anomaly
    :type M: float
    :param e: eccentricity
    :type e: float
    :param ROUNDS: steps
    :type ROUNDS: float
    :return: eccentric anomaly
    :rtype: float
    """

    E = M
    # e = np.rad2deg(e) * u.deg
    e = e * u.rad
    for i in range(ROUNDS):
        E = M + e * np.sin(E)
    return E


def r_E(a: float, e: float, E: float) -> float:
    """
    r = a * (1 - e * cos(E))

    :param a: semi-major axis
    :type a: float
    :param e: eccentricity
    :type e: float
    :param E: eccentric anomaly
    :type E: float
    :return: radius vector
    :rtype: float
    """
    return a * (1 - e * np.cos(E))


def angle_E(e: float, E: float) -> float:
    """
    theta = 2 * arctan(sqrt((1 + e) / (1 - e)) * tan(E / 2))

    :param e: eccentricity
    :type e: float
    :param E: eccentric anomaly
    :type E: float
    :return: theta angle
    :rtype: float
    """
    return 2 * np.arctan(np.sqrt((1 + e) / (1 - e)) * np.tan(E / 2))


def E_angle(e: float, angle: float) -> float:
    """
    E = 2 * arctan(np.sqrt((1 - e) / (1 + e)) * tan(theta / 2))

    :param angle: true anomaly
    :type angle: float
    :param e: eccentricity
    :type e: float
    :return: eccentric anomaly
    :rtype: float
    """
    return 2 * np.arctan(np.sqrt((1 - e) / (1 + e)) * np.tan(angle / 2))


def r1(a: float, e: float, angle: float, mu: float) -> float:
    """
    r. = sqrt(mu / a / (1 - e * e)) * e * sin(angle)

    :param a: semi-major axis
    :type a: float
    :param e: eccentricity
    :type e: float
    :param angle: theta angle
    :type angle: float
    :param mu: G * (m1 + m2)
    :type mu: float
    :return: r.
    :rtype: float
    """
    r1 = mu / a / (1 - e * e)
    r1 = np.sqrt(r1)
    r1 *= e * np.sin(angle)
    return r1


def r_angle1(a: float, e: float, r: float, mu: float) -> float:
    """
    rthetha.= sqrt(mu * a * (1 - e * e)) / r

    :param a: semi-major axis
    :type a: float
    :param e: eccentricity
    :type e: float
    :param r: radius vector
    :type r: float
    :param mu: G * (m1 + m2)
    :type m1: float
    :return: r theta.
    :rtype: float
    """
    r_angle1 = mu * a * (1 - e * e)
    r_angle1 = np.sqrt(r_angle1)
    r_angle1 /= r
    return r_angle1


def delta_t_t0_Mn(M: float, n: float) -> float:
    """

    :param M: mean anomaly
    :type M: float
    :param n: mean movement
    :type n: float
    :return: t - t0
    :rtype: float
    """
    return M / n


def delta_t_t0_aeangle(a: float, e: float, angle: float, mu: float) -> float:
    """

    :param a: semi-major axis
    :type a: float
    :param e: eccentricity
    :type e: float
    :param angle: true anomaly
    :type angle: float
    :param mu: G * (m1 + m2)
    :type mu: float
    :return: t - t0
    :rtype: float
    """
    E = E_angle(e, angle)
    M = angle_M_eE(e, E)
    return delta_t_t0_Mn(M, n(a, mu))


def t0(M_r: float, t_r: float, n: float) -> float:
    """
    t0 = t_r - M_r / n

    :param M_r: mean anomaly at t_r
    :type M_r: float
    :param t_r: reference time
    :type t_r: float
    :param n: mean movement
    :type n: float
    :return: t0
    :rtype: float
    """
    return t_r - M_r / n

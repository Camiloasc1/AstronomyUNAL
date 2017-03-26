from typing import List, Tuple

import numpy as np

from astropy import units as u

ROUNDS = 15


def solve(a: float, e: float, W: float, w: float, i: float, M_r: float, mu: float, t_r: float, t: float) -> Tuple[
    Tuple[float, float, float], Tuple[float, float, float]]:
    """
    Solve the 2-bodies problem with the orbital elements

    :param a: semi-major axis
    :type a: float
    :param e: eccentricity
    :type e: float
    :param W: longitude of the ascending node
    :type W: float
    :param w: argument of periapsis
    :type w: float
    :param i: inclination
    :type i: float
    :param M_r: mean anomaly at tr
    :type M_r: float
    :param mu: G * (m1 + m2)
    :type mu: float
    :param t_r: reference time
    :type t_r: float
    :param t: time
    :type t: float
    :return: position and velocity in (heliocentric ecliptical) rectangular coordinates
    :rtype: (float, float, float), (float, float, float)
    """
    if e < 1.:  # Ellipse
        r, angle, r1, r_angle1 = solve_ellipse(a, e, M_r, mu, t_r, t)
    elif e > 1.:  # Hyperbola
        r, angle, r1, r_angle1 = solve_hyperbola(a, e, M_r, mu, t_r, t)
    else:  # Parable
        r, angle, r1, r_angle1 = solve_parable(a, e, M_r, mu, t_r, t)

    # x = y = z = r
    # x *= np.cos(W) * np.cos(w + angle) - np.sin(W) * np.sin(w + angle) * np.cos(i)
    # y *= np.sin(W) * np.cos(w + angle) - np.cos(W) * np.sin(w + angle) * np.cos(i)
    # z *= np.sin(w + angle) * np.sin(i)

    v = np.array([[r.to(u.au).value * np.cos(angle)],
                  [r.to(u.au).value * np.sin(angle)],
                  [0]])
    x, y, z = np.matmul(S(W, w, i), v).flatten() * u.au

    v1 = np.array([[r1.to(u.au / u.d).value * np.cos(angle) - r_angle1.to(u.au / u.d).value * np.sin(angle)],
                   [r1.to(u.au / u.d).value * np.sin(angle) + r_angle1.to(u.au / u.d).value * np.cos(angle)],
                   [0]])
    x1, y1, z1 = np.matmul(S(W, w, i), v1).flatten() * u.au / u.d

    # print(S(W, w, i))
    # print(v1)
    # print(x1, y1, z1)

    return (x, y, z), (x1, y1, z1)


def solve_ellipse(a: float, e: float, M_r: float, mu: float, t_r: float, t: float) -> Tuple[float, float, float, float]:
    """
    Solve inside the orbital plane for the ellipse case

    :param a: semi-major axis
    :type a: float
    :param e: eccentricity
    :type e: float
    :param M_r: mean anomaly at tr
    :type M_r: float
    :param mu: G * (m1 + m2)
    :type mu: float
    :param t_r: reference time
    :type t_r: float
    :param t: time
    :type t: float
    :return: radius vector, theta angle, r., r theta.
    :rtype: (float, float, float, float)
    """
    import CelestialMechanics.orbits.ellipse as ellipse

    n = ellipse.n(a, mu)
    M = ellipse.angle_M(M_r, t_r, n, t)
    E = ellipse.solve_E(M, e, ROUNDS)
    r = ellipse.r_E(a, e, E)
    angle = ellipse.angle_E(e, E)
    r1 = ellipse.r1(a, e, angle, mu)
    r_angle1 = ellipse.r_angle1(a, e, r, mu)
    return r, angle, r1, r_angle1


def solve_parable(a: float, e: float, M_r: float, mu: float, t_r: float, t: float) -> Tuple[float, float, float, float]:
    """
    Solve inside the orbital plane for the parable case

    :param a: semi-major axis
    :type a: float
    :param e: eccentricity
    :type e: float
    :param M_r: mean anomaly at tr
    :type M_r: float
    :param mu: G * (m1 + m2)
    :type mu: float
    :param t_r: reference time
    :type t_r: float
    :param t: time
    :type t: float
    :return: radius vector, theta angle, r., r theta.
    :rtype: (float, float, float, float)
    """
    pass


def solve_hyperbola(a: float, e: float, M_r: float, mu: float, t_r: float, t: float) -> Tuple[float, float, float, float]:
    """
    Solve inside the orbital plane for the hyperbola case

    :param a: semi-major axis
    :type a: float
    :param e: eccentricity
    :type e: float
    :param M_r: mean anomaly at tr
    :type M_r: float
    :param mu: G * (m1 + m2)
    :type mu: float
    :param t_r: reference time
    :type t_r: float
    :param t: time
    :type t: float
    :return: radius vector, theta angle, r., r theta.
    :rtype: (float, float, float, float)
    """
    pass


def SolveE(a: float, e: float, M: float, mu: float, dt: float) -> float:
    if e < 1.:  # Ellipse
        import CelestialMechanics.orbits.ellipse as ellipse
        return ellipse.solve_E(M, e, ROUNDS)
    elif e > 1.:  # Hyperbola
        Mh = np.sqrt(mu / a ** 3) * dt
        F = np.arcsinh(1 / e * (2 * Mh))
        for i in range(ROUNDS):
            F = np.arcsinh(1 / e * (Mh + F))
        return F
    else:  # Parable
        C = np.sqrt(mu / 2 / a ** 3) * dt
        return C


def S(W: float, w: float, i: float) -> List[List[float]]:
    """
    The rotation matrix S

    :param W: longitude of the ascending node
    :type W: float
    :param w: argument of periapsis
    :type w: float
    :param i: inclination
    :type i: float
    :return: rotation matrix S
    :rtype: [[float, float, float], [float, float, float], [float, float, float]]
    """

    s11 = np.cos(W) * np.cos(w) - np.sin(W) * np.sin(w) * np.cos(i)
    s12 = -np.cos(W) * np.sin(w) - np.sin(W) * np.cos(w) * np.cos(i)
    s13 = np.sin(W) * np.sin(i)
    s21 = np.sin(W) * np.cos(w) + np.cos(W) * np.sin(w) * np.cos(i)
    s22 = -np.sin(W) * np.sin(w) + np.cos(W) * np.cos(w) * np.cos(i)
    s23 = -np.cos(W) * np.sin(i)
    s31 = np.sin(i) * np.sin(w)
    s32 = np.sin(i) * np.cos(w)
    s33 = np.cos(i)
    return [[s11, s12, s13],
            [s21, s22, s23],
            [s31, s32, s33]]

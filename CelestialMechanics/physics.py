from astropy import units as u


def GravitationalForce(r, m1, m2):
    from astropy.constants import G
    from numpy.linalg import norm

    r = r.to(u.m)
    r = r.value / norm(r) ** 3 / u.m ** 2
    return G * m1 * m2 * r

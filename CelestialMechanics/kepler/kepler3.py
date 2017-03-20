import numpy as np
from astropy import constants as astroconst
from CelestialMechanics.kepler.constants import *


def T(a, m1, m2):
    T = (4 * np.pi * np.pi) / (mu_gm1m2(m1, m2)) * a * a * a
    T = np.sqrt(T)

    T = float(T)
    return T


def T_sun(a, m2_over_m1):
    T = 2 * np.pi * np.sqrt(a * a * a)
    T = T / mu_sun(m2_over_m1)

    T = float(T)
    return T


def a(T, m1, m2):
    a = T * T * mu_gm1m2(m1, m2) / (4 * np.pi * np.pi)
    a = np.power(a, 1. / 3.)

    a = float(a)
    return a


def a_sun(T, m1_over_m2):
    a = T * T * mu_sun(m1_over_m2) * mu_sun(m1_over_m2) / (4 * np.pi * np.pi)
    a = np.power(a, 1. / 3.)

    a = float(a)
    return a


def n(T):
    return 360 / T


def n_sun(a, m1_over_m2):
    n = mu_sun(m1_over_m2)
    n = n / np.sqrt(a * a * a)
    n = np.rad2deg(n)

    n = float(n)
    return n


def mu_sun(m1_over_m2):
    return K * np.sqrt(1. + m1_over_m2)


def mu_na(n, a):
    return n * n * a * a * a


def mu_gm1m2(m1, m2):
    return astroconst.G.value * (m1 + m2)

import numpy as np
from astropy import units as u

from CelestialMechanics import util

def eclip2ecu_(r):
    x, y, z = eclip2ecu(r[0], r[1], r[2])
    return [x.value, y.value, z.value] * x.unit

def eclip2ecu(x, y, z):
    E = 23.43927944 * u.deg
    s11 = 1
    s12 = 0
    s13 = 0
    s21 = 0
    s22 = np.cos(E)
    s23 = -np.sin(E)
    s31 = 0
    s32 = np.sin(E)
    s33 = np.cos(E)
    S = [[s11, s12, s13],
         [s21, s22, s23],
         [s31, s32, s33]] * (u.d / u.d)
    geo = [[x.value],
           [y.value],
           [z.value]] * x.unit
    x, y, z = util.matmul(S, geo).flatten()
    return [x, y, z] * x.unit


def ecu2eclip_(r):
    x, y, z = ecu2eclip(r[0], r[1], r[2])
    return [x.value, y.value, z.value] * x.unit


def ecu2eclip(x, y, z):
    E = 23.43927944 * u.deg
    s11 = 1
    s12 = 0
    s13 = 0
    s21 = 0
    s22 = np.cos(E)
    s23 = -np.sin(E)
    s31 = 0
    s32 = np.sin(E)
    s33 = np.cos(E)
    S = [[s11, s12, s13],
         [s21, s22, s23],
         [s31, s32, s33]] * (u.d / u.d)
    S = np.transpose(S)
    geo = [[x.value],
           [y.value],
           [z.value]] * x.unit
    x, y, z = util.matmul(S, geo)
    return (x, y, z)

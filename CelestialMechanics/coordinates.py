import numpy as np
from astropy import units as u


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
         [s31, s32, s33]]
    geo = np.array([[x.to(u.au).value],
                    [y.to(u.au).value],
                    [z.to(u.au).value]])
    x, y, z = np.matmul(S, geo).flatten() * u.au
    return (x, y, z)


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
         [s31, s32, s33]]
    S = np.transpose(S)
    geo = np.array([[x.to(u.au).value],
                    [y.to(u.au).value],
                    [z.to(u.au).value]])
    x, y, z = np.matmul(S, geo).flatten() * u.au
    return (x, y, z)

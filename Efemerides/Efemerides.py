import math


def JD(day, month, year, UT):
    return 367.0 * year - math.floor(7.0 / 4.0 * (year + math.floor((month + 9.0) / 12.0))) + math.floor(
        (275.0 * month) / 9.0) + day + 1721013.5 + UT / 24.0


def T(JD):
    return (JD - 2451545.0) / 36525.0


def TSMG(T):
    return 24110.54841 + 8640184.812866 * T + 0.093104 * T * T - 6.2E-6 * T * T * T


def TSMGUT(TSMG, UT):
    return TSMG + UT * 1.002737909


def UT(Tof, HH):
    return Tof - HH


def conv(d, m, s, H):
    return (d + m / 60.0 + s / 3600.0) * (1.0 if H == 'W' else -1.0)


jd = JD(29, 8, 1996, 0)
t = T(jd)
tsmg = TSMG(t) / 3600.0 % 24
tsmgut = TSMGUT(tsmg, UT(12, -5)) % 24
delta = conv(74, 5, 0, 'E') / 15.0

print 'Fecha Juliana JD=', jd
print 'Fraccion de Siglo T=', t
print 'TSMG 0h UT=', tsmg, '=', TSMG(t)
print 'TSMG UT=', tsmgut, '=', TSMGUT(tsmg, UT(12, -5))
print 'Longitud =', conv(74, 5, 0, 'E'), ', Correcion por Longitud =', delta
print 'Ascension Recta =', tsmgut + delta

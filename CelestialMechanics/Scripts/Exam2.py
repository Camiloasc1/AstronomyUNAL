from astropy import units as u
from astropy.time import Time
from CelestialMechanics.mu import mu_sun
from CelestialMechanics.orbital_determination.orbital_determination import solve, solve_gauss
from CelestialMechanics import util
import numpy as np
from CelestialMechanics.canonical_units.canonical_units import x1

# 1

r = [0.2864083741928, 2.6239808837979, -0.0040602075566] * u.au
r1 = [-0.0095334851080, 0.0031642149156, 0.0014689459610] * u.au / u.d
t_r = Time('2014-05-23T00:00:00', format='isot', scale='tt').jd * u.d
mu = mu_sun(0.)

a, e, i, W, w, M_r, t0 = solve(r, r1, mu, t_r)
print('a', a)
print('e', e)
print('i', i.to(u.deg))
print('W', W.to(u.deg))
print('w', w.to(u.deg))
print(Time(t0, format='jd', scale='utc').isot)

# 2
print()

r_1 = [1.345274, 3.454566, 4.467456] * u.au
r_2 = [2.567474, 2.565789, 1.025774] * u.au

h_ = util.cross(r_1, r_2)
h_ /= util.norm(h_).value
h = util.norm(h_)

i = np.arccos(h_[2] / h)
W = np.arctan2(h_[0], -h_[1])
print('i', i.to(u.deg))
print('W', W.to(u.deg).value % 360)

# 3
print()

m1 = 0.99904612
m2 = 0.000953879
x = x1(m2) * 5.202 * u.au
print(x)
print(x.to(u.km))

# 4
print()

rt = 6378140 * u.m
re = rt + 8848 * u.m
s = 8000 * u.m / u.s
angleImpact = 330 * u.deg
angleEverest = 35 * u.deg
h = 50983105851.22712 * u.m * u.m / u.s

r_1 = [rt.value * np.cos(360 * u.deg - angleImpact), rt.value * np.sin(360 * u.deg - angleImpact), 0] * u.m
r_2 = [rt.value * np.cos(angleImpact), rt.value * np.sin(angleImpact), 0] * u.m
r_3 = [re.value * np.cos(angleEverest), re.value * np.sin(angleEverest), 0] * u.m

t1 = Time('2017-05-30T00:00:00', format='isot', scale='utc').jd * u.d
a, e, i, W, w, M_r, t0 = solve_gauss(r_1, r_2, r_3, mu_sun(0), t1)

print(a)  # not exact as the sample problem
print(e)  # not exact as the sample problem
v = np.arcsin(h / rt / s).to(u.deg)
print(v)

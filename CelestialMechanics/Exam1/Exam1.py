import numpy as np
from astropy import units as u
from astropy.time import Time

from CelestialMechanics.kepler import constants
from CelestialMechanics.kepler import kepler3
from CelestialMechanics.mu import mu_sun, mu_gm1m2
from CelestialMechanics.orbital_elements.orbital_elements import solve, solve_ellipse, solve_hyperbola, solve_parable
from CelestialMechanics.orbits import ellipse, parable, hyperbola

# 1
print("# 1")

a = 2.6785123 * u.au
e = 0.2543422
i = 0 * u.deg
W = 0 * u.deg
w = 0 * u.deg
M_r = 0 * u.deg

T = kepler3.T_sun(a, 0)
print('T', T)
mu = mu_sun(0)
t_r = Time('2006-06-30T00:00:00Z', format='isot', scale='utc').jd * u.d
t = Time('2006-05-05T00:00:00Z', format='isot', scale='utc').jd * u.d
r, angle, r1, r_angle1 = solve_ellipse(a, e, M_r, mu, t_r, t)
print('theta', angle.to(u.deg))
print('theta', angle.to(u.deg) % (360 * u.deg))
print('r', r)
(x, y, z), (x1, y1, z1) = solve(a, e, W, w, i, M_r, mu, t_r, t)
print('r', np.sqrt(x * x + y * y + z * z))

# 2
print("# 2")

R_earth = 6378.14 * u.km
m1 = 5.97E24 * u.kg
m2 = 0 * u.kg

q = R_earth + 300 * u.km
Q = R_earth + 500 * u.km
print('q', q)
print('Q', Q)

a, e = ellipse.ae(q, Q)
print('a', a)
print('e', e)
T = kepler3.T(a, m1, m2)
print('T', T.to(u.h))
n = ellipse.n(a, mu_gm1m2(m1, m2))
print('n', n.to(u.deg / u.d))
vq = ellipse.v(q, a, m1, m2)
print('vq', vq.to(u.km / u.s))
vq, vQ = ellipse.vqQ(a, e, m1, m2)
print('vq', vq.to(u.km / u.s))

# 4
print("# 4")

t0 = Time('2014-05-28T00:00:00Z', format='isot', scale='utc').jd * u.d
print('t0', t0)
t1 = Time('2014-06-21T00:00:00Z', format='isot', scale='utc').jd * u.d
print('t1', t1)
L_r = 245.26107
w = 102.9777
print('M_r', L_r - w)

# 5
print("# 5")
t0 = Time('2014-01-03T00:00:00Z', format='isot', scale='utc').jd * u.d + 0.633 * u.d

a = 1 * u.au
e = 0.01669
r = 1 * u.au
mu = mu_sun(1 / constants.Earth_Moon)

angles = ellipse.angle(a, e, r)
print('theta1', angles[0])
print('theta2', angles[1])

delta_t1_t0 = ellipse.delta_t_t0_aeangle(a, e, angles[0], mu) % (1 * u.yr).to(u.d)  # module 1 year
print('t1', Time(t0 + delta_t1_t0, format='jd', scale='utc').isot)
delta_t2_t0 = ellipse.delta_t_t0_aeangle(a, e, angles[1], mu) % (1 * u.yr).to(u.d)  # module 1 year
print('t2', Time(t0 + delta_t2_t0, format='jd', scale='utc').isot)

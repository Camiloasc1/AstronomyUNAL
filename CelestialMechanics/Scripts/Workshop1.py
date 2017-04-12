import numpy as np

from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord

from CelestialMechanics.coordinates import ecu2eclip
from CelestialMechanics.mu import mu_sun
from CelestialMechanics.orbital_determination.orbital_determination import solve_gauss


def SkyCoord2xyz(c):
    return [c.cartesian.x.to(u.au).value, c.cartesian.y.to(u.au).value, c.cartesian.z.to(u.au).value] * u.au


c_1 = SkyCoord('02h54m23.18s', '+17d03m19.53s', distance=2.236612957 * u.au)
c_2 = SkyCoord('03h08m30.88s', '+18d05m15.88s', distance=2.265868176 * u.au)
c_3 = SkyCoord('03h25m34.90s', '+19d14m01.58s', distance=2.300057744 * u.au)

# To x,y,z
c_1 = ecu2eclip(SkyCoord2xyz(c_1))
c_2 = ecu2eclip(SkyCoord2xyz(c_2))
c_3 = ecu2eclip(SkyCoord2xyz(c_3))

t1 = Time('2017-04-02T00:00:00Z', format='isot', scale='utc').jd * u.d
t2 = Time('2017-04-07T00:00:00Z', format='isot', scale='utc').jd * u.d
t3 = Time('2017-04-13T00:00:00Z', format='isot', scale='utc').jd * u.d

print(t1, t2, t3)

# http://vo.imcce.fr/webservices/miriade/?forms
# Earth 	2017-04-02T00:00:00.00 	-0.9769221914620 	-0.2113582577129 	0.0000161609786
# Earth 	2017-04-07T00:00:00.00 	-0.9565107703851 	-0.2948601443041 	0.0000171051404
# Earth 	2017-04-13T00:00:00.00 	-0.9227685743776 	-0.3921175067709 	0.0000179127328

earth_1 = [-0.9769221914620, -0.2113582577129, 0.0000161609786] * u.au
earth_2 = [-0.9565107703851, -0.2948601443041, 0.0000171051404] * u.au
earth_3 = [-0.9227685743776, -0.3921175067709, 0.0000179127328] * u.au

r_1 = earth_1 + c_1
r_2 = earth_2 + c_2
r_3 = earth_3 + c_3

a, e, i, W, w, M_r, t0 = solve_gauss(r_1, r_2, r_3, mu_sun(0), t1)
print('a', a)
print('e', e)
print('i', i)
print('W', W)
print('w', w)
print('M_r', M_r)
print('t0', t0)

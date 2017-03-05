import  math

from astropy import units as u
from astropy.table.table import Table
from matplotlib import pylab, pyplot
import astropy
import matplotlib


# listStar = Table.read('test.dat', format = 'ascii')
# # print  listStar
# print listStar.keys()
#
# NGC3227 = Table.read('spectra.dat', format = 'ascii')
# # print NGC3227
# pylab.plot(NGC3227['col1'], NGC3227['col2'], c = 'k')
# pylab.show()
# print NGC3227.keys()
#
# pyplot.plot(NGC3227['col1'], NGC3227['col2'], c = 'k')
# pyplot.show()
# afgranadosc@unal.edu.co

def calcAbsMagnitude(apparentMagnitude, distance):
    return apparentMagnitude - 5 * (math.log10(distance)) + 5

pi = 7.282E-6 * u.radian
d = pi.to(u.arcsec).to(u.parsec, equivalencies = u.parallax())
print pi
print pi.to(u.arcsec)
print d
print d.to(u.au)
print d.to(u.km)
ma = calcAbsMagnitude(11.05, d.value)
print ma

print

pi = 1.745E-7 * u.radian
d = pi.to(u.arcsec).to(u.parsec, equivalencies = u.parallax())
print pi
print pi.to(u.arcsec)
print d
print d.to(u.au)
print d.to(u.km)
mb = calcAbsMagnitude(-0.7, d.value)
print mb

print

print ma - mb
print 100.0 ** ((ma - mb) / 5.0)

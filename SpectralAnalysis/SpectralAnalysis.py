import math

from astropy import units as u
from astropy.modeling import models, fitting
from astropy.table.table import Table
from matplotlib import pyplot
from scipy.integrate import simps
import numpy


class ClimbingAgent(object):
    def __init__(self, world, model, position, movementRange, minimumHeight, hillUp = True):
        self.world = world
        self.model = model
        self.position = position
        self.movementRange = movementRange
        self.minimumHeight = minimumHeight
        self.hillUp = hillUp

    def __str__(self):
        return "Agent [at: " + str(self.position) + " value: " + str(self.value) + "]"

    def climb(self):
        for i in xrange(self.position - self.movementRange, self.position + self.movementRange + 1):
            if i >= len(self.world):
                continue
            if (self.hillUp and self.world[i][1] > self.value[1]) or (not self.hillUp and self.world[i][1] < self.value[1]):
                self.position = i

    @property
    def value(self):
        return self.world[self.position]

    @property
    def isAlive(self):
        maxVal, minVal = self.position, self.position
        currentValue = self.value[1]
        spectedValue = self.model(self.value[0])

        for i in xrange(self.position - self.movementRange, self.position + self.movementRange + 1):
            if i >= len(self.world):
                continue
            if self.world[i][1] > self.world[maxVal][1]:
                maxVal = i
            if self.world[i][1] < self.world[minVal][1]:
                minVal = i

        if self.hillUp :
            if self.world[maxVal][1] > currentValue or currentValue < spectedValue:
                return False
        else:
            if self.world[minVal][1] < currentValue or currentValue > spectedValue:
                return False
        if abs(currentValue - spectedValue) < self.minimumHeight:
            return False
        return True

def read(name):
    dataTable = Table.read(name + '.xml', format = 'votable')
    freq = numpy.array(dataTable['Frequency']) * u.Hz
    wavelength = freq.to(u.AA, equivalencies = u.spectral())
    flux = numpy.array(dataTable['Flux']) * u.Watt / (u.m * u.m * u.Hz)
#     dataTable.add_column(Column(wavelength, name = 'Wave Length'), index = 1)
    return freq, wavelength, flux


def fitModel(maxWidth, wavelength, valuesMap, continuous, agent, l):
    # Select only the values in +-maxWidth from l
    x = filter(lambda agent:(l - maxWidth) < agent and agent < (l + maxWidth), wavelength.value)
    y = map(lambda agent:valuesMap[agent], x)
    model = models.Gaussian1D(amplitude = agent.value[1] - continuous(agent.value[0]), mean = l, stddev = 1.)
    # model = models.Lorentz1D(amplitude = agent.value[1] - continuous(agent.value[0]), x_0 = l, fwhm = 1.)
    fitter = fitting.LevMarLSQFitter()
    function = fitter(model, x, y - continuous(x))
    # Plot in green the fitted model
    x = numpy.linspace(l - maxWidth, l + maxWidth, 1000)
    pyplot.plot(x, function(x) + continuous(x), 'g-', label = str(l))
    pyplot.fill_between(x, continuous(x), function(x) + continuous(x), color = 'g', alpha = .5)
    fwhm = function.parameters[2] * (2 * math.sqrt(2 * math.log(2)))
    flux = simps(function(x), x)
    ew = simps(function(x) / continuous(x), x)
    print 'Measured Wavelength:', agent.value[0], '- Measured Flux:', agent.value[1], '- Well Known Wavelength:', l, '- FWHM:', fwhm, '- Integrated Flux:', flux, '- Equivalent Width:', ew
    return fwhm, flux, ew

def Analyze(name, delta, minimumHeight, maxWidth):
    # Read file
    _, wavelength, flux = read(name)
#     wavelength = wavelength / (z + 1) # This not work properly
    values = zip(wavelength.value, flux.value)
    valuesMap = dict(values)
#     WKL = Table.read('WellKnownLambdas.dat', format = 'ascii')

    # Continuous (by linear least squares)
    model = models.Linear1D(1, 0)
    fitter = fitting.LinearLSQFitter()
    continuous = fitter(model, wavelength.value, flux.value)

    # Agents
    # Max (hill up)
    agents = [ClimbingAgent(values, continuous, i, delta, minimumHeight, True) for i in xrange(delta, len(values), delta * 2)]
    # Min (hill down)
    # agents = agents + [ClimbingAgent(values, continuous, i, delta, minimumHeight, False) for i in xrange(delta, len(values), delta * 2)]
    for a in agents:
        a.climb()
    agents = filter(lambda a: a.isAlive, agents)

    #
    print 'Points found for ', name
    found = {}
    for a in agents:
        if a.value[0] not in found:
            found[a.value[0]] = fitModel(maxWidth, wavelength, valuesMap, continuous, a, a.value[0])

    # Plot
    pyplot.title('Spectrum for ' + name)
    pyplot.ylabel('Flux')
    pyplot.xlabel('Wave Length')

    pyplot.plot(wavelength, flux, 'b-', label = 'Spectrum')
    pyplot.plot(wavelength, continuous(wavelength.value), 'k-', label = 'Continuous')
    pyplot.plot([a.value[0] for a in agents], [a.value[1] for a in agents], 'r*', label = 'Interesting Points')
    pyplot.show()

# Taken by Palomar 200in
Analyze('NGC4321', 1, 0.25E-29, 10)
Analyze('NGC5195', 1, 0.15E-28, 10)
Analyze('NGC4579', 1, 0.25E-28, 10)
Analyze('NGC4303', 1, 0.25E-28, 10)
Analyze('NGC4826', 1, 0.25E-28, 10)
# Assigned by the professor
Analyze('NGC931', 3, 1E-29, 10)
Analyze('NGC5929', 2, 0.25E-28, 10)


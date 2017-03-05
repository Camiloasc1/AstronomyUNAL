from astropy import units as u
from astropy.table.table import Table
from matplotlib import pyplot
import numpy

rarange = (12 * u.hourangle, 12.5 * u.hourangle)
nonNullColumns = ['H5', 'H11', 'H32', 'H34', 'H37', 'H40']

def readStars():
    try:
        head = ['H' + str(i) for i in xrange(78)]
        return Table.read('stars.db', names = head, format = 'ascii')
    except IOError:
        print "Range to search in Hipparcos:"
        for i in rarange:
            print '\t', i, '=', i.to(u.deg)
    return None

def filterStars(stars):
    invalid = []
    for row, star in enumerate(stars):
        for c in nonNullColumns:
            if type(star[c]) != numpy.float64:
                invalid.append(row)
                break
        else:
            if star['H11'] < 0.0:
                invalid.append(row)
    stars.remove_rows(invalid)
    return stars

def calcAbsMagnitude(apparentMagnitude, distance):
    return apparentMagnitude - 5 * (numpy.log10(distance)) + 5

def plot(colourIndex, apparentMagnitude, distance, label1, label2, plotColour):
    pyplot.title('HR Diagram for stars from ' + str(rarange[0].value) + 'h to ' + str(rarange[1].value) + 'h')
    pyplot.ylabel('Absolute Magnitude ' + label1)
    pyplot.xlabel(label2)
    pyplot.gca().invert_yaxis()
    pyplot.scatter(colourIndex, calcAbsMagnitude(apparentMagnitude, distance.value) , c = plotColour, marker = '.')
    pyplot.show()

stars = filterStars(readStars())

# appMag = numpy.array(stars['H5'])
b = numpy.array(stars['H32'])
v = numpy.array(stars['H34'])
bv = numpy.array(stars['H37'])
vi = numpy.array(stars['H40'])
parallax = numpy.array(stars['H11']) * u.milliarcsecond
distance = parallax.to(u.parsec, equivalencies = u.parallax())

plot(bv, b, distance, 'B', 'B - V', 'b')
plot(vi, b, distance, 'B', 'V - I', 'b')
plot(bv, v, distance, 'V', 'B - V', 'r')
plot(vi, v, distance, 'V', 'V - I', 'r')


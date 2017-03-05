from matplotlib import pyplot
import numpy as np

def f(a, b, x):
    return a * np.power(x, b)

x = np.array([28., 30., 33., 35., 38.])
y = np.array([-2410., -3033., -3895., -4491., -5717.])

# x = np.log10(x)
# y = np.log10(-y)

b = 2
a = np.sum(np.power(x, b) * y) / np.sum(np.power(x, 2 * b))
print a
# a = (np.sum(y) - np.sum(x) * b) / 5.
# print a
# a = (np.sum(x * y) - np.sum(x * x) * b) / np.sum(x)
# print a
# a = -np.power(10, a)
# x = np.power(10, x)
# y = -np.power(10, y)
x2 = np.array([28., 30., 33., 35., 38., 37., 27.])
print a, b, x, y
print f(a, b, x2)

pyplot.plot(x, y, 'b-')
pyplot.plot(x2, f(a, b, x2), 'r.')
pyplot.show()

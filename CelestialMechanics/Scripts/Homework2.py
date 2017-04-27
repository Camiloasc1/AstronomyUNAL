import matplotlib.pyplot as plt
from scipy.optimize import brentq
import numpy as np

DELTA = 0.01

m2 = 1 / 82

f = lambda x: x - (1 - m2) * (x + m2) / (abs(x + m2) ** 3) - m2 * (x - 1 + m2) / (abs(x - 1 + m2) ** 3)

segments = [[-2, -m2 - DELTA],
            [-m2 + DELTA, 1 - m2 - DELTA],
            [1 - m2 + DELTA, 2]]

roots = [brentq(f, segments[i][0], segments[i][1]) for i in range(len(segments))]

x = [np.linspace(segments[i][0], segments[i][1], round(1 / DELTA)) for i in range(len(segments))]
y = [f(x[i]) for i in range(len(segments))]

plt.plot(roots, [0, 0, 0], marker="o", linestyle='')
for i in range(len(segments)):
    plt.plot(x[i], y[i])

plt.axhline(color="grey", ls="--", zorder=-1)
plt.axvline(color="grey", ls="--", zorder=-1)

plt.axvline(color="red", ls=":", zorder=-1, x=-m2)
plt.axvline(color="red", ls=":", zorder=-1, x=1 - m2)

plt.axis([-2, 2, -10, 10])
plt.show()

print(roots)

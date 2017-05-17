from typing import List

from scipy.optimize import brentq
import numpy as np


def m1(m2: float) -> float:
    return 1 - m2


def x1(m2: float) -> float:
    return -m2


def x2(m2: float) -> float:
    return 1 - m2


def lagrange_points(m2: float) -> List[float]:
    """
    The 3 colinear Lagrange points
    
    :param m2: mass 2
    :type m2: float
    :return: L1, L2, L3
    :rtype: list
    """
    DELTA = 0.01
    f = lambda x: x - (1 - m2) * (x + m2) / (abs(x + m2) ** 3) - m2 * (x - 1 + m2) / (abs(x - 1 + m2) ** 3)

    segments = [[-2, x1(m2) - DELTA],
                [x1(m2) + DELTA, x2(m2) - DELTA],
                [x2(m2) + DELTA, 2]]

    return [brentq(f, segments[i][0], segments[i][1]) for i in range(len(segments))]


def T(T):
    return T / 2 / np.pi

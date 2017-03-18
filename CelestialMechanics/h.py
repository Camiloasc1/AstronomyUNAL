import numpy as np
from numpy.linalg import norm
from typing import List, Tuple


def h_and_angle(r: List[float], r1: List[float]) -> Tuple[List[float], float, float, float]:
    """
    Calculate h= r x r' (cross product)

    :param r: the position
    :type r: list[float]
    :param r1: the velocity
    :type r1: list[float]
    :return: the h vector, the angle, the supplementary angle, the angle in the right quadrant based on norm(r1)
    :rtype: (list[float], float, float, float)
    """

    # The h constant vector
    h = np.cross(r, r1).tolist()

    # Angle between r and r'
    angle = norm(h) / (norm(r) * norm(r1))
    angle = np.arcsin(angle)
    angle = np.rad2deg(angle)

    h = list(h)
    angle = float(angle)
    return h, angle, 180. - angle, angle if norm(r1) < 0 else 180. - angle

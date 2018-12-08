# Defines weight matrix
# data =  points number, point = analysed point, r = minimum radius
import numpy as np
from src.helpers import weight as w


def W(data, point, r, derivative=None):
    weight = []
    if not derivative:
        for index, row in enumerate(data):
            d2d = row[0:2]
            weight.append(w.gaussian_with_radius(np.subtract(d2d, point), r))
        return np.diag(weight)
    else:
        for index, row in enumerate(data):
            d2d = row[0:2]
            weight.append(w.gaussian_with_radius(np.subtract(d2d, point), r, derivative))
        return np.diag(weight)

# Defines weight matrix
# data =  points number, point = analysed point, r = minimum radius
import numpy as np
from src.helpers import weight as w


def W(data, point, r, derivative=None):
    W = []
    if not derivative:
        for index, row in enumerate(data):
            d2d = row[0:2]  # (x,y)
            leftZeroes = np.zeros([1, index])  # 0 0 0 ... 0 } index times
            rightZeroes = np.zeros([1, len(data) - index - 1])
            weightfunction = w.gaussian_with_radius(np.subtract(d2d, point), r)
            newRow = np.concatenate([leftZeroes, [[weightfunction]], rightZeroes], axis=1)[0]
            W.append(newRow)
        return W
    else:
        for index, row in enumerate(data):
            d2d = row[0:2]  # (x,y)
            leftZeroes = np.zeros([1, index])  # 0 0 0 ... 0 } index times
            rightZeroes = np.zeros([1, len(data) - index - 1])
            weight = w.gaussian_with_radius(np.subtract(d2d, point), r, derivative)
            newRow = np.concatenate([leftZeroes, [[weight]], rightZeroes], axis=1)[0]
            W.append(newRow)
        return W

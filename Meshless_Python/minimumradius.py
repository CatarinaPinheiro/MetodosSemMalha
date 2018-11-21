import numpy as np
from numpy import linalg as la


def get_radius(data, point, m):
    distances = []
    for dat in data:
        dif = np.subtract(point, dat[0:2])
        dist = la.norm(dif)
        distances.append(dist)
    distances = sorted(distances)

    # Distances from the coordinates of the analyzed point to the contour points

    return distances[m]


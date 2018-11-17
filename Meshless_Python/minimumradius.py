import numpy as np
from numpy import linalg as la


def get_radius(data, point, m, contour_point=None):
    distances = []
    for dat in data:
        dif = np.subtract(point, dat[0:2])
        dist = la.norm(dif)
        distances.append(dist)
    distances = sorted(distances)

    # Distances from the coordinates of the analyzed point to the contour points
    point_dist = list(map(lambda x: (la.norm(np.subtract(point, x))), contour_point))
    if all(x > distances[m] for x in point_dist):
        return distances[m]
    else:
        return distances[int(1.5*(m+1))]


def minimum_radius(data, point, m):
    distances = []
    for dat in data:
        dif = np.subtract(point, dat[0:2])
        dist = la.norm(dif)
        distances.append(dist)
    distances = sorted(distances)
    return distances[m]
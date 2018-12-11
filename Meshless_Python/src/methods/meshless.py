from scipy.spatial import Delaunay
import numpy as np
from numpy import linalg as la
from src.helpers import unique_rows


def get_radius(self, point, m):
    distances = []
    for dat in self.data:
        dif = np.subtract(point, dat[0:2])
        dist = la.norm(dif)
        distances.append(dist)
    distances = sorted(distances)

    return distances[m]


class MeshlessMethod:
    def __init__(self, data):
        self.data = data

    @property
    def boundary_data(self):
        boundary_data_initial = []
        data_array = np.array(self.data)
        x = data_array[:, 0]
        y = data_array[:, 1]
        x = x.flatten()
        y = y.flatten()
        points2D = np.vstack([x, y]).T
        tri = Delaunay(points2D)
        boundary = (points2D[tri.convex_hull]).flatten()
        bx = boundary[0:-2:2]
        by = boundary[1:-1:2]

        for i in range(len(bx)):
            boundary_data_initial.append([bx[i], by[i]])

        boundary_data = unique_rows(boundary_data_initial)

        return boundary_data

    @property
    def domain_data(self):
        boundary_list = [[x, y] for x, y in self.boundary_data]
        return [x for x in self.data if x not in boundary_list]

########################################


def boundary_function(point):
    if point[0] == 0:
        return 0
    elif point[0] == 10:
        return 100
    elif point[1] == 0:
        return 0
    elif point[1] == 10:
        return 0
    else:
        raise ValueError("point not in boundary")


def domain_function():  # domain_function(point)
    return 0

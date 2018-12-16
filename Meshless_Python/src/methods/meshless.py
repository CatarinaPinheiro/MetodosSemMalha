from scipy.spatial import Delaunay
import numpy as np
from numpy import linalg as la
from src.helpers import unique_rows
from src.helpers import duration
from src.methods import mls2d as mls


def get_radius(self, point, m):
    distances = []
    for dat in self.data:
        dif = np.subtract(point, dat[0:2])
        dist = la.norm(dif)
        distances.append(dist)
    distances = sorted(distances)

    return distances[m]


class MeshlessMethod:
    def __init__(self, data, basis, domain_function, domain_operator, boundary_operator, boundary_function):
        self.data = data
        self.basis = basis
        self.domain_function = domain_function
        self.domain_operator = domain_operator
        self.boundary_operator = boundary_operator
        self.boundary_operator = boundary_operator
        self.m2d = mls.MovingLeastSquares2D(self.data, self.basis)

    def solve(self):
        lphi = []
        b = []

        for i, d in enumerate(self.data):

            duration.duration.start("%d/%d" % (i, len(self.data)))
            self.m2d.point = d

            if d in self.domain_data:

                def integration_element(integration_point, i):
                    if found:
                        return value[i]
                    else:
                        self.m2d.point = integration_point
                        phi = self.m2d.numeric_phi
                        value = phi.eval(d)[0] * self.integration_weight(d, integration_point, radius)
                        return value[i]

                lphi.append(
                    [self.integration(d, radius, lambda p: integration_element(p, i)) for i in range(len(self.data))])
                b.append(self.integration(d, radius, self.domain_function))

            elif d in self.boundary_data:

                lphi.append(self.boundary_operator(self.m2d.numeric_phi, d).eval(d)[0])

                b.append(self.boundary_function(self.m2d.point))
            duration.duration.step()

        return la.solve(lphi, b)

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

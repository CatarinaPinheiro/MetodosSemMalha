import numpy as np
from numpy import linalg as la


class Meshless:
    def __init__(self):
        pass

    def get_radius(self, point, m):
        distances = []
        for dat in self.data:
            dif = np.subtract(point, dat[0:2])
            dist = la.norm(dif)
            distances.append(dist)
        distances = sorted(distances)

        return distances[m]

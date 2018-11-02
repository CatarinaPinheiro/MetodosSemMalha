import numpy as np
from numpy import linalg as la


def gaussian_with_radius(dist, r, derivative=None):
    c = 100
    exp1 = np.exp(-((la.norm(dist) / c) ** 2))  # Exponencial 1
    exp2 = np.exp(-((r / c) ** 2))  # Exponencial 2
    # Função Peso normal Gaussiana com raio:
    if not derivative:
        if la.norm(dist) <= r:
            weight = (exp1 - exp2) / (1 - exp2)
            return weight
        else:
            return 0

    # Derivadas dessa função peso:
    else:
        c1 = 1 / (1 - exp2)

        if derivative['var'] == 'x':
            axis = dist[0]
        elif derivative['var'] == 'y':
            axis = dist[1]
        elif derivative['var'] == 'xy':
            axis = dist[0] * dist[1]
        else:
            axis = input('axis =')

        # Primeira derivada da função peso
        d1 = -2 * c1 * exp1 * axis / (c ** 2)

        d2 = -2 * c1 * (c ** 2 - 2 * axis ** 2) * exp1 / (c ** 4)  # Segunda derivada da função peso

        dxy = 4 * c1 * axis * exp1 / (c ** 4)

        if derivative['order'] == 1:
            return d1
        elif derivative['order'] == 2 and derivative['var'] != 'xy':
            return d2
        elif derivative['order'] == 2 and derivative['var'] == 'xy':
            return dxy

import coefficients as cf
import numpy as np


def interpol(data, base,point,contour_point, derivative=None):
    c = cf.coefficients(data, point, base,contour_point, derivative)
    u = np.array([z for x, y, z in data])
   # for i in range(len(interpol_data)):
       # print('Ponto',i, 'Coef',coefficients(data,point,base))
    return np.matmul(c, u)

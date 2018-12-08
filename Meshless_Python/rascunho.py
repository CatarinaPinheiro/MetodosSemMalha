import numpy as np
from test.helpers import weightmatrix_test as wm
from numpy import linalg as la
import basismatrix as bm
import basis as b
import coefficients as c
import minimumradius as mr

r_data = [[0, 1], [0.5, 1.5], [1, 2], [1.5, 2.5], [2, 3], [2.5, 3.5], [3, 4], [3.5, 4.5], [4, 5]]
contour_point = [[0, 0], [4, 4]]

radius_approximation = mr.get_radius(r_data, r_data[4], 3, contour_point)
print(radius_approximation)

base = b.pde_basis(1)[0]
P = bm.create_basis(base, r_data, r_data[4], radius_approximation)
Pt = np.transpose(P)
Peso = wm.W(r_data, r_data[4], radius_approximation)

print('P = ', P)
print('Pt = ', Pt)
print('W = ', Peso)

A = Pt @ Peso @ P

print('A = ', A)

B_matrix = Pt @ Peso

pt = bm.create_basis(base, [r_data[4]])
print('pt = ', pt)

phi = pt @ la.inv(A) @ B_matrix

print('phi = ', phi)

coef = c.coefficients(r_data, r_data[4], 1, contour_point)

print('coefficients = ', coef)

dptd_ = bm.create_basis(b.pde_basis(1)[1], [r_data[4]])
dWd_ = wm.W(r_data, r_data[4], radius_approximation, {
    'order': 1,
    'var': 'x'
})
dAd_ = Pt @ dWd_ @ P
dBd_ = Pt @ dWd_
invA = la.inv(A)

d1 = dptd_ @ invA @ B_matrix - pt @ invA @ dAd_ @ invA @ B_matrix + pt @ invA @ dBd_
print('derivative = ', d1)

pde_differential_x = {  # Information about pde equation
    'order': 1,  # order of pde
    'var': 'x',  # Variable
    'base1': b.pde_basis(1)[1],
    'base2': b.pde_basis(1)[3]
}

dcoef_1 = c.coefficients(r_data, r_data[4], 1, contour_point, pde_differential_x)

print('coefficients derivative = ', dcoef_1)

dptd_2 = bm.create_basis(b.pde_basis(1)[4], [r_data[4]])
d2Wd_ = wm.W(r_data, r_data[4], radius_approximation, {
    'order': 2,
    'var': 'y'
})
d2Bd_2 = Pt @ d2Wd_
d2Ad_2 = d2Bd_2 @ P

d2 = dptd_2 @ invA @ B_matrix + np.array([[
    2]]) @ pt @ invA @ dAd_ @ invA @ dAd_ @ invA @ B_matrix - pt @ invA @ d2Ad_2 @ invA @ B_matrix + pt @ invA @ d2Bd_2 - np.array(
    [[2]]) @ dptd_ @ invA @ dAd_ @ invA @ B_matrix + np.array([[2]]) @ dptd_ @ invA @ dBd_ - np.array(
    [[2]]) @ pt @ invA @ dAd_ @ invA @ dBd_

print('d2 =', d2)

pde_differential_xx = {  # Information about pde equation
    'order': 2,  # order of pde
    'var': 'x',  # Variable
    'base1': b.pde_basis(1)[1],
    'base2': b.pde_basis(1)[4]
}

coef_2 = c.coefficients(r_data, r_data[4], 1, contour_point, pde_differential_xx)
print('coef_2 = ', coef_2)


print(np.log(2))

print('phi = ', phi)
print('coefficients = ', coef)
import numpy as np
from src.helpers import weightmatrix as wm
import src.basis as b
import scipy as sp


class MovingLeastSquares2D:
    def __init__(self, data, basis_order, m,  derivative = None):
        self.basis_order = basis_order
        self.data = data
        self.point = np.zeros(np.shape(data[0]))
        self.derivative = derivative
        self.m = m


    @property
    def r_min(self):
        distances = [np.linalg.norm(np.subtract(d, self.point)) for d in self.data]
        return np.sort(distances)[self.m + 1]

    def AB(self, r):
        P = b.d2d.polynomial_basis.PolynomialBasis(self.basis_order).matrix(self.data, self.point, r)
        B = sp.transpose(P) @ wm.W(self.data, self.point, r, self.derivative)
        A = B @ P
        return A, B

    def compute_phi(self):
        pt = b.d2d.polynomial_basis.PolynomialBasis(self.basis_order).matrix([self.point])

        ri = self.r_min
        while np.linalg.det(
                np.array(self.AB(ri))) < 1e-9:
            ri *= 1.05
            print(ri)
        A, B = self.AB(ri)

        print(A)

        if self.derivative is None:
            invA = sp.sparse.linalg.splu(A)

            self.phi = pt @ invA @ B
            return self.phi
        else:
            # For dx or dy

            P = b.d2d.polynomial_basis.PolynomialBasis(self.basis_order).matrix(self.data, self.point, ri)
            Pt = np.transpose(P)
            dptd_ = b.d2d.polynomial_basis.PolynomialBasis(self.basis_order).matrix([self.point], derivative=self.derivative)
            dWd_ = wm.W(self.data, self.point, ri, self.derivative)
            dAd_ = Pt @ dWd_ @ P
            dBd_ = Pt @ dWd_
            invA = la.inv(A)  # A inverse
            (dist, r, derivative=None):
            # For dx² or dy²

            dptd_2 = bm.create_basis(derivative, [self.point])
            d2Wd_ = wm.W(self.data, self.point, ri, {
                'order': 2,
                'var': self.derivative['var']
            })
            d2Bd_2 = Pt @ d2Wd_
            d2Ad_2 = d2Bd_2 @ P

            # For dxy and dyx:
            dptd_x = bm.create_basis(b.pde_basis(basis_order)[1], [point])
            dptd_y = bm.create_basis(b.pde_basis(basis_order)[2], [point])

            dxyWd_ = wm.W(data, point, r, derivative)
            dxWd_ = wm.W(data, point, r, 'x')
            dyWd_ = wm.W(data, point, r, 'y')

            dxyBd = Pt @ dxyWd_
            dxBd = Pt @ dxWd_
            dyBd = Pt @ dyWd_
            dxyAd = dxyBd @ P
            dxAd = dxBd @ P
            dyAd = dyBd @ P

            # First derivative:

            d1 = dptd_ @ invA @ B - pt @ invA @ dAd_ @ invA @ B + pt @ invA @ dBd_

            # Second derivative:

            d2 = dptd_2 @ invA @ B + np.array([[
                2]]) @ pt @ invA @ dAd_ @ invA @ dAd_ @ invA @ B - pt @ invA @ d2Ad_2 @ invA @ B + pt @ invA @ d2Bd_2 - np.array(
                [[2]]) @ dptd_ @ invA @ dAd_ @ invA @ B + np.array([[2]]) @ dptd_ @ invA @ dBd_ - np.array(
                [[2]]) @ pt @ invA @ dAd_ @ invA @ dBd_

            # Derivative for xy:

            dxy = dptd_2 @ invA @ B - dptd_y @ invA @ dxAd @ invA @ B + dptd_y @ invA @ dxBd - dptd_x @ invA @ dyAd @ invA @ B + pt @ invA @ dxAd @ invA @ dyAd @ invA @ B - pt @ invA @ dxyAd @ invA @ B + pt @ invA @ dyAd @ invA @ dxAd @ invA @ B - pt @ invA @ dyAd @ invA @ dxBd + dptd_x @ invA @ dyBd - pt @ invA @ dxAd @ invA @ dyBd + pt @ invA @ dxyBd

            if derivative['order'] == 1:
                return d1
            elif derivative['order'] == 2 and derivative['var'] != 'xy':
                return d2
            elif derivative['order'] == 2 and derivative['var'] == 'xy':
                return dxy
        break


    @property
    def numeric_phi(self):
        dict = {
            'x': self.point[0],
            'y': self.point[1]
        }

        pt = np.array([sp.Matrix(self.basis).evalf(subs=dict)], dtype=np.float64)

        ri = self.r_min
        while np.linalg.det(np.array(self.AB(ri)[0].evalf(subs=dict), dtype=np.float64)) < 1e-6:
            ri *= 1.05
        A, B = self.AB(ri)

        return pt @ np.linalg.inv(np.array(A.evalf(subs=dict), dtype=np.float64)) @ np.array(B.evalf(subs=dict),
                                                                                             dtype=np.float64)

    def set_point(self, point):
        self.point = point

    def approximate(self, u):
        return self.numeric_phi @ u


import src.basis.d2d.polynomial_basis as pb
import numpy as np
from numpy import linalg as la
import minimumradius as mr
from test.helpers import weightmatrix_test as wm


def coefficients(data, point, basis_order, contour_point, derivative=None):
    basis = b.pde_basis(basis_order)[0]
    m = len(basis)
    r = mr.get_radius(data, point, m, contour_point)

    while True:
        P = pb.Polynomial(2)
            bm.create_basis(basis, data, point, r)  # Basis matrix
        Pt = np.transpose(P)
        weight_ = wm.W(data, point, r)
        pt = bm.create_basis(basis, [point])
        B = Pt @ weight_
        A = B @ P
        _, det = la.slogdet(A)
        determinante = det
        print(determinante)
        if det < np.log(1e-6) and m < len(data) - 1:
            r *= 1.05
            continue
        else:
            pass

        if not derivative:
            return pt @ la.inv(A) @ B

        else:

            # For dx or dy

            dptd_ = bm.create_basis(derivative, [point])
            dWd_ = wm.W(data, point, r, derivative)
            dAd_ = Pt @ dWd_ @ P
            dBd_ = Pt @ dWd_
            invA = la.inv(A)  # A inverse

            # For dx² or dy²

            dptd_2 = bm.create_basis(derivative, [point])
            d2Wd_ = wm.W(data, point, r, {
                'order': 2,
                'var': derivative['var']
            })
            d2Bd_2 = Pt @ d2Wd_
            d2Ad_2 = d2Bd_2 @ P

            # For dxy and dyx:
            dptd_x = bm.create_basis(b.pde_basis(basis_order)[1], [point])
            dptd_y = bm.create_basis(b.pde_basis(basis_order)[2], [point])

            dxyWd_ = wm.W(data, point, r, derivative)
            dxWd_ = wm.W(data, point, r, 'x')
            dyWd_ = wm.W(data, point, r, 'y')

            dxyBd = Pt @ dxyWd_
            dxBd = Pt @ dxWd_
            dyBd = Pt @ dyWd_
            dxyAd = dxyBd @ P
            dxAd = dxBd @ P
            dyAd = dyBd @ P

            # First derivative:

            d1 = dptd_ @ invA @ B - pt @ invA @ dAd_ @ invA @ B + pt @ invA @ dBd_

            # Second derivative:

            d2 = dptd_2 @ invA @ B + np.array([[
                2]]) @ pt @ invA @ dAd_ @ invA @ dAd_ @ invA @ B - pt @ invA @ d2Ad_2 @ invA @ B + pt @ invA @ d2Bd_2 - np.array(
                [[2]]) @ dptd_ @ invA @ dAd_ @ invA @ B + np.array([[2]]) @ dptd_ @ invA @ dBd_ - np.array(
                [[2]]) @ pt @ invA @ dAd_ @ invA @ dBd_

            # Derivative for xy:

            dxy = dptd_2 @ invA @ B - dptd_y @ invA @ dxAd @ invA @ B + dptd_y @ invA @ dxBd - dptd_x @ invA @ dyAd @ invA @ B + pt @ invA @ dxAd @ invA @ dyAd @ invA @ B - pt @ invA @ dxyAd @ invA @ B + pt @ invA @ dyAd @ invA @ dxAd @ invA @ B - pt @ invA @ dyAd @ invA @ dxBd + dptd_x @ invA @ dyBd - pt @ invA @ dxAd @ invA @ dyBd + pt @ invA @ dxyBd

            if derivative['order'] == 1:
                return d1
            elif derivative['order'] == 2 and derivative['var'] != 'xy':
                return d2
            elif derivative['order'] == 2 and derivative['var'] == 'xy':
                return dxy
        break


import numpy as np
from src.helpers import weightmatrix as wm
import src.basis.d2d.polynomial_basis as b
from scipy.sparse import linalg as sla


class MovingLeastSquares2D:
    def __init__(self, data, basis_order, m, derivative=None):
        self.basis_order = basis_order
        self.basis = b.PolynomialBasis(self.basis_order)
        self.data = data
        self.derivative = derivative
        self.m = m

    @property
    def r_min(self):
        distances = [np.linalg.norm(np.subtract(d, self.point)) for d in self.data]
        return np.sort(distances)[self.m + 1]

    def AB(self, r):
        P = self.basis.matrix(self.data, self.point, r)
        B = np.transpose(P) @ wm.W(self.data, self.point, r, self.derivative)
        A = B @ P
        return A, B

    @property
    def numeric_phi(self):
        pt = self.basis.matrix([self.point])

        ri = self.r_min

        while True:
            if np.linalg.det(np.array(self.AB(ri))) < 1e-9:
                ri *= 1.05
                continue

            A, B = self.AB(ri)

            print(A)

            if self.derivative is None:
                invA = sla.splu(A)
                return pt @ invA @ B

            else:
                # For dx or dy

                P = self.basis.matrix(self.data, self.point, ri)
                Pt = np.transpose(P)

                dptd_ = self.basis.matrix([self.point], derivative={'order': 1, 'var': self.derivative})
                dWd_ = wm.W(self.data, self.point, ri, {
                    'order': 1,
                    'var': self.derivative})
                dAd_ = Pt @ dWd_ @ P
                dBd_ = Pt @ dWd_
                invA = sla.splu(A)  # A inverse
                # For dx² or dy²

                dptd_2 = self.basis.matrix([self.point], derivative={'order': 2, 'var': self.derivative})
                d2Wd_ = wm.W(self.data, self.point, ri, {
                    'order': 2,
                    'var': self.derivative
                })
                d2Bd_2 = Pt @ d2Wd_
                d2Ad_2 = d2Bd_2 @ P

                # For dxy and dyx:
                dptd_x = self.basis.matrix([self.point], derivative={'order': 1, 'var': 'x'})
                dptd_y = self.basis.matrix([self.point], derivative={'order': 1, 'var': 'y'})

                dxyWd_ = wm.W(self.data, self.point, ri, {
                    'order': 2,
                    'var': self.derivative})
                dxWd_ = wm.W(self.data, self.point, ri, {
                    'order': 1,
                    'var': 'x'})
                dyWd_ = wm.W(self.data, self.point, ri, {
                    'order': 1,
                    'var': 'y'})

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
                    2]]) @ pt @ invA @ dAd_ @ invA @ dAd_ @ invA @ B - pt @ invA @ d2Ad_2 @ invA @ B +\
                    pt @ invA @ d2Bd_2 - np.array([[2]]) @ dptd_ @ invA @ dAd_ @ invA @ B +\
                    np.array([[2]]) @ dptd_ @ invA @ dBd_ - np.array([[2]]) @ pt @ invA @ dAd_ @ invA @ dBd_

                # Derivative for xy:
                dxy = dptd_2 @ invA @ B - dptd_y @ invA @ dxAd @ invA @ B + dptd_y @ invA @ dxBd - \
                    dptd_x @ invA @ dyAd @ invA @ B + pt @ invA @ dxAd @ invA @ dyAd @ invA @ B - \
                    pt @ invA @ dxyAd @ invA @ B + pt @ invA @ dyAd @ invA @ dxAd @ invA @ B - \
                    pt @ invA @ dyAd @ invA @ dxBd + dptd_x @ invA @ dyBd - pt @ invA @ dxAd @ invA @ dyBd +\
                    pt @ invA @ dxyBd

                if self.derivative == 'x' or self.derivative == 'y':
                    return d1
                elif self.derivative == 'xx' or self.derivative == 'yy':
                    return d2
                elif self.derivative == 'xy':
                    return dxy
                else:
                    return input('Phi Derivative: ')
            break

    def set_point(self, point):
        self.point = point

    def approximate(self, u):
        return self.numeric_phi @ u

    # @property
    # def numeric_phi(self):
    #     dict = {
    #         'x': self.point[0],
    #         'y': self.point[1]
    #     }
    #
    #     pt = np.array([Matrix(self.basis).evalf(subs=dict)], dtype=np.float64)
    #
    #     ri = self.r_min
    #     while np.linalg.det(np.array(self.AB(ri)[0].evalf(subs=dict), dtype=np.float64)) < 1e-6:
    #         ri *= 1.05
    #     A, B = self.AB(ri)
    #
    #     return pt @ np.linalg.inv(np.array(A.evalf(subs=dict), dtype=np.float64)) @ np.array(B.evalf(subs=dict),
    #                                                                                          dtype=np.float64)

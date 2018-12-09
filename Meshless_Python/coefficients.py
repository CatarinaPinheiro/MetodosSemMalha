import basis as b
import basismatrix as bm
import numpy as np
from numpy import linalg as la
import minimumradius as mr
import weightmatrix as wm


def coefficients(data, point, basis_order, derivative=None):
    basis = b.pde_basis(basis_order)[0]
    m = len(basis)
    r = mr.get_radius(data, point, m)

    while True:
        P = bm.create_basis(basis, data, point, r)  # Basis matrix
        Pt = np.transpose(P)
        weight_ = wm.W(data, point, r)
        pt = bm.create_basis(basis, [point])
        B = Pt @ weight_
        A = B @ P
        _, det = la.slogdet(A)
        determinante = det
        if det < np.log(1e-6) and m < len(data) - 1:
            r *= 1.05
            continue
        else:
            pass

        if not derivative:
            return pt @ la.inv(A) @ B

        else:

            # For dx or dy

            dptd_ = bm.create_basis(derivative['base1'], [point])
            dWd_ = wm.W(data, point, r, {
                'order': 1,
                'var': derivative['var']
            })
            dAd_ = Pt @ dWd_ @ P
            dBd_ = Pt @ dWd_
            invA = la.inv(A)  # A inverse

            # For dx² or dy²

            dptd_2 = bm.create_basis(derivative['base2'], [point])
            d2Wd_ = wm.W(data, point, r, {
                'order': 2,
                'var': derivative['var']
            })
            d2Bd_2 = Pt @ d2Wd_
            d2Ad_2 = d2Bd_2 @ P

            # For dxy and dyx:
            dptd_x = bm.create_basis(b.pde_basis(basis_order)[1], [point])
            dptd_y = bm.create_basis(b.pde_basis(basis_order)[2], [point])

            dxyWd_ = wm.W(data, point, r, {
                'order': 2,
                'var': derivative['var']
            })
            dxWd_ = wm.W(data, point, r, {
                'order': 1,
                'var': 'x'
            })
            dyWd_ = wm.W(data, point, r, {
                'order': 1,
                'var': 'y'
            })

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

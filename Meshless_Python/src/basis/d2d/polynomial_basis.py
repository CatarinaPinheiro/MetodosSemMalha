import numpy as np
from scipy import linalg as la


class PolynomialBasis:  # Defines basis configuration of pde.
    def __init__(self, basis_order):
        
        if basis_order == 1:  # Linear Basis
            self.basis = ["1", "x", "y"]
            self.basis_x = ["0", "1", "0"]  # First derivative in x
            self.basis_y = ["0", "0", "1"]
            self.basis_xx = ["0", "0", "0"]  # Second derivative in x
            self.basis_yy = ["0", "0", "0"]
            self.basis_xy = ["0", "0", "0"]
        elif basis_order == 2:  # Quadratic Basis
            self.basis = ["1", "x", "x**2", "y", "y**2", "x*y"]
            self.basis_x = ["0", "1", "2*x", "0", "0", "y"]
            self.basis_y = ["0", "0", "0", "1", "2*y", "x"]
            self.basis_xx = ["0", "0", "2", "0", "0", "0"]
            self.basis_yy = ["0", "0", "0", "0", "2", "0"]
            self.basis_xy = ["0", "0", "0", "0", "0", "1"]
        elif basis_order == 3:  # Cubic Basis
            self.basis = ["1", "x", "y", "x**2", "x*y", "y**2", "x**3", "(x**2)*y", "x*(y**2)", "y**3"]
            self.basis_x = ["0", "1", "0", "2*x", "y", "0", "3*(x**2)", "2*x*y", "(y**2)", "0"]
            self.basis_y = ["0", "0", "1", "0", "x", "2*y", "0", "(x**2)", "2*x*y", "3*(y**2)"]
            self.basis_xx = ["0", "0", "0", "2", "0", "0", "6*x", "2*y", "0", "0"]
            self.basis_yy = ["0", "0", "0", "0", "0", "2", "0", "0", "2*x", "6*y"]
            self.basis_xy = ["0", "0", "0", "0", "1", "0", "0", "2*x", "2*y", "0"]
    
        else:
            input('Pde order: ')

    def matrix(self,  data, point=None, r=None,  derivative=None):
        if derivative is None:
            basis = self.basis
        if derivative['order'] == 1:
            if derivative['var'] == 'x':
                basis = self.basis_x
            elif derivative == 'y':
                basis = self.basis_y
        elif derivative['order'] == 2:
            if derivative == 'xx':
                basis = self.basis_xx
            elif derivative == 'yy':
                basis = self.basis_yy
            elif derivative == 'xy':
                basis = self.basis_xy
        else:
            basis = input('Basis: ')

        P = []  # P = Basis function
        for dat in data:
            row = []
            for b in basis:
                [x, y] = dat[0:2]
                if point is not None and la.norm(np.subtract(point, [x, y])) > r:
                    row.append(0)
                else:
                    row.append(eval(b))
            P.append(row)
        return P

import gaussian_quadrature as gq
import numpy as np
integral = gq.polar_gauss_integral([0,0],1,lambda p: 1)
print("error = %f" % (np.pi-integral))
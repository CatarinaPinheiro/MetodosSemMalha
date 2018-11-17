import galerkin_potential as cp
import result_vector_potential as rvp
import nodalcoordinates as nc
from scipy import linalg as la

size_ix = 1  # Initial size in y domain
size_iy = 1  # Initial size in x domain
size_x = 6  # Size in x domain
size_y = 6  # Size in y domain

k = 4  # Domain division
kc = 4

# Domain:

domain = nc.cartesianproduct(nc.inside(nc.segments(size_ix, size_x, k)), nc.inside(nc.segments(size_iy, size_y, k)))

# Contour:

contour_x1 = [[size_ix, y] for y in nc.segments(size_iy, size_y, kc)]
contour_xf = [[size_x, y] for y in nc.segments(size_iy, size_y, kc)]
contour_y1 = [[x, size_iy] for x in nc.inside(nc.segments(size_ix, size_x, kc))]
contour_yf = [[x, size_y] for x in nc.inside(nc.segments(size_ix, size_x, kc))]

contour = contour_x1 + contour_xf + contour_y1 + contour_yf

pde_domain_count = len(domain)
pde_contour_count = len(contour)  # (kc-1)*4


# Domínio + Contorno

pde_data = domain + contour

print(len(pde_data))

phi = cp.lphi(pde_data, size_ix, size_x, size_iy, size_y, 1, contour)

b = rvp.b_matrix(pde_data, size_ix, size_x, size_iy, size_y)

answer = la.solve(phi, b)
print(answer)

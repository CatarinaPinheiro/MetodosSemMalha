import coefficients as c
import basis as b
import pdefunc_potential as pdefunc
import numpy as np
import time

pde_base_x = b.pde_basis(1)[1]
pde_base_y = b.pde_basis(1)[2]
pde_base_xx = b.pde_basis(1)[3]
pde_base_yy = b.pde_basis(1)[4]

pde_differential_x = {
    'order': 2,
    'var': 'x',
    'base1': pde_base_x,
    'base2': pde_base_xx
}

pde_differential_y = {
    'order': 2,
    'var': 'y',
    'base1': pde_base_y,
    'base2': pde_base_yy
}

pde_contour_conditions = {
    'top': {
        'kind': 'dirichlet',
        'value': 100
    },
    'bottom': {
        'kind': 'dirichlet',
        'value': 10
    },
    'left': {
        'kind': 'neumann',
        'order': 1,
        'var': 'x',
        'base1': pde_base_x,
        'base2': pde_base_xx,
        'value': 0
    },
    'right': {
        'kind': 'neumann',
        'order': 1,
        'var': 'x',
        'base1': pde_base_x,
        'base2': pde_base_xx,
        'value': 0
    }
}


def lphi(data, first_x, last_x, first_y, last_y, basis_order):
    result = []
    for point in data:
        clazz = pdefunc.pde_contour_class(point, first_x, last_x, first_y, last_y)
        start_time = time.time()
        if pdefunc.pde_contour_class(point, first_x, last_x, first_y, last_y) is None:
            cx = c.coefficients(data, point, basis_order, pde_differential_x)
            cy = c.coefficients(data, point, basis_order, pde_differential_y)
            result.append(np.add(cx[0], cy[0]))
        else:
            cond = pde_contour_conditions[clazz]
            if cond['kind'] == 'dirichlet':
                cd = c.coefficients(data, point, basis_order)
                result.append(cd[0])
            elif cond['kind'] == 'neumann':
                cn = c.coefficients(data, point, basis_order, cond)
                result.append(cn[0])
        print(time.time() - start_time)
    return result

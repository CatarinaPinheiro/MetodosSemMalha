import coefficients as cf
import basis as b
from src.models import pdefunction as pdefunc
from enum import Enum


class ContourConditionKind(Enum):
    NEUMANN = 1
    DIRICHLET = 2


pde_base_x = b.pde_basis(2)[1]
pde_base_y = b.pde_basis(2)[2]
pde_base_xx = b.pde_basis(2)[3]
pde_base_yy = b.pde_basis(2)[4]
pde_base_xy = b.pde_basis(2)[5]

pde_differential_y2 = {
    'order': 2,
    'var': 'y',
    'base1': pde_base_y,
    'base2': pde_base_yy
}

pde_differential_x2 = {
    'order': 2,
    'var': 'x',
    'base1': pde_base_x,
    'base2': pde_base_xx
}

pde_differential_xy = {
    'order': 2,
    'var': 'xy',
    'base1': pde_base_x,
    'base2': pde_base_xy
}

# Objects with contour conditions

pde_contour_conditions_x = {
    'left': {
        'kind': ContourConditionKind.DIRICHLET  # Technical features about top object
    },
    'right': {
        'kind': ContourConditionKind.NEUMANN,
        'order': 1,
        'var': 'x',
        'base1': pde_base_x,
        'base2': pde_base_xx
    },
    'top': {
        'kind': ContourConditionKind.DIRICHLET  # ContourConditionKind.NEUMANN,
        # 'order': 1,
        # 'var': 'x',
        # 'base1': pde_base_x,
        # 'base2': pde_base_xx
    },
    'bottom': {
        'kind': ContourConditionKind.DIRICHLET  # ContourConditionKind.NEUMANN,
        # 'order': 1,
        # 'var': 'x',
        # 'base1': pde_base_x,
        # 'base2': pde_base_xx
    }
}

pde_contour_conditions_y = {
    'left': {
        'kind': ContourConditionKind.DIRICHLET  # ContourConditionKind.NEUMANN,
        # 'order': 1,
        # 'var': 'y',
        # 'base1': pde_base_y,
        # 'base2': pde_base_yy
    },
    'right': {
        'kind': ContourConditionKind.NEUMANN,  # ContourConditionKind NEUMANN
        'order': 1,
        'var': 'y',
        'base1': pde_base_y,
        'base2': pde_base_yy
    },
    'top': {
        'kind': ContourConditionKind.DIRICHLET  # ContourConditionKind.NEUMANN,
        # 'order': 1,
        # 'var': 'y',
        # 'base1': pde_base_y,
        # 'base2': pde_base_yy
    },
    'bottom': {
        'kind': ContourConditionKind.DIRICHLET  # ContourConditionKind.NEUMANN,
        # 'order': 1,
        # 'var': 'y',
        # 'base1': pde_base_y,
        # 'base2': pde_base_yy
    }
}


def lphi(young_mod, poisson_ratio, data, first_x, last_x, first_y, last_y, basis_order, contour):
    E = young_mod  # Young's modulus (MPa)
    ny = poisson_ratio  # Poisson's ratio
    const_1 = 1 / (2 - ny)
    const = 1 / (1 - ny)
    G = E / (2 * (1 + ny))

    w = 1  # Deformation in z axis (Plane stress)

    result = []
    for point in data:
        clazz = pdefunc.pde_contour_class(point, first_x, last_x, first_y, last_y)
        if clazz is None:  # Domain

            if w == 0:

                cx_2 = cf.coefficients(data, point, basis_order, contour, pde_differential_x2)
                cy_2 = cf.coefficients(data, point, basis_order, contour, pde_differential_y2)
                cxy = cf.coefficients(data, point, basis_order, contour, pde_differential_xy)

                row_1 = []
                row_2 = []

                for j in range(len(data)):
                    row_1.append(G * (cx_2[0][j] + cy_2[0][j] + const_1 * cx_2[0][j]))
                    row_1.append(G * const_1 * cxy[0][j])
                    row_2.append(G * const_1 * cxy[0][j])
                    row_2.append(G * (cy_2[0][j] + cx_2[0][j] + const_1 * cy_2[0][j]))

                result.append(row_1)
                result.append(row_2)
            else:
                cx_2 = cf.coefficients(data, point, basis_order, contour, pde_differential_x2)
                cy_2 = cf.coefficients(data, point, basis_order, contour, pde_differential_y2)
                cxy = cf.coefficients(data, point, basis_order, contour, pde_differential_xy)

                row_1 = []
                row_2 = []

                for j in range(len(data)):
                    row_1.append(G * (cx_2[0][j] + cy_2[0][j] + const * cx_2[0][j]))
                    row_1.append(G * const * cxy[0][j])
                    row_2.append(G * const * cxy[0][j])
                    row_2.append(G * (cy_2[0][j] + cx_2[0][j] + const * cy_2[0][j]))

                result.append(row_1)
                result.append(row_2)

        else:  # Contour

            cond_x = pde_contour_conditions_x[clazz]
            cond_y = pde_contour_conditions_y[clazz]
            if cond_x['kind'] == ContourConditionKind.DIRICHLET:
                if cond_y['kind'] == ContourConditionKind.DIRICHLET:
                    cx = cf.coefficients(data, point, basis_order, contour)
                    cy = cf.coefficients(data, point, basis_order, contour)
                    row_1 = []
                    row_2 = []
                    for j in range(len(data)):
                        row_1.append(cx[0][j])
                        row_1.append(0)
                        row_2.append(0)
                        row_2.append(cy[0][j])

                    result.append(row_1)
                    result.append(row_2)
                elif cond_y['kind'] == ContourConditionKind.NEUMANN:
                    cx = cf.coefficients(data, point, basis_order, contour)
                    cy = cf.coefficients(data, point, basis_order, contour, cond_y)
                    row_1 = []
                    row_2 = []
                    for j in range(len(data)):
                        row_1.append(cx[0][j])
                        row_1.append(0)
                        row_2.append(0)
                        row_2.append(cy[0][j])

                    result.append(row_1)
                    result.append(row_2)
                elif point == [0, 0]:
                    cx = cf.coefficients(data, point, basis_order, contour)
                    cy = cf.coefficients(data, point, basis_order, contour)
                    row_1 = []
                    row_2 = []
                    for j in range(len(data)):
                        row_1.append(cx[0][j])
                        row_1.append(0)
                        row_2.append(0)
                        row_2.append(cy[0][j])

                    result.append(row_1)
                    result.append(row_2)

            elif cond_x['kind'] == ContourConditionKind.NEUMANN:
                if cond_y['kind'] == ContourConditionKind.NEUMANN:
                    cx = cf.coefficients(data, point, basis_order, contour, cond_x)
                    cy = cf.coefficients(data, point, basis_order, contour, cond_y)
                    row_1 = []
                    row_2 = []
                    for j in range(len(data)):
                        row_1.append(cx[0][j])
                        row_1.append(0)
                        row_2.append(0)
                        row_2.append(cy[0][j])

                    result.append(row_1)
                    result.append(row_2)
                elif cond_y['kind'] == ContourConditionKind.DIRICHLET:
                    cx = cf.coefficients(data, point, basis_order, contour, cond_x)
                    cy = cf.coefficients(data, point, basis_order, contour)
                    row_1 = []
                    row_2 = []
                    for j in range(len(data)):
                        row_1.append(cx[0][j])
                        row_1.append(0)
                        row_2.append(0)
                        row_2.append(cy[0][j])

                    result.append(row_1)
                    result.append(row_2)
            else:
                print("não deveria passar por aqui")

    return result

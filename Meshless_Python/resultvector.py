from src.models import pdefunction as pdefunc


def b_matrix(data, first_x, last_x, first_y, last_y):
    bb = []
    for p in data:
        clazz = pdefunc.pde_contour_class(p, first_x, last_x, first_y, last_y)
        if clazz is None:
            bb.append(pdefunc.pde_domain_func_plane_stress_x(p))
            bb.append(pdefunc.pde_domain_func_plane_stress_y(p))
        else:
            bb.append(pdefunc.pde_contour_func_plane_stress_x(p, clazz))
            bb.append(pdefunc.pde_contour_func_plane_stress_y(p, clazz))
    return bb

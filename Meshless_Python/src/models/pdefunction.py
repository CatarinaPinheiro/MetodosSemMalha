# Defines pde class:
def pde_contour_class(p, first_x, last_x, first_y, last_y):  # p = point
    if p[1] == first_y:
        return 'bottom'
    elif p[1] == last_y:
        return 'top'
    elif p[0] == first_x:
        return 'left'
    elif p[0] == last_x:
        return 'right'
    else:
        return None


# Returns a function in contour:

def pde_contour_func_plane_stress_x(p, clazz):
    E = 1
    P = 1

    if clazz == 'bottom':
        return p[0]
    elif clazz == 'top':
        return p[0]
    elif clazz == 'left':
        return 0
    elif clazz == 'right':
        return P/E


def pde_contour_func_plane_stress_y(p, clazz):
    ny = 0.25
    P = 1
    E = 1
    if clazz == 'bottom':
        return 5/4
    elif clazz == 'top':
        return -5/4
    elif clazz == 'left':
        if p[1] == 0:
            return 0
        else:
            return -p[1]/4
    elif clazz == 'right':
        return -ny*P/E


# Returns a function in domain:

def pde_domain_func_plane_stress_x(p):
    return 0


def pde_domain_func_plane_stress_y(p):
    return 0

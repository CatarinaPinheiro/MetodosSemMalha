from src.helpers import weight as w
from test.helpers import weightmatrix_test as wm
import numpy as np

# Error Test - Gaussian with radius function:

point_test = [[-100, -100], [0, 0], [40, 40], [80, 80]]
radius = 80


def test_gaussian_out():
    assert w.gaussian_with_radius(point_test[0], radius) == 0


def test_gaussian_maximum():
    assert w.gaussian_with_radius(point_test[1], radius) == 1


def test_gaussian_medium():
    assert (w.gaussian_with_radius(point_test[2], radius) - 0.420675747852) <= 0.00001


def test_gaussian_limit():
    assert w.gaussian_with_radius(point_test[3], radius) == 0


def test_matrix_weight():
    a = np.array(wm.W(point_test, point_test[1], 80))
    b = np.array([[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0.420675747852, 0], [0, 0, 0, 0]])
    assert np.allclose(a, b) is True


# Error test- First Derivative:


def test_gaussian_derivative_x_out():
    assert (w.gaussian_with_radius(point_test[0], radius, {'order': 1, 'var': 'x'}) - 0.0057259621009) <= 0.0000001


def test_gaussian_derivative_y_out():
    assert (w.gaussian_with_radius(point_test[0], radius, {'order': 1, 'var': 'y'}) - 0.0057259621009) <= 0.0000001


def test_gaussian_derivative_x_maximum():
    assert w.gaussian_with_radius(point_test[1], radius, {'order': 1, 'var': 'x'}) == 0


def test_gaussian_derivative_y_maximum():
    assert w.gaussian_with_radius(point_test[1], radius, {'order': 1, 'var': 'y'}) == 0


def test_gaussian_derivative_x_medium():
    assert (w.gaussian_with_radius(point_test[2], radius, {'order': 1, 'var': 'x'}) + 0.0122891880563) <= 0.0000001


def test_gaussian_derivative_y_medium():
    assert (w.gaussian_with_radius(point_test[2], radius, {'order': 1, 'var': 'y'}) + 0.0122891880563) <= 0.0000001


def test_gaussian_derivative_x_limit():
    assert (w.gaussian_with_radius(point_test[3], radius, {'order': 1, 'var': 'x'}) + 0.00941088536234) <= 0.0000001


def test_gaussian_derivative_y_limit():
    assert (w.gaussian_with_radius(point_test[3], radius, {'order': 1, 'var': 'y'}) + 0.00941088536234) <= 0.0000001


def test_matrix_derivative_x_weight():
    a = np.array(wm.W(point_test, point_test[1], 80, {'order': 1, 'var': 'x'}))
    b = np.array([[0.0057259621009, 0, 0, 0], [0, 0, 0, 0], [0, 0, -0.0122891880563, 0], [0, 0, 0, -0.00941088536234]])
    assert np.allclose(a, b) is True


# Second Derivative:


def test_gaussian_derivative_xx_out():
    assert (w.gaussian_with_radius(point_test[0], radius, {'order': 2, 'var': 'x'}) - 0.000057259621009) <= 0.0000001


def test_gaussian_derivative_yy_out():
    assert (w.gaussian_with_radius(point_test[0], radius, {'order': 2, 'var': 'y'}) - 0.000057259621009) <= 0.0000001


def test_gaussian_derivative_xx_maximum():
    assert (w.gaussian_with_radius(point_test[1], radius, {'order': 2, 'var': 'x'}) + 0.000423094551838) <= 0.0000001


def test_gaussian_derivative_yy_maximum():
    assert (w.gaussian_with_radius(point_test[1], radius, {'order': 2, 'var': 'y'}) + 0.000423094551838) <= 0.0000001


def test_gaussian_derivative_xx_medium():
    assert (w.gaussian_with_radius(point_test[2], radius, {'order': 2, 'var': 'x'}) + 0.000208916196958) <= 0.0000001


def test_gaussian_derivative_yy_medium():
    assert (w.gaussian_with_radius(point_test[2], radius, {'order': 2, 'var': 'y'}) + 0.000208916196958) <= 0.0000001


def test_gaussian_derivative_xx_limit():
    assert (w.gaussian_with_radius(point_test[3], radius, {'order': 2, 'var': 'x'}) - 0.0000329380987683) <= 0.000001


def test_gaussian_derivative_yy_limit():
    assert (w.gaussian_with_radius(point_test[3], radius, {'order': 2, 'var': 'y'}) - 0.0000329380987683) <= 0.0000001


def test_matrix_derivative_xx_weight():
    a = np.array(wm.W(point_test, point_test[1], 80, {'order': 2, 'var': 'x'}))
    b = np.array([[0.000057259621009, 0, 0, 0], [0, -0.000423094551838, 0, 0], [0, 0, -0.000208916196958, 0],
                  [0, 0, 0, 0.0000329380987683]])
    assert np.allclose(a, b) is True


# Derivative xy


def test_gaussian_derivative_xy_out():
    assert (w.gaussian_with_radius(point_test[0], radius, {'order': 2, 'var': 'xy'}) - 0.000114519242018) <= 0.0000001


def test_gaussian_derivative_xy_maximum():
    assert w.gaussian_with_radius(point_test[1], radius, {'order': 2, 'var': 'xy'}) == 0


def test_gaussian_derivative_xy_medium():
    assert (w.gaussian_with_radius(point_test[2], radius, {'order': 2, 'var': 'xy'}) - 0.0000983135044506) <= 0.0000001


def test_gaussian_derivative_xy_limit():
    assert (w.gaussian_with_radius(point_test[3], radius, {'order': 2, 'var': 'xy'}) - 0.000150574165798) <= 0.0000001


def test_matrix_derivative_xy_weight():
    a = np.array(wm.W(point_test, point_test[1], 80, {'order': 2, 'var': 'xy'}))
    b = np.array(
        [[0.000114519242018, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0.0000983135044506, 0], [0, 0, 0, 0.000150574165798]])
    assert np.allclose(a, b) is True

import unittest
import numpy as np
import test.test_functions as tf
import coefficients as c


class TestMovingLeastSquare2D(unittest.TestCase):
    def template(self, example):
        point = np.array([2.5, 2.5])

        data = example.data
        example.point = point
        approx = c.coefficients(data, point, 2) @ np.array(data[:, 2])

        self.assertEqual(np.round(approx, 3)[0], np.round(example.eval(point), 3))

    def test_polynomial_phi(self):
        self.template(tf.PolynomialExample(10))

    def test_linear_phi(self):
        self.template(tf.LinearExample(10))

    def test_exponential_phi(self):
        self.template(tf.ExponentialExample(10))

    def test_trigonometric_phi(self):
        self.template(tf.TrigonometricExample(10))


if __name__ == '__main__':
    unittest.main()

from src.helpers import weight as w
import matplotlib.pyplot as plt
import numpy as np

# Plot Gaussian with radius and derivatives
fig = plt.figure("peso -gaussiana com raio")
xs = np.arange(-100, 100, 10)
ys = np.arange(-100, 100, 10)
X, Y = np.meshgrid(xs, ys)

zs = np.array([w.gaussian_with_radius([x, y], 80) for x, y in zip(np.ravel(X), np.ravel(Y))])
Z = zs.reshape(X.shape)

mesh_plot = fig.add_subplot(111, projection='3d')
mesh_plot.plot_surface(X, Y, Z)

fig = plt.figure("peso -gaussiana com raio - Primeira derivada")
xs = np.arange(-100, 100, 10)
ys = np.arange(-100, 100, 10)
X, Y = np.meshgrid(xs, ys)

zs = np.array([w.gaussian_with_radius([x, y], 80, {'order': 1, 'var': 'x'}) for x, y in zip(np.ravel(X), np.ravel(Y))])
Z = zs.reshape(X.shape)

mesh_plot = fig.add_subplot(111, projection='3d')
mesh_plot.plot_surface(X, Y, Z)

fig = plt.figure("peso -gaussiana com raio - Segunda derivada")
xs = np.arange(-100, 100, 10)
ys = np.arange(-100, 100, 10)
X, Y = np.meshgrid(xs, ys)

zs = np.array([w.gaussian_with_radius([x, y], 80, {'order': 2, 'var': 'x'}) for x, y in zip(np.ravel(X), np.ravel(Y))])
Z = zs.reshape(X.shape)

mesh_plot = fig.add_subplot(111, projection='3d')
mesh_plot.plot_surface(X, Y, Z)
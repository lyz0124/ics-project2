import numpy as np
import sounddevice as sd
import matplotlib.pyplot as plt
from math import pi
import time

# Constants

# -- basic 
E = 21e10 # Pascal, young's modulus
# circular section case 
r = 0.5e-3 # m, radius of section
I = pi * r**4 / 4 # m^4, area moment of inertia
rou = 7874 # kg/m^3, density
A = pi * r**2 # m^2 (=1e6 mm^2), sectional area
L = 4e-2 # m (=1e-2 cm), length of beam

# -- damping

# -- initial conditions
F = 3 # Newton

# -- simulation
time = 10
clockmax = 1000
dt = time / clockmax
nx = 100
dx = L / nx
X = np.arange(0, L, dx)

# Setup (Initial conditions)
U = np.zeros(nx)
H = F/(6*E*I) * X**3
assert U.shape == X.shape
# -- plot initial condition
# plt.plot(X, H)
# plt.axis('equal')
# plt.show()

# Iteration
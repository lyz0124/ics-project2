import numpy as np
import sounddevice as sd
import matplotlib.pyplot as plt
from math import pi
import time

E_IRON = 21e10 # Pascal, young's modulus
# circular section case 
RADIUS = 0.5e-3 # m, radius of section
ROU_IRON = 7874 # kg/m^3, density
LENGTH = 4e-2 # m (=1e-2 cm), length of beam

class Beam:
    def __init__(self, E, I, rou, A, L):
        self.E = E # Pascal, young's modulus
        self.I = I # m^4, area moment of inertia
        self.rou = rou # kg/m^3, density
        self.A = A # m^2 (=1e6 mm^2), sectional area
        self.L = L # m (=1e-2 cm), length of beam

    def set_simulation(self, time, clockmax, nx):
        self.time = time
        self.dt = time / clockmax
        self.clockmax = clockmax
        self.nx = nx
        self.dx = self.L / nx
        self.X = np.arange(0, self.L, self.dx)

    def set_initial_state(self, F):
        self.H = F/(6*self.E*self.I) * self.X**3
        self.U = np.zeros(self.nx)

    @property
    def get_X(self):
        return self.X

    @property
    def get_H(self):
        return self.H

    @property
    def get_U(self):
        return self.U

    def vibrate():
        pass

class CircularBeam(Beam):
    def __init__(self, E, r, rou, L):
        self.E = E # Pascal, young's modulus
        self.I = pi * r**4 / 4 # m^4, area moment of inertia
        self.rou = rou # kg/m^3, density
        self.A = pi * r**2 # m^2 (=1e6 mm^2), sectional area
        self.L = L # m (=1e-2 cm), length of beam

if __name__ == "__main__":
    b = CircularBeam(E_IRON, RADIUS, ROU_IRON, LENGTH)
    b.set_simulation(2, 1000, 100)
    b.set_initial_state(3)
    plt.plot(b.get_X, b.get_H)
    plt.axis("equal")
    plt.show()
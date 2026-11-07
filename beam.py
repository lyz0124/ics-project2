import numpy as np
import sounddevice as sd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from math import pi
import time

E_IRON = 21e10 # Pascal, young's modulus
# circular section case 
RADIUS = 0.5e-3 # m, radius of section
ROU_IRON = 7874 # kg/m^3, density
LENGTH = 4e-2 # m (=1e-2 cm), length of beam

def my_derivative_4th(y, dx):
    """Fourth derivative for a cantilever beam.

    Clamped at x = 0 (y = 0, dy/dx = 0) and free at x = L
    (d2y/dx2 = 0, d3y/dx3 = 0).  Ghost nodes are used at both ends so the
    boundary points are updated by the same physics as the interior instead
    of being frozen.
    """
    d4y = np.zeros_like(y)
    d4y[2:-2] = (y[4:] - 4*y[3:-1] + 6*y[2:-2] - 4*y[1:-3] + y[:-4]) / (dx**4)

    # Free end: mirror the two ghost nodes outside the beam so that both the
    # bending moment and the shear force vanish at the last node.
    y_n = 2*y[-1] - y[-2]          # ghost node at x = L + dx
    y_n1 = 2*y_n - 2*y[-2] + y[-3] # ghost node at x = L + 2*dx
    d4y[-2] = (y_n - 4*y[-1] + 6*y[-2] - 4*y[-3] + y[-4]) / (dx**4)
    d4y[-1] = (y_n1 - 4*y_n + 6*y[-1] - 4*y[-2] + y[-3]) / (dx**4)

    # Clamped end: y[0] = 0 and the even ghost node y[-1] = y[1] enforces
    # dy/dx = 0.  Node 0 itself stays fixed, so d4y[0] keeps its zero value.
    d4y[1] = (y[3] - 4*y[2] + 7*y[1] - 4*y[0]) / (dx**4)
    return d4y

class Beam:
    def __init__(self, E, I, rou, A, L):
        self.E = E # Pascal, young's modulus
        self.I = I # m^4, area moment of inertia
        self.rou = rou # kg/m^3, density
        self.A = A # m^2 (=1e6 mm^2), sectional area
        self.L = L # m (=1e-2 cm), length of beam

    def set_simulation(self, dt, clockmax, dx, nskip):
        self.dt = dt
        self.clockmax = clockmax
        self.dx = dx
        self.nx = int(self.L / dx)
        self.X = np.arange(self.nx) * self.dx
        self.nskip = nskip

        self.Hsave = np.zeros((clockmax//nskip+1, self.nx))
        self.tsave = np.zeros(clockmax//nskip+1)

    def set_initial_state(self, F):
        # F in Newtons
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

    def vibrate(self):
        for clock in range(self.clockmax):
            t = clock * self.dt
            constant = (self.E*self.I) / (self.rou*self.A)
            self.U -= self.dt * constant * my_derivative_4th(self.H, self.dx)
            self.H += self.dt * self.U
            # self.H[-1] = 2 * self.H[-2] - self.H[-3]
            # self.H[-2] = 2 * self.H[-3] - self.H[-4]

            # assert self.U[1] == 0 and self.U[-2] == 0
            assert not np.isnan(self.get_U[2])

            if clock % self.nskip == 0:
                self.tsave[clock // self.nskip] = t
                self.Hsave[clock // self.nskip] = self.H

class CircularBeam(Beam):
    def __init__(self, E, r, rou, L):
        self.E = E # Pascal, young's modulus
        self.I = pi * r**4 / 4 # m^4, area moment of inertia
        self.rou = rou # kg/m^3, density
        self.A = pi * r**2 # m^2 (=1e6 mm^2), sectional area
        self.L = L # m (=1e-2 cm), length of beam

def gen_sim_parameters(beam: Beam, nx=50, time=0.1, safety_factor=0.4, nframes=300):
    """
    - nx: Number of spatial grid points.
    - safety_factor: CFL safety multiplier (< 0.5) to guarantee numerical stability.
    - nframes: Target number of frames to save for rendering/animation.
    """
    dx = beam.L / nx
    # CFL limit: dt <= 0.5 * dx^2 * sqrt((rou * A) / (E * I))
    cfl_dt = 0.5 * (dx**2) * np.sqrt((beam.rou * beam.A) / (beam.E * beam.I))
    dt = cfl_dt * safety_factor
    clockmax = max(1, int(np.ceil(time / dt)))
    nskip = max(1, clockmax // nframes)

    return dt, clockmax, dx, nskip

class BeamPlot:
    def __init__(self, beam: Beam):
        self.beam = beam
        self.fig, self.ax = plt.subplots()
        self.line = self.ax.plot(beam.X, beam.Hsave[0])[0]
        self.ax.axis("equal")

    def update(self, frame):
        self.line.set_ydata(self.beam.Hsave[frame])
        return (self.line,)

    def animate(self):
        self.ani = animation.FuncAnimation(fig=self.fig, func=self.update, 
                                      frames=self.beam.Hsave.shape[0], interval=30,
                                      blit=True)
        plt.show()

if __name__ == "__main__":
    b = CircularBeam(E_IRON, RADIUS, ROU_IRON, LENGTH)
    dt, clockmax, dx, nskip = gen_sim_parameters(b)
    b.set_simulation(dt, clockmax, dx, nskip)
    b.set_initial_state(3)

    b.vibrate()

    # print(b.get_H)
    # plt.plot(b.get_X, b.get_H)
    # plt.axis("equal")
    # plt.show()
    np.savetxt("vibbeam.csv", b.Hsave, delimiter=",")

    p = BeamPlot(b)
    p.animate()


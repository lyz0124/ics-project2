import numpy as np
from beam import *
import matplotlib.pyplot as plt

Hsave = np.load("output/b_Hsave.npy")
sound = np.reshape(Hsave[:, -2], -1)
soundsc(sound)
plt.plot(sound)
plt.show()


b = CircularBeam(E_IRON, RADIUS, ROU_IRON, LENGTH)
dt, clockmax, dx, nskip = gen_sim_parameters(b)
b.set_simulation(dt, clockmax, dx, nskip)
b.set_initial_state(3)
b.set_damping(8)
b.Hsave = Hsave

p = BeamPlot(b)
p.animate()
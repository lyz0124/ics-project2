import numpy as np
import matplotlib.pyplot as plt
from math import pi

from beam import *

def gen_circular_beam(freq) -> Beam:
    b = CircularBeam()
    dt, clockmax, dx, nskip = gen_sim_parameters(b)
    b.set_simulation(dt, clockmax, dx, nskip)
    b.set_initial_state(3)
    b.set_damping(8)
    return b

if __name__ == "__main__":
    pass
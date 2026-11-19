import numpy as np
import matplotlib.pyplot as plt
from math import pi, sqrt

from beam import *

BETA1L = 1.875104  # first root of cos(bL)cosh(bL) = -1 (cantilever)
FORCE = 3
GAMMA = 8
MIDC_FREQ = 256

def gen_beam_length(beam: Beam, freq):
    E, I, rou, A = beam.E, beam.I, beam.rou, beam.A
    L = BETA1L * sqrt(sqrt(E*I/(rou*A)) / (2*pi*freq))
    return L

def gen_circular_beam(freq, *, time=1.5, samplerate=8000, force=FORCE,
                      gamma=GAMMA, nx=50):
    beam = CircularBeam(E_IRON, RADIUS, ROU_IRON, LENGTH0)
    beam.L = gen_beam_length(beam, freq)
    dt, clockmax, dx, nskip = gen_sim_parameters(
        beam, nx=nx, time=time, samplerate=samplerate)
    beam.set_simulation(dt, clockmax, dx, nskip)
    beam.set_initial_state(force)
    beam.set_damping(gamma)
    return beam

if __name__ == "__main__":
    n = 12
    freq = []
    for i in range(n):
        freq.append(MIDC_FREQ*2**(i/12))
    beams: list[Beam] = []
    for i in range(n):
        beams.append(gen_circular_beam(freq[i]))

    for i in range(n):
        beams[i].vibrate()
        sound = beams[i].playsound()
        np.save(f"output/sound_{freq[i]:.2f}.npy", sound)
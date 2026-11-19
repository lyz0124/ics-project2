import numpy as np
import matplotlib.pyplot as plt
from math import pi, sqrt
import mido

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

def midi_to_numpy(midi_path: str):
    """Converts a single-channel MIDI file into NumPy arrays of notes, start times,

    durations, and amplitudes (velocities).

    Args:
        midi_path: Path to the .mid file.

    Returns:
        notes: 1D NumPy array of MIDI note numbers (0-127, e.g., 60 = Middle C).
        times: 1D NumPy array of note start times in seconds.
        durations: 1D NumPy array of note durations in seconds.
        amplitudes: 1D NumPy array of normalized note velocities (0.0 to 1.0).
    """
    mid = mido.MidiFile(midi_path)
    # Default MIDI tempo is 500,000 microseconds per beat (120 BPM)
    tempo = 500000
    ticks_per_beat = mid.ticks_per_beat
    # Track active notes to calculate durations: {note_number: (start_time_sec, velocity)}
    active_notes = {}
    completed_notes = []
    current_time_sec = 0.0

    for msg in mido.merge_tracks(mid.tracks):
        if msg.time > 0:
            current_time_sec += mido.tick2second(msg.time, ticks_per_beat, tempo)
        if msg.type == "set_tempo":
            tempo = msg.tempo
        elif msg.type == "note_on" and msg.velocity > 0:
            active_notes[msg.note] = (current_time_sec, msg.velocity)
        elif msg.type == "note_off" or (msg.type == "note_on" and msg.velocity == 0):
            if msg.note in active_notes:
                start_time, velocity = active_notes.pop(msg.note)
                duration = current_time_sec - start_time
                amplitude = velocity / 127.0
                completed_notes.append((msg.note, start_time, duration, amplitude))
    completed_notes.sort(key=lambda x: x[1])
    
    data = np.array(completed_notes)
    notes = data[:, 0].astype(np.int32)
    times = data[:, 1]
    durations = data[:, 2]
    amplitudes = data[:, 3]

    return notes, times, durations, amplitudes

if __name__ == "__main__":
    notes, times, durations, amplitudes = midi_to_numpy("ssg.mid")
    # print(notes, times, durations, amplitudes)
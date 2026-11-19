import numpy as np
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor, as_completed
from math import pi, sqrt
import time
import mido

from beam import *
from soundplayer import save_wav

BETA1L = 1.875104  # first root of cos(bL)cosh(bL) = -1 (cantilever)
FORCE = 3
GAMMA = 8
MIDC_FREQ = 256
SAMPLERATE = 8000
MAX_DURATION = 8.0  # a plucked beam is inaudible long before this
NX = 20  # spatial grid; <0.3% pitch error but ~7x faster than NX=50
MAX_SIM_TIME = 200.0  # only notes starting within this many seconds are
                      # simulated; 0 (or negative) keeps the whole score

def get_freq(note):
    return MIDC_FREQ * 2**((note-60)/12)

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


def pluck_waveform(note, *, duration=MAX_DURATION, samplerate=SAMPLERATE,
                   force=FORCE, gamma=GAMMA, nx=NX, position=None,
                   progress=True):
    """Simulate one pluck of the beam tuned to `note`, sampled near the free end.

    Every pluck uses the same force and starts from rest, so this waveform is
    identical each time the same beam is plucked.  It is therefore computed
    once per pitch and reused for all notes of that pitch.
    """
    beam = gen_circular_beam(get_freq(note), time=duration, samplerate=samplerate,
                             force=force, gamma=gamma, nx=nx)
    beam.vibrate(desc=f"beam {note:g}", position=position, progress=progress)
    return beam.sound


def simulate_beams(pitches, *, duration=MAX_DURATION, samplerate=SAMPLERATE,
                   force=FORCE, nx=NX, workers=None):
    """Simulate every beam in parallel, one worker process per pitch.

    The integrator in `Beam.vibrate` is a Python-level loop over small NumPy
    arrays, so it holds the GIL from start to finish: threads would serialise
    (and measure ~2x slower than a plain loop due to hand-off overhead).
    Processes do give real parallelism, at the cost of pickling each
    waveform back to the parent.  Results are returned in `pitches` order.
    """
    if len(pitches) == 0:
        return []

    workers = min(len(pitches), workers or len(pitches))
    results: list[np.ndarray] = [np.empty(0)] * len(pitches)
    start = time.perf_counter()

    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {
            pool.submit(pluck_waveform, p, duration=duration,
                        samplerate=samplerate, force=force, nx=nx,
                        progress=False): i
            for i, p in enumerate(pitches)
        }
        for done, future in enumerate(as_completed(futures), start=1):
            i = futures[future]
            results[i] = future.result()
            print(f"[{done}/{len(pitches)}] beam {pitches[i]:g} done "
                  f"after {time.perf_counter() - start:.1f} s", flush=True)

    return results


def main(midi_path="ssg.mid", out_path="output/music.wav", *,
         max_duration=MAX_DURATION, samplerate=SAMPLERATE, nx=NX,
         force=FORCE, max_sim_time=MAX_SIM_TIME, play=True):
    notes, times, _, _ = midi_to_numpy(midi_path)
    # MIDI note-off times and velocities are deliberately ignored: a beam
    # rings on its own and is always plucked with the same force.

    # 0. Cap the score: only plucks that happen inside the first `max_sim_time`
    #    seconds are simulated, everything later is dropped before any physics
    #    runs.  A non-positive cap means "no cap": every note is kept.
    if max_sim_time > 0:
        keep = times < max_sim_time
        cap_msg = f"{int(notes.size - int(keep.sum()))} notes after " \
                  f"{max_sim_time:g} s skipped"
        no_note_msg = (f"No notes start within the first {max_sim_time:g} s "
                       f"of {midi_path}")
    else:
        keep = np.ones(notes.size, dtype=bool)
        cap_msg = "no time cap"
        no_note_msg = f"No notes found in {midi_path}"

    notes, times = notes[keep], times[keep]
    if notes.size == 0:
        raise SystemExit(no_note_msg)

    # 1. An array of beams, one beam per distinct pitch in the score, each
    #    simulated in its own worker process.
    pitches = np.unique(notes)
    beams = simulate_beams(pitches, duration=max_duration,
                           samplerate=samplerate, force=force, nx=nx)
    beam_of = {int(p): i for i, p in enumerate(pitches)}

    # Note onsets for each beam, in chronological order (midi_to_numpy sorts
    # the notes by start time already).
    onsets = [[] for _ in pitches]
    for note, start in zip(notes, times):
        onsets[beam_of[int(note)]].append(float(start))

    # 2. Render every note's ring-out into a shared, polyphonic waveform.  A
    #    note is NOT cut off by its MIDI duration: it keeps ringing until the
    #    next note on the same beam re-plucks it, or until max_duration has
    #    elapsed (when the damping has made it effectively silent).
    first = float(times.min())
    n_samples = int(np.ceil((float(times.max()) - first + max_duration) * samplerate))
    mix = np.zeros(n_samples)

    for wave, starts in zip(beams, onsets):
        for k, start in enumerate(starts):
            next_note = starts[k + 1] if k + 1 < len(starts) else start + max_duration
            stop = min(next_note, start + max_duration)
            i0 = int(round((start - first) * samplerate))
            count = int(round((stop - start) * samplerate))
            mix[i0:i0 + count] += wave[:count]

    peak = np.max(np.abs(mix))
    if peak > 0:
        mix /= peak

    if play:
        soundsc(mix, samplerate)
    path = save_wav(mix, out_path, samplerate)
    print(f"Rendered {len(notes)} notes on {len(pitches)} beams in parallel "
          f"({mix.size / samplerate:.1f} s, {cap_msg}) -> {path}")
    return mix


if __name__ == "__main__":
    main()
import wave

import numpy as np
from beam import *
import matplotlib.pyplot as plt
from pathlib import Path

SAMPLERATE = 8000  # must match the samplerate used to generate the .npy files

# Hsave = np.load("output/b_Hsave.npy")
# sound = np.reshape(Hsave[:, -2], -1)
# soundsc(sound)
# plt.plot(sound)
# plt.show()


# b = CircularBeam(E_IRON, RADIUS, ROU_IRON, LENGTH0)
# dt, clockmax, dx, nskip = gen_sim_parameters(b)
# b.set_simulation(dt, clockmax, dx, nskip)
# b.set_initial_state(3)
# b.set_damping(8)
# b.Hsave = Hsave

# p = BeamPlot(b)
# p.animate()

def load_sounds() -> list[np.ndarray]:
    """Load all output/sound_*.npy files, ordered by their frequency."""
    soundfiles = Path("output").glob("sound_*.npy")
    # sort by the frequency encoded in the filename (sound_256.00.npy -> 256.0)
    soundfiles = sorted(soundfiles, key=lambda p: float(p.stem.removeprefix("sound_")))
    return [np.load(soundfile).ravel() for soundfile in soundfiles]


def normalize(sound: np.ndarray) -> np.ndarray:
    """Scale a sound so its peak absolute amplitude is 1.0 (like soundsc)."""
    peak = np.max(np.abs(sound))
    return sound / peak if peak > 0 else sound


def save_wav(sound: np.ndarray, path: str | Path = "output/all_sounds.wav",
             samplerate: int = SAMPLERATE) -> Path:
    """Write a float array in [-1, 1] to a 16-bit mono PCM .wav file."""
    path = Path(path)
    pcm = (np.clip(sound, -1.0, 1.0) * np.iinfo(np.int16).max).astype(np.int16)
    with wave.open(str(path), "wb") as w:
        w.setnchannels(1)
        w.setsampwidth(2)  # 16-bit
        w.setframerate(samplerate)
        w.writeframes(pcm.tobytes())
    return path


if __name__ == "__main__":
    sounds = load_sounds()
    if not sounds:
        raise SystemExit("No sound files found in output/")

    # concatenate every note into a single continuous signal
    waveform = np.concatenate([normalize(sound) for sound in sounds])

    outfile = save_wav(waveform)
    duration = waveform.size / SAMPLERATE
    print(f"Saved {len(sounds)} sounds ({duration:.2f} s) to {outfile}")

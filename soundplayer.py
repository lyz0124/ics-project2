import numpy as np
from beam import soundsc
import matplotlib.pyplot as plt

Hsave = np.load("output/b_Hsave.npy")
sound = np.reshape(Hsave[:, -2], -1)
soundsc(sound)
plt.plot(sound)
plt.show()
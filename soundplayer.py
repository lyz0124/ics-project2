import numpy as np
from beam import soundsc

sound = np.load("b_Xsave.npy")
soundsc(np.reshape(sound[:, -2], -1))
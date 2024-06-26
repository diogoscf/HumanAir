import numpy as np
import matplotlib.pyplot as plt

alpha = np.arange(-10, 21, 1)
Cm = [
    0.5423781,
    0.4849711,
    0.426917,
    0.3682909,
    0.3091593,
    0.2495903,
    0.1896587,
    0.1294319,
    0.06898292,
    0.008385102,
    -0.05228833,
    -0.1129614,
    -0.1735526,
    -0.2339882,
    -0.2942099,
    -0.3541424,
    -0.4137129,
    -0.4728568,
    -0.5314991,
    -0.5895758,
    -0.6470197,
    -0.7037585,
    -0.7597447,
    -0.8149012,
    -0.8691561,
    -0.9224438,
    -0.9746959,
    -1.025863,
    -1.075897,
    -1.124741,
    -1.172447,
]

plt.figure(figsize=(10, 7))
plt.plot(alpha, Cm)
plt.xlabel(r"$\alpha$, [deg]", fontsize=12, loc="right")
plt.ylabel(r"$C_m$, [-]", fontsize=12, loc="top")
plt.axhline(color="black", lw=0.5)
plt.axvline(color="black", lw=0.5)
plt.savefig("Stability_FlyingWing.svg")
plt.savefig("Stability_FlyingWing.png")
plt.show()

print((Cm[1] - Cm[0]) / (alpha[1] - alpha[0]) * 180 / np.pi)

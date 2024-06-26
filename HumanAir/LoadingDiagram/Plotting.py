import matplotlib.pyplot as plt
import numpy as np


def Plotx(xvalue, ylst, name: str, colour=None, linestyle=None):  # pragma: no cover
    xlst = np.ones(len(ylst)) * xvalue
    if colour and linestyle:
        plt.plot(xlst, ylst, label=name, color=colour, linestyle=linestyle)
    else:
        plt.plot(xlst, ylst, label=name)


def Ploty(xlst, ylst, name: str, colour=None, linestyle=None):  # pragma: no cover
    if colour and linestyle:
        plt.plot(xlst, ylst, label=name, color=colour, linestyle=linestyle)
    else:
        plt.plot(xlst, ylst, label=name)

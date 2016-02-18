import numpy as np
import pylab as pl
from math import log

class Radar(object):

    def __init__(self, fig, titles, labels, rect=None):
        if rect is None:
            rect = [0.13, 0.15, 0.72, 0.72]

        self.n = len(titles)
        self.angles = np.arange(90, 90+360, 360.0/self.n)
        self.axes = [fig.add_axes(rect, projection="polar", label="axes%d" % i) 
                         for i in range(self.n)]

        self.ax = self.axes[0]
        self.ax.set_thetagrids(self.angles, labels=titles, fontsize=14, frac=1.2)

        for ax in self.axes[1:]:
            ax.patch.set_visible(False)
            ax.grid("off")
            ax.xaxis.set_visible(False)

        for ax, angle, label in zip(self.axes, self.angles, labels):
            ax.set_rgrids(range(1, 4), angle=angle, labels=label)
            #ax.set_rgrids(range(1, 4), angle=angle)
            ax.spines["polar"].set_visible(False)
            #ax.set_ylim(0, 3)
            #ax.set_yscale('log')

    def plot(self, values, *args, **kw):
        angle = np.deg2rad(np.r_[self.angles, self.angles[0]])
        values = np.r_[values, values[0]]
        self.ax.plot(angle, values, *args, **kw)

def lerp(A, B, alpha):
	return A * (1.0 - alpha) + B * alpha

def normalizeValue(data, max, min):
	result = data;

	for i in range(len(data)):
		result[i] = ((data[i] - min[i]) * 3) / (max[i] - min[i])

	result[3] = 3.0 - result[3]

	return result


fig = pl.figure(figsize=(6, 6))

titles = ['ERROR', 'EUC', 'MSE', 'PSNR']

axMax = [4.000000E-009, 10.000000e-009, 2.000000e-021, 380.0]
axMin = [1.000000E-030, 1.000000e-015, 1.000000e-030, 300.0]

dataOri = [3.405148E-009, 9.21961e-009, 1.29702e-021, 305.2]
dataDev = [2.110028E-010, 5.93178e-011, 5.36896e-026, 349.031]

newDataOri = normalizeValue(dataOri, axMax, axMin)
newDataDev = normalizeValue(dataDev, axMax, axMin)

#for i in range(len(newDataOri) - 1):
#	newDataOri[i] = log(newDataOri[i])

#for i in range(len(newDataOri) - 1):
#	newDataOri[i] = log(newDataOri[i])

labels = [
    [lerp(axMin[0], axMax[0], 0.333), lerp(axMin[0], axMax[0], 0.666), lerp(axMin[0], axMax[0], 1.0)],
    [lerp(axMin[1], axMax[1], 0.333), lerp(axMin[1], axMax[1], 0.666), lerp(axMin[1], axMax[1], 1.0)],
    [lerp(axMin[2], axMax[2], 0.333), lerp(axMin[2], axMax[2], 0.666), lerp(axMin[2], axMax[2], 1.0)],
    [lerp(axMax[3], axMin[3], 0.333), lerp(axMax[3], axMin[3], 0.666), lerp(axMax[3], axMin[3], 1.0)]
]

radar = Radar(fig, titles, labels)
radar.plot(newDataOri,  "-", lw=2, color="b", alpha=0.4, label="first")
radar.plot(newDataDev,"-", lw=2, color="r", alpha=0.4, label="second")
#radar.ax.legend()
pl.show()
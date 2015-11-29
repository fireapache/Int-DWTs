import numpy as np
import pylab as pl

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
            ax.spines["polar"].set_visible(False)
            ax.set_ylim(0, 3)

    def plot(self, values, *args, **kw):
        angle = np.deg2rad(np.r_[self.angles, self.angles[0]])
        values = np.r_[values, values[0]]
        self.ax.plot(angle, values, *args, **kw)

fig = pl.figure(figsize=(5, 5))

titles = ['ERROR', 'EUC', 'MSE', 'PSNR']

axMax = [20, 15, 10, 5]
axMin = [5, 4, 3, 2]

dataOri = [10, 9, 8, 7]
dataDev = [9, 8, 7, 6]

newDataOri = dataOri
newDataDev = dataDev



labels = [
    [axMax[0] * 0.5, axMax[0] * 0.75, axMax[0]],
    [axMax[1] * 0.5, axMax[1] * 0.75, axMax[1]],
    [axMax[2] * 0.5, axMax[2] * 0.75, axMax[2]],
    [axMax[3] * 0.5, axMax[3] * 0.75, axMax[3]]
]

radar = Radar(fig, titles, labels)
radar.plot([1, 3, 2, 5],  "-", lw=2, color="b", alpha=0.4, label="first")
radar.plot([2.3, 2, 3, 3],"-", lw=2, color="r", alpha=0.4, label="second")
#radar.ax.legend()
pl.show()
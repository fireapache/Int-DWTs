"""
Bar chart demo with pairs of bars grouped for easy comparison.
"""
import numpy as np
import matplotlib.pyplot as plt


n_groups = 3

# fig, ax = plt.subplots(nrows=3, sharex=True, figsize=(5, 6))
fig = plt.figure(figsize=(5, 6))
fig.subplots_adjust(hspace=0.01)

ax = []

ax.append(fig.add_subplot(311))
ax.append(fig.add_subplot(312, sharex=ax[0]))
ax.append(fig.add_subplot(313))

# ax.append(fig.add_subplot(337))
# ax.append(fig.add_subplot(338))
# ax.append(fig.add_subplot(339))

index = np.arange(n_groups)
bar_width = 0.6

opacity = 1.0

# colors1 = (1.0, 0.6, 0.6)
# colors2 = (1.0, 0.6, 0.6)
colors1 = [(1.0, 0.6, 0.6), (1.0, 1.0, 0.0), (0.0, 1.0, 1.0)]
colors2 = [(0.6, 0.6, 1.0), (1.0, 0.6, 0.8), (0.6, 1.0, 0.6)]

error_config = {'ecolor': '0.3'}

plt.sca(ax[0])

developed = (27.8399, 58.12, 45.0949)
std_dev = (0.0, 0.0, 0.0)

plt.ylabel('Performance (%)')
plt.title('Decimated 2D Standard HWT')
plt.xticks(index, ('Decomposition', 'Composition', 'Decomp & Comp'))
plt.grid(True)
plt.ylim([0.0,70.0])
# ax[0].xaxis.set_visible(False)

Performance = ax[0].bar(index, developed, bar_width,
                 alpha=opacity,
                 color=colors1,
                 yerr=std_dev,
                 error_kw=error_config,
                 label='Developed',
                 align='center')

plt.xlim([min(index) - 0.5, max(index) + 0.5])

# plt.legend(loc='upper left', prop={'size':11})

plt.sca(ax[1])

developedError = (99.824766, 93.548387, 98.285199)

plt.ylabel('Accuracy (%)')
plt.xticks(index, ('Decomposition', 'Composition', 'Decomp & Comp'))
plt.grid(True)
plt.ylim([0.0,110.0])

Errors = ax[1].bar(index, developedError, bar_width,
                 alpha=opacity,
                 color=colors1,
                 #yerr=std_ori,
                 error_kw=error_config,
                 label='Developed',
                 align='center')

plt.sca(ax[2])

MetricResults = (99.8564, 99.9998, 19.1594)

plt.ylabel('Quality (%)')
plt.xticks(index, ('EUC', 'MSE', 'PSNR'))
plt.grid(True)
plt.ylim([0.0,110.0])
plt.xlim([min(index) - 0.5, max(index) + 0.5])

Metrics = ax[2].bar(index, MetricResults, bar_width,
                 alpha=opacity,
                 color=colors2,
                 #yerr=std_ori,
                 error_kw=error_config,
                 label='Developed',
                 align='center')

plt.tight_layout()
#plt.show()
fig.savefig("test2S.pdf", format='pdf')
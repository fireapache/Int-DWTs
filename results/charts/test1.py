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
rcolor = (1.0, 0.6, 0.6)
bcolor = (0.6, 0.6, 1.0)
error_config = {'ecolor': '0.3'}

plt.sca(ax[0])

developed = (18.9141, 37.2692, 29.08)
std_dev = (0.0, 0.0, 0.0)

plt.ylabel('Speedup (%)')
plt.title('Decimated 1D HWT')
plt.xticks(index, ('Decomposition', 'Composition', 'Decomp & Comp'))
plt.grid(True)
# plt.ylim([0.0,110.0])
# ax[0].xaxis.set_visible(False)

Speedup = ax[0].bar(index, developed, bar_width,
                 alpha=opacity,
                 color=rcolor,
                 yerr=std_dev,
                 error_kw=error_config,
                 label='Developed',
                 align='center')

plt.xlim([min(index) - 0.5, max(index) + 0.5])

# plt.legend(loc='upper left', prop={'size':11})

plt.sca(ax[1])

developedError = (99.764151, 95.833333, 95.107914)

plt.ylabel('Error Gain (%)')
plt.xticks(index, ('Decomposition', 'Composition', 'Decomp & Comp'))
plt.grid(True)
plt.ylim([0.0,110.0])

Errors = ax[1].bar(index, developedError, bar_width,
                 alpha=opacity,
                 color=rcolor,
                 #yerr=std_ori,
                 error_kw=error_config,
                 label='Developed',
                 align='center')

plt.sca(ax[2])

developedError = (99.4716, 99.9972, 13.8997)

plt.ylabel('Metrics Gain (%)')
plt.xticks(index, ('EUC', 'MSE', 'PSNR'))
plt.grid(True)
plt.ylim([0.0,110.0])
plt.xlim([min(index) - 0.5, max(index) + 0.5])

Metrics = ax[2].bar(index, developedError, bar_width,
                 alpha=opacity,
                 color=rcolor,
                 #yerr=std_ori,
                 error_kw=error_config,
                 label='Developed',
                 align='center')

plt.tight_layout()
plt.show()

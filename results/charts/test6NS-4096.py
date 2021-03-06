"""
Bar chart demo with pairs of bars grouped for easy comparison.
"""
import numpy as np
import matplotlib.pyplot as plt


n_groups = 3

# fig, ax = plt.subplots(nrows=3, sharex=True, figsize=(5, 6))
fig = plt.figure(figsize=(3, 4))
fig.subplots_adjust(hspace=0.01)

ax = []

ax.append(fig.add_subplot(211))
ax.append(fig.add_subplot(212))

#ax.append(fig.add_subplot(311))
#ax.append(fig.add_subplot(312, sharex=ax[0]))
#ax.append(fig.add_subplot(313))

patterns1 = ('/', '\\', '-', '|')
patterns2 = ('o', 'O', '.', '*')

index = np.arange(n_groups)
bar_width = 0.6

opacity = 1.0
# colors1 = (1.0, 0.6, 0.6)
# colors2 = (1.0, 0.6, 0.6)
colors1 = [(1.0, 0.6, 0.6), (1.0, 1.0, 0.0), (0.0, 1.0, 1.0)]
colors2 = [(0.6, 0.6, 1.0), (1.0, 0.6, 0.8), (0.6, 1.0, 0.6)]
error_config = {'ecolor': '0.3'}

plt.sca(ax[0])

developed = (-21.9542, -15.203, -17.3513)
std_dev = (0.0, 0.0, 0.0)

plt.ylabel('')
plt.title('TWD 2D Nao Padrao #3')
plt.xticks(index, ('D', 'C', 'D & C'))
plt.grid(True)
plt.ylim([-25.0,5.0])
# ax[0].xaxis.set_visible(False)

Performance = ax[0].bar(index, developed, bar_width,
                 alpha=opacity,
                 color=colors1,
                 yerr=std_dev,
                 error_kw=error_config,
                 label='Developed',
                 align='center')

# for bar, pattern in zip(Speedup, patterns1):
#     bar.set_hatch(pattern)

plt.xlim([min(index) - 0.5, max(index) + 0.5])
plt.axhline(y=0, xmin=min(index) - 0.5, xmax=max(index) + 0.5, color='black')

# plt.legend(loc='upper left', prop={'size':11})

# plt.sca(ax[1])

# developedError = (-3.425046E+004, -2.140426E+003, -2.036283E+007)

# plt.ylabel('Accuracy (%)')
# plt.xticks(index, ('D', 'Composition', 'Decomp & Comp'))
# plt.grid(True)
# plt.ylim([-0.3E+008,0.0])

# Errors = ax[1].bar(index, developedError, bar_width,
#                  alpha=opacity,
#                  color=colors1,
#                  #yerr=std_ori,
#                  error_kw=error_config,
#                  label='Developed',
#                  align='center')

plt.sca(ax[1])

MetricResults = (43.8987, 68.5265, 1.73043)


plt.ylabel('')
plt.xticks(index, ('EUC', 'MSE', 'PSNR'))
plt.grid(True)
plt.ylim([0.0,80.0])
plt.xlim([min(index) - 0.5, max(index) + 0.5])

Metrics = ax[1].bar(index, MetricResults, bar_width,
                 alpha=opacity,
                 color=colors2,
                 #yerr=std_ori,
                 error_kw=error_config,
                 label='Developed',
                 align='center')

plt.tight_layout()
#plt.show()
fig.savefig("test6NS-4096.pdf", format='pdf')
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
ax.append(fig.add_subplot(337))
ax.append(fig.add_subplot(338))
ax.append(fig.add_subplot(339))

index = np.arange(n_groups)
bar_width = 0.35

opacity = 1.0
rcolor = (1.0, 0.6, 0.6)
bcolor = (0.6, 0.6, 1.0)
error_config = {'ecolor': '0.3'}

plt.sca(ax[0])

original = (0.00109071, 0.00125803, 0.00236344)
std_ori = (3.71057e-005, 4.4377e-005, 0.000101865)

developed = (0.000805147, 0.00079437, 0.00157319)
std_dev = (4.34346e-005, 2.91344e-005, 5.13007e-005)

plt.ylabel('Time')
plt.title('Decimated 1D HWT')
plt.xticks(index + bar_width, ('Decomposition', 'Composition', 'Decomp & Comp'))
ax[0].xaxis.set_visible(False)

rects1 = ax[0].bar(index, original, bar_width,
                 alpha=opacity,
                 color=bcolor,
                 yerr=std_ori,
                 error_kw=error_config,
                 capsize=6.0,
                 label='Original')

rects2 = ax[0].bar(index + bar_width, developed, bar_width,
                 alpha=opacity,
                 color=rcolor,
                 yerr=std_dev,
                 error_kw=error_config,
                 label='Developed')

plt.legend(loc='upper left', prop={'size':11})

plt.sca(ax[1])

originalError = (3.128662E-010, 6.332994E-008, 3.405148E-009)
developedError = (1.136868E-012, 5.587935E-009, 2.110028E-010)

plt.ylabel('Error')
plt.xticks(index + bar_width, ('Decomposition', 'Composition', 'Decomp & Comp'))
plt.ylim([1.000000E-013,7.000000E-007])
plt.yscale('log')

rectError1 = ax[1].bar(index, originalError, bar_width,
                 alpha=opacity,
                 color=bcolor,
                 #yerr=std_ori,
                 error_kw=error_config,
                 label='Original')

rectError2 = ax[1].bar(index + bar_width, developedError, bar_width,
                 alpha=opacity,
                 color=rcolor,
                 #yerr=std_ori,
                 error_kw=error_config,
                 label='Developed')

plt.sca(ax[2])

oriEUC = (9.21961e-009)
devEUC = (5.93178e-011)

plt.ylabel('EUC')
plt.xticks(np.arange(0))
plt.ylim([1.000000E-012,7.000000E-008])
plt.yscale('log')

rectEUC1 = ax[2].bar(1, oriEUC, bar_width,
                 alpha=opacity,
                 color=bcolor,
                 #yerr=std_ori,
                 #error_kw=error_config,
                 label='Original')

rectEUC2 = ax[2].bar(1, devEUC, bar_width,
                 alpha=opacity,
                 color=rcolor,
                 #yerr=std_ori,
                 #error_kw=error_config,
                 label='Original')

plt.sca(ax[3])

oriMSE = (1.29702e-021)
devMSE = (5.36896e-026)

ylim = [1.000000E-027,7.000000E-021]

plt.ylabel('MSE')
plt.xticks(np.arange(0))
plt.ylim(ylim)
plt.yscale('log')

rectMSE1 = ax[3].bar(1, oriMSE, bar_width,
                 alpha=opacity,
                 color=bcolor,
                 #yerr=std_ori,
                 #error_kw=error_config,
                 label='Original')

rectMSE2 = ax[3].bar(1, devMSE, bar_width,
                 alpha=opacity,
                 color=rcolor,
                 #yerr=std_ori,
                 #error_kw=error_config,
                 label='Original')

plt.sca(ax[4])

oriPSNR = (305.2)
devPSNR = (349.031)

ylim = [280,360]

plt.ylabel('PSNR')
plt.xticks(np.arange(0))
plt.ylim(ylim)
plt.yticks(np.arange(min(ylim), max(ylim), 15.0))
#plt.yscale('log')

rectPSNR2 = ax[4].bar(1, devPSNR, bar_width,
                 alpha=opacity,
                 color=rcolor,
                 #yerr=std_ori,
                 #error_kw=error_config,
                 label='Original')

rectPSNR1 = ax[4].bar(1, oriPSNR, bar_width,
                 alpha=opacity,
                 color=bcolor,
                 #yerr=std_ori,
                 #error_kw=error_config,
                 label='Original')

plt.tight_layout()
plt.show()

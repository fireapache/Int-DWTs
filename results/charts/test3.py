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

original = (0.0273379, 0.0115317, 0.0383284)
std_ori = (0.000339537, 7.12352e-05, 0.000255277)

developed = (0.0415496, 0.0115378, 0.0523903)
std_dev = (0.000298855, 9.94945e-05, 0.000295155)

plt.ylabel('Time')
plt.title('Undecimated 1D HWT')
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

plt.legend(loc='upper center', prop={'size':11})

plt.sca(ax[1])

originalError = (1.891749E-010, 0.000000, 2.124580E-009)
developedError = (5.820766E-011, 0.000000, 4.074536E-010)

plt.ylabel('Error')
plt.xticks(index + bar_width, ('Decomposition', 'Composition', 'Decomp & Comp'))
plt.ylim([1.000000E-012,7.000000E-009])
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

oriEUC = (4.7575e-09)
devEUC = (2.16581e-09)

plt.ylabel('EUC')
plt.xticks(np.arange(0))
plt.ylim([8.000000E-010,1.000000E-008])
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

oriMSE = (2.15852e-23)
devMSE = (4.47342e-24)

ylim = [1.000000E-024,1.000000E-022]

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

oriPSNR = (347.07)
devPSNR = (353.906)

ylim = [345,357.5]

plt.ylabel('PSNR')
plt.xticks(np.arange(0))
plt.ylim(ylim)
plt.yticks(np.arange(min(ylim), max(ylim), 5.0))
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

"""
Bar chart demo with pairs of bars grouped for easy comparison.
"""
import numpy as np
import matplotlib.pyplot as plt


n_groups = 3

fig, ax = plt.subplots(nrows=2)

index = np.arange(n_groups)
bar_width = 0.35

opacity = 0.4
error_config = {'ecolor': '0.3'}

plt.sca(ax[0])

original = (0.179035, 0.327291, 0.18877)
std_ori = (0.00635899, 0.00506325, 0.00703837)

developed = (0.212338, 0.317583, 0.225165)
std_dev = (0.00659189, 0.00760844, 0.00625787)

plt.ylabel('Times')
plt.title('Times of methods for the decimated 1D HWT')
plt.xticks(index + bar_width, ('Decomposition', 'Composition', 'Decomp & Comp'))

rects1 = ax[0].bar(index, original, bar_width,
                 alpha=opacity,
                 color='b',
                 yerr=std_ori,
                 error_kw=error_config,
                 label='Original')

rects2 = ax[0].bar(index + bar_width, developed, bar_width,
                 alpha=opacity,
                 color='r',
                 yerr=std_dev,
                 error_kw=error_config,
                 label='Developed')

plt.legend(loc='upper left')

plt.sca(ax[1])

originalError = (9.822543E-011, 0.000000, 3.028617E-010)
developedError = (7.275958E-012, 0.000000, 3.819878E-011)

plt.ylabel('Error')
plt.title('Calculation errors of the decimated 1D HWT')
plt.xticks(index + bar_width, ('Decomposition', 'Composition', 'Decomp & Comp'))
plt.ylim([0,4.000000E-010])

rectError1 = ax[1].bar(index, originalError, bar_width,
                 alpha=opacity,
                 color='b',
                 #yerr=std_ori,
                 error_kw=error_config,
                 label='Original')

rectError2 = ax[1].bar(index, developedError, bar_width,
                 alpha=opacity,
                 color='r',
                 #yerr=std_ori,
                 error_kw=error_config,
                 label='Developed')

plt.tight_layout()
plt.show()

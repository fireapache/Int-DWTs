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

original = (0.00154292, 0.000443668, 0.00200258)
std_ori = (.000123006, 1.57083e-005, 0.000140022)

developed = (0.00195723, 0.000456063, 0.00264747)
std_dev = (7.71919e-005, 7.54541e-005, 0.00034877)

plt.ylabel('Times')
plt.title('Times of methods for the decimated 1D HWT')
plt.xticks(index + bar_width, ('Decomposition', 'Composition', 'Decomp & Comp'))
plt.ylim([0,0.0035])

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

originalError = (1.891749E-010, 0.000000, 2.124580E-009)
developedError = (5.820766E-011, 0.000000, 4.074536E-010)

plt.ylabel('Error')
plt.title('Calculation errors of the decimated 1D HWT')
plt.xticks(index + bar_width, ('Decomposition', 'Composition', 'Decomp & Comp'))
plt.ylim([0,7.000000E-008])

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

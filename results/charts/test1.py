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

original = (0.00109071, 0.00125803, 0.00236344)
std_ori = (3.71057e-005, 4.4377e-005, 0.000101865)

developed = (0.000805147, 0.00079437, 0.00157319)
std_dev = (4.34346e-005, 2.91344e-005, 5.13007e-005)

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

originalError = (3.128662E-010, 6.332994E-008, 3.405148E-009)
developedError = (1.136868E-012, 5.587935E-009, 2.110028E-010)

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

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

original = (0.0498059, 0.0558483, 0.105519)
std_ori = (0.00165441, 0.00158507, 0.00226186)

developed = (0.0443343, 0.0341321, 0.0796046)
std_dev = (0.00203959, 0.002759, 0.00450897)

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

originalError = (6.082246E-012, 1.443550E-008, 2.519300E-010)
developedError = (1.065814E-014, 9.313226E-010, 4.320100E-012)

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

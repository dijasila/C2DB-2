import matplotlib.pyplot as plt
import matplotlib.cbook as cbook
import matplotlib.image as mpimg
import os
import sys
p = os.path.abspath('/home/niflheim/fafb/papers/c2db2/c2db-version-2/')

if p not in sys.path:
    sys.path.append(p)

from rcparams import rcParams, columnwidth

img = mpimg.imread('dimensionality_structures.png')

fig, ax = plt.subplots()
figsize = (columnwidth, columnwidth*75/79.502)
fig.set_size_inches(*figsize)

im = ax.imshow(img)

plt.text(240, 910, '$N_{2D}^{initial} = 1$', verticalalignment='center', horizontalalignment='center')
plt.text(240, 1010, '$N_{2D}^{final} = 2$', verticalalignment='center', horizontalalignment='center')
plt.text(720, 910, '$N_{2D}^{initial} = 1$', verticalalignment='center', horizontalalignment='center')
plt.text(720, 1010, '$N_{2D}^{final} = 1$', verticalalignment='center', horizontalalignment='center')

plt.axis('off')
plt.tight_layout()
plt.savefig('dimensionality.png', dpi=300)
plt.show()

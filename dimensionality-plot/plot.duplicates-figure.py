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

plt.text(750, 180, '$N_{2D}^{initial} = 1$')
plt.text(750, 280, '$N_{2D}^{final} = 1$')
plt.text(750, 600, '$N_{2D}^{initial} = 1$')
plt.text(750, 700, '$N_{2D}^{final} = 2$')

plt.axis('off')
plt.savefig('dimensionality.png', dpi=300)

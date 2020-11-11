import matplotlib
matplotlib.use('TkAgg')
from gpaw import GPAW
from ase.dft.bandgap import bandgap
import matplotlib.pyplot as plt
import numpy
from pathlib import Path

def get_edge():
    calc = GPAW('../../defects.pristine_sc/gs.gpw')
    gap, p1, p2 = bandgap(calc)
    evbm = calc.get_eigenvalues(spin=p1[0], kpt=p1[1])[p1[2]]
    ecbm = calc.get_eigenvalues(spin=p2[0], kpt=p2[1])[p2[2]]
    return evbm, ecbm, gap

def band_edge(energy, edge, color, offset=2, ax=None):
    if edge == 'vbm':
        eoffset = energy - offset
        elabel = energy - offset/2
    elif edge == 'cbm':
        eoffset = energy + offset
        elabel = energy + offset/2

    ax.plot([0,1],[energy]*2, color=color, lw=2,zorder=1)
    ax.fill_between([0,1],[energy]*2,[eoffset]*2, color=color, alpha=0.7)
    ax.text(0.5, elabel, edge.upper(), color='w', fontsize=18, ha='center', va='center')

class Level:
    " Class to draw a single defect state level in the gap"

    def __init__(self, energy, size=0.05, ax=None):
        self.size = size
        self.energy = energy
        self.ax = ax

    def draw(self, spin, deg):
        """ Method to draw the defect state according to the 
          spin and degeneracy"""

        relpos = [[1/4,1/8],[3/4,5/8]][spin][deg-1]
        pos = [relpos - self.size, relpos + self.size]
        self.relpos = relpos
        self.spin = spin
        self.deg = deg

        if deg == 1:
            self.ax.plot(pos, [self.energy] * 2, '-k')

        if deg == 2:
            newpos = [p + 1/4 for p in pos]
            self.ax.plot(pos, [self.energy] * 2, '-k')
            self.ax.plot(newpos, [self.energy] * 2, '-k')

    def add_occupation(self, length):
        " Draw an arrow if the defect state if occupied"

        updown = [1,-1][self.spin]
        self.ax.arrow(self.relpos, self.energy - updown*length/2, 0, updown*length, head_width=0.01, head_length=length/5, fc='k', ec='k')
        if self.deg == 2:
            self.ax.arrow(self.relpos + 1/4, self.energy - updown*length/2, 0, updown*length, head_width=0.01, head_length=length/5, fc='k', ec='k')

fig, ax = plt.subplots()

evbm, ecbm, gap = get_edge()
print(evbm, ecbm, gap)
# Draw bands edge
band_edge(evbm, 'vbm', 'C0', offset=gap/5, ax=ax)
band_edge(ecbm, 'cbm', 'C1', offset=gap/5, ax=ax)
# Loop over eigenvalues to draw the level
calc = GPAW('gs.gpw')
nband = calc.get_number_of_bands()
ef = calc.get_fermi_level()

for s in range(calc.get_number_of_spins()):
    for n in range(nband):
        ene = calc.get_eigenvalues(spin=s, kpt=0)[n]
        occ = calc.get_occupation_numbers(spin=s, kpt=0)[n]
        enenew = calc.get_eigenvalues(spin=s, kpt=0)[n+1]
        print(n, ene, occ)
        lev = Level(ene, ax=ax)
        if (ene >= evbm + 0.05 and ene <= ecbm - 0.05):
            if abs(enenew - ene) <= 0.03:
                lev.draw(spin=s, deg=2)
            elif abs(eneold - ene) <= 0.03:
                continue
            else:
                lev.draw(spin=s, deg=1)
            if ene <= ef:
                lev.add_occupation(length=gap/10)
        if ene >= ecbm:
            break
        eneold = ene

fontsize_ticks = 16
fontsize_labels = 18

ax.plot([0,1],[ef]*2, '--k')
#ax.plot([0.5]*2,[-1,3], '--k')
#ax.plot([0.25]*2,[-1,3], '--k')
#ax.plot([0.75]*2,[-1,3], '--k')
ax.set_xlim(0,1)
ax.set_ylim(evbm-gap/5,ecbm+gap/5)
ax.set_xticks([])
plt.yticks(fontsize=fontsize_ticks)
ax.set_ylabel('Energy (eV)', fontsize=fontsize_labels)
plt.tight_layout()
#filename = '/home/niflheim/smanti/5-Update/ks-plot/' + calc.atoms.get_chemical_formula()
#filename = calc.atoms.get_chemical_formula()
filename = 'ks'
plt.savefig(filename + '.png', dpi=150)
plt.close()
plt.show()

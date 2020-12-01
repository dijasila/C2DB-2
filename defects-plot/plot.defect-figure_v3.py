import matplotlib.pyplot as plt
from asr.core import read_json
import numpy as np
from pathlib import Path
from gpaw import GPAW
from ase.dft.bandgap import bandgap
from matplotlib.gridspec import GridSpec

import os
import sys
p = os.path.abspath('/home/niflheim/fafb/papers/c2db2/c2db-version-2/')

if p not in sys.path:
    sys.path.append(p)

from rcparams import rcParams, textwidth

def f(x, a, b):
    return a * x + b


def order_transitions(trans_dict):
    translist = []
    reflist = ['0/1', '1/2', '2/3', '-1/0', '-2/-1', '-3/-2']
    for element in trans_dict:
        if element in reflist:
            translist.append(element)

    ordered_dict = {}
    for element in reflist:
        if element in trans_dict:
            ordered_dict[element] = trans_dict[element]

    print(f'Ordered dictionary: {ordered_dict}')

    return ordered_dict


def get_b(x, y, a):
    return y - a * x


def plot_formation_energies(ax):
    defnamelist = ['vC', 'CSi']
    c = ['C2', 'C4']
    data = read_json(f'results-asr.sj_analyze_{defnamelist[0]}.json')
    vbm = data['pristine']['vbm'] - data['pristine']['evac']
    cbm = data['pristine']['cbm'] - data['pristine']['evac']
    gap = cbm - vbm
    eform = data['eform']

    ax.fill_betweenx([-10, 30], vbm - 10, vbm, color='C0', alpha=0.5)
    ax.fill_betweenx([-10, 30], cbm + 10, cbm, color='C1', alpha=0.5)
    ax.axhline(0, color='black', linestyle='dotted')

    ax.set_xlim(vbm - 0.1 * gap, cbm + 0.1 * gap)
    ax.set_ylim(-0.1, eform + 0.2 * eform)

    for l in [0, 1]:
        data = read_json(f'results-asr.sj_analyze_{defnamelist[l]}.json')

        eform = data['eform']

        transitions = data['transitions']
        #ax[l].set_ylim(-0.1, 9.2)
        e_m = transitions["-1/0"][0] - transitions["-1/0"][1] - transitions["-1/0"][2]
        e_p = transitions["0/1"][0] - transitions["0/1"][1] - transitions["0/1"][2]
        ax.plot([vbm, cbm], [eform, eform], color=c[l])
        # if e_m < cbm and e_m > vbm:
        #     ax[l].axvline(e_m, color='black', linestyle='-.')
        # if e_p < cbm and e_p > vbm:
        #     ax[l].axvline(e_p, color='black', linestyle='-.')

        transitions = order_transitions(transitions)

        enlist = []
        for element in transitions:
            enlist.append(transitions[element][0] - transitions[element][1] - transitions[element][2])

        # ax2 = ax[l].twiny()

        tickslist = []
        labellist = []
        for i, element in enumerate(transitions):
            energy = transitions[element][0] - transitions[element][1] - transitions[element][2]
            enlist.append(energy)
            name = element
            if energy > vbm and energy < cbm:
                ax.axvline(energy, linewidth=0.5, linestyle='-.', color='grey')
                if name.split('/')[1].startswith('0') and name.split('/')[0].startswith('-'):
                    y1 = eform
                    y2 = eform
                elif name.split('/')[0].startswith('0'):
                    y1 = eform
                    y2 = eform
                elif not name.split('/')[0].startswith('-'):
                    y2 = None
                else:
                    y1 = None
                if name.split('/')[0].startswith('-'):
                    tickslist.append(energy)
                    labellist.append(name)
                    a = float(name.split('/')[0])
                    b = get_b(enlist[i], y2, a)
                    if y1 is None:
                        y1 = f(enlist[i], a, b)
                    y2 = f(enlist[i + 1], a, b)
                    print(enlist[i], enlist[i+1], y1, y2, a)
                    ax.plot([vbm, cbm], [f(vbm, a, b), f(cbm, a, b)], color=c[l])
                    ax.plot([enlist[i], enlist[i + 1]], [y1, y2], color=c[l], marker='x')
                elif not name.split('/')[0].startswith('-'):
                    tickslist.append(energy)
                    labellist.append(name)
                    a = float(name.split('/')[1])
                    b = get_b(enlist[i], y1, a)
                    if y2 is None:
                        y2 = f(enlist[i], a, b)
                    y1 = f(enlist[i + 1], a, b)
                    print(enlist[i], enlist[i+1], y1, y2, a)
                    ax.plot([vbm, cbm], [f(vbm, a, b), f(cbm, a, b)], color=c[l])
                    ax.plot([enlist[i + 1], enlist[i]], [y1, y2], color=c[l], marker='x')
    ax.set_xlabel('$E_F$ [eV]')
    # ax[l].set_xticks([-4.7, -3.7])
    ax.set_ylim(-0.1, 3.0)
    ax.tick_params(axis='both')
    legendlist = ['$\mathrm{V_{C}}$', '$\mathrm{C_{Si}}$']
    for i, color in enumerate(c):
        ax.plot([], [], label=legendlist[i], color=color)
    ax.legend()
    # ax2.set_xlim(ax[l].get_xlim())
    # ax2.set_xticks(tickslist)
    # ax2.set_xticklabels(labellist)

    # ax[1].set_yticks([])
    # ax[1].set_yticklabels([])
    ax.set_ylabel('Formation energy [eV]')
    ypos = ax.get_ylim()[0] + 0.5 * (ax.get_ylim()[1] - ax.get_ylim()[0])
    xdiff = abs(vbm - ax.get_xlim()[0]) * 0.5
    ax.text(vbm - xdiff, ypos, 'VBM', ha='center', va='center', rotation=90, color='w')
    ax.text(cbm + xdiff, ypos, 'CBM', ha='center', va='center', rotation=90, color='w')
    ax.axvline(vbm, color='C0', lw=2, zorder=1)
    ax.axvline(cbm, color='C1', lw=2, zorder=1)
    # ax[0].text(ax[0].get_xlim()[0] + 0.5 * (ax[0].get_xlim()[1] - ax[0].get_xlim()[0]), 
    #            ax[0].get_ylim()[0] + 0.2 * (ax[0].get_ylim()[1] - ax[0].get_ylim()[0]),
    #            '$\mathrm{V_C}$', ha='center', va='center')
    # ax[1].text(ax[1].get_xlim()[0] + 0.5 * (ax[1].get_xlim()[1] - ax[1].get_xlim()[0]), 
    #            ax[1].get_ylim()[0] + 0.2 * (ax[1].get_ylim()[1] - ax[1].get_ylim()[0]),
    #            '$\mathrm{C_{Si}}$', ha='center', va='center')


def get_edge():
    calc = GPAW('pristine_gs.gpw')
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
    ax.fill_between([0,1],[energy]*2,[eoffset]*2, color=color, alpha=0.5)
    ax.text(0.5, elabel, edge.upper(), color='w', ha='center', va='center')


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


def plot_ks(ax):
    evbm, ecbm, gap = get_edge()
    # Draw bands edge
    band_edge(evbm, 'vbm', 'C0', offset=gap/5, ax=ax)
    band_edge(ecbm, 'cbm', 'C1', offset=gap/5, ax=ax)
    # Loop over eigenvalues to draw the level
    calc = GPAW('defect_gs.gpw')
    nband = calc.get_number_of_bands()
    ef = calc.get_fermi_level()
    for s in range(calc.get_number_of_spins()):
        for n in range(nband):
            ene = calc.get_eigenvalues(spin=s, kpt=0)[n]
            occ = calc.get_occupation_numbers(spin=s, kpt=0)[n]
            enenew = calc.get_eigenvalues(spin=s, kpt=0)[n+1]
            lev = Level(ene, ax=ax)
            if (ene >= evbm - 0.15 and ene <= ecbm + 0.15):
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
    ax.plot([0,1],[ef]*2, '--k')
    #ax.plot([0.5]*2,[-1,3], '--k')
    #ax.plot([0.25]*2,[-1,3], '--k')
    #ax.plot([0.75]*2,[-1,3], '--k')
    ax.set_xlim(0,1)
    ax.set_xticks([], [])
    # ax.set_xticks([0,1])
    # ax.set_xticklabels([])
    ax.set_ylim(evbm-gap/5,ecbm+gap/5)
    #ax.set_yticks([0, 1])
    #ax.set_yticklabels([0,1])
    ax.set_ylabel('Energy [eV]')
    delta = 0.02
    ax.text(0.30 + delta, -2.74, '$a_2$', horizontalalignment='left')
    ax.text(0.198 - delta, -3.14, '$a_1$', horizontalalignment='right')
    ax.text(0.80 + delta, -1.32, '$a_2$', horizontalalignment='left')
    ax.text(0.70 - delta, -1.91, '$a_1$', horizontalalignment='right')
    ax.text(0.05, 0.95, '$\mathrm{V_{Si}}$', horizontalalignment='left', verticalalignment='center')
    # ax.text(-2.7, 1, '$a_1$', verticalalignment='center')
    #filename = '/home/niflheim/smanti/5-Update/ks-plot/' + calc.atoms.get_chemical_formula()
    #filename = calc.atoms.get_chemical_formula()


def parabola(x):
    return x**2

def plot_cc(ax):
    x1 = np.linspace(-1, 1, 100)
    x2 = np.linspace(-0.8, 1.2, 100)
    delta = 0.05
    #arrow_params = {'shape': 'full', 'color': 'black', 'width': 0.005, 'length_includes_head': True, 'head_length': 0.1}
    arrow_params = {'color': 'black', 'length_includes_head': True, 'head_width': 0.03, 'head_length': 0.09,
                    'overhang': 0.3}
    ax.arrow(0.6, 2, -0.6, -2, **arrow_params)
    #ax.arrow(0, 0 + delta, 0.6, 2 - delta, **arrow_params)
    ax.arrow(0, 0, 0, 2.35, **arrow_params)
    ax.arrow(1.2, 2 + delta, 0, 0.35 - delta, **arrow_params)
    ax.arrow(1.2, 2.35 - delta, 0, -0.35 + delta, **arrow_params)
    ax.arrow(-0.6, 0 + delta, 0, 0.35 - delta, **arrow_params)
    ax.arrow(-0.6, 0.35 - delta, 0, -0.35 + delta, **arrow_params)
    ax.text(-0.8, 2.5, '$D_{\mathrm{excited}}$', color='C1')
    ax.text(1.1, 0.5, '$D_{\mathrm{ground}}$', color='C0')
    ax.text(1.2 + delta, 2.1725, '$\lambda_{\mathrm{exc}}^{\mathrm{reorg}}$', verticalalignment='center', horizontalalignment='left')
    ax.text(0.55 + delta, 1.35, '$E_{\mathrm{ZPL}} = 3.84$ eV', verticalalignment='center', horizontalalignment='left')
    ax.text(0.0 - delta, 1.4025, 'Excitation', verticalalignment='center', horizontalalignment='right', rotation=90)
    ax.text(-0.6 - delta, 0.1725, '$\lambda_{\mathrm{gs}}^{\mathrm{reorg}}$', verticalalignment='center', horizontalalignment='right')
    ax.text(-0.76 - delta, 2.90, '$\mathrm{V_{Si}}$', horizontalalignment='right', verticalalignment='center')

    # draw vibrational levels
    for element in np.arange(0.2, 0.8, 0.2):
        ax.plot([-np.sqrt(element), np.sqrt(element)],
                [element, element], color='grey')
        ax.plot([-np.sqrt(element) + 0.6, np.sqrt(element) + 0.6],
                [element + 2, element + 2], color='grey')

    ax.plot(x1, parabola(x1))
    ax.plot(x1 + 0.6, parabola(x1)+2)
    ax.set_xticks([], [])
    ax.set_yticks([], [])
    # ax.set_yticks([0,2])
    # ax.set_yticklabels([0, "$E_{exc}$"])
    ax.set_yticklabels([])
    # plt.yticks([],[])
    ax.set_xlabel('$Q$')
    # plt.ylabel("$E - E_{ground}$", fontsize=font_labels)
    ax.set_ylabel("Energy")
    ax.axhline(0, linestyle='dotted', color='black')
    ax.axhline(2, linestyle='dotted', color='black')


def plot_structures(ax):
    import matplotlib.image as mpimg

    # img = mpimg.imread('supercell.png')
    # imgplot = ax.imshow(img, interpolation='gaussian')
    #circle = plt.Circle((1580, 940), 80, color='r', fill=False)
    #ax.add_artist(circle)
    ax.axis('off')

def append_label(ax, letter='a'):
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    print(xlim, ylim)
    if letter == 'a':
        ax.text(xlim[0] - 0.12*abs(xlim[0] - xlim[1]), ylim[1] + 0.1*abs(ylim[0] - ylim[1]), letter, horizontalalignment='right', verticalalignment='top', weight='bold')
    else:
        ax.text(xlim[0] - 0.12*abs(xlim[0] - xlim[1]), ylim[1] + 0.1*abs(ylim[0] - ylim[1]), letter, horizontalalignment='right', verticalalignment='top', weight='bold')


# fig, axs = plt.subplots(ncols=4, nrows=4)
#
# # set subplot for structure visualisation
# gs = axs[0, 0].get_gridspec()
# for ax in axs[0:2, 0]:
#     ax.remove()
# for ax in axs[0:2, 1]:
#     ax.remove()
# axstruc = fig.add_subplot(gs[:2, :2])
# 
# # set subplot for KS state
# gs = axs[2, 0].get_gridspec()
# for ax in axs[2:, 0]:
#     ax.remove()
# for ax in axs[2:, 1]:
#     ax.remove()
# axks = fig.add_subplot(gs[2:, 2:])
# 
# # set subplot for exc plot
# gs = axs[0, 2].get_gridspec()
# for ax in axs[2:, 2]:
#     ax.remove()
# for ax in axs[2:, 3]:
#     ax.remove()
# axexc = fig.add_subplot(gs[2:, :2])
# 
# # set subplot for formation plot
# gs = axs[0, 3].get_gridspec()
# for ax in axs[0:2, 2]:
#     ax.remove()
# for ax in axs[0:2, 3]:
#     ax.remove()
# axform1 = fig.add_subplot(gs[0:2, 2:3])
# axform2 = fig.add_subplot(gs[0:2, 3:4])

fig = plt.figure(constrained_layout=True)
figsize = (textwidth, 5.)
fig.set_size_inches(*figsize)
gs = GridSpec(6, 8, figure=fig)
# axform1 = fig.add_subplot(gs[:3, 4:6])
# axform2 = fig.add_subplot(gs[:3, 6:])
axform = fig.add_subplot(gs[:3, 4:])
axexc = fig.add_subplot(gs[3:, :4])
axks = fig.add_subplot(gs[3:, 4:])
axim = fig.add_subplot(gs[:3, :4])
# axform_help = fig.add_subplot(gs[:3, 4:])
# axform_help.axis('off')


plot_formation_energies(axform)
plot_ks(axks)
plot_cc(axexc)
plot_structures(axim)

# labelsp
append_label(axim, 'a')
append_label(axform, 'b')
append_label(axexc, 'c')
append_label(axks, 'd')

# plt.tight_layout()
plt.savefig('defects_v3.pdf')
plt.savefig('defects_v3.png', dpi=300)
plt.show()


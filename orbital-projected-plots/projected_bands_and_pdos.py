"""Plot orbital projected band structure along with pdos for MoS2."""

# Import functionality for plot methodology
import numpy as np
from collections import defaultdict

import matplotlib.patheffects as path_effects
from matplotlib.lines import Line2D

from ase.spectrum.band_structure import BandStructure, BandStructurePlot
from asr.projected_bandstructure import (get_yl_ordering, get_bs_sampling,
                                         get_pie_markers, get_pie_slice)
from asr.pdos import get_ordered_syl_dict, get_yl_colors

# Import c2db-version-2 functionality
# Add c2db-version-2 root dir to path
import os
import sys
p = os.path.abspath('../')
if p not in sys.path:
    sys.path.append(p)

from rcparams import plotter


# ---------- Modified asr plot methodology ---------- #


@plotter()
def projected_bs_pbe(row, ax,
                     npoints=40, markersize=36., res=64):
    """Produce the projected band structure on given matplotlib axis.

    Plot the projection weight fractions as pie charts.

    Parameters
    ----------
    npoints : int,
        number of pie charts per band
    markersize : float
        size of pie chart markers
    res : int
        resolution of the pie chart markers (points around the circumference)
    """

    # Extract projections data
    data = row.data.get('results-asr.projected_bandstructure.json')
    weight_skni = data['weight_skni']
    yl_i = data['yl_i']

    # Get color indeces
    c_i = get_yl_ordering(yl_i, data['symbols'])

    # Extract band structure data
    d = row.data.get('results-asr.bandstructure.json')
    path = d['bs_nosoc']['path']
    ef = d['bs_nosoc']['efermi']

    # If a vacuum energy is available, use it as a reference
    ref = row.get('evac', d.get('bs_nosoc').get('efermi'))
    if row.get('evac') is not None:
        label = r'$E - E_\mathrm{vac}$ [eV]'
    else:
        label = r'$E - E_\mathrm{F}$ [eV]'

    # Determine plotting window based on band gap
    gaps = row.data.get('results-asr.gs.json', {}).get('gaps_nosoc', {})
    if gaps.get('vbm'):
        emin = gaps.get('vbm') - 3
    else:
        emin = ef - 3
    if gaps.get('cbm'):
        emax = gaps.get('cbm') + 3
    else:
        emax = ef + 3

    # Take bands with energies in range
    e_skn = d['bs_nosoc']['energies']
    inrange_skn = np.logical_and(e_skn > emin, e_skn < emax)
    inrange_n = np.any(np.any(inrange_skn, axis=1), axis=0)
    e_skn = e_skn[:, :, inrange_n]
    weight_skni = weight_skni[:, :, inrange_n, :]

    # Use band structure objects to plot outline
    bs = BandStructure(path, e_skn - ref, ef - ref)
    # Use colors if spin-polarized
    if e_skn.shape[0] == 2:
        spincolors = ['0.8', '0.4']
    else:
        spincolors = ['0.8'] * e_skn.shape[0]
    style = dict(
        colors=spincolors,
        ls='-',
        zorder=0)
    bsp = BandStructurePlot(bs)
    bsp.plot(ax=ax, show=False, emin=emin - ref, emax=emax - ref,
             ylabel=label, **style)

    xcoords, k_x = get_bs_sampling(bsp, npoints=npoints)

    # Generate energy and weight arrays based on band structure sampling
    ns, nk, nb = e_skn.shape
    s_u = np.array([s for s in range(ns) for n in range(nb)])
    n_u = np.array([n for s in range(ns) for n in range(nb)])
    e_ux = e_skn[s_u[:, np.newaxis],
                 k_x[np.newaxis, :],
                 n_u[:, np.newaxis]] - ref
    weight_uxi = weight_skni[s_u[:, np.newaxis],
                             k_x[np.newaxis, :],
                             n_u[:, np.newaxis], :]
    # Plot projections
    for e_x, weight_xi in zip(e_ux, weight_uxi):

        # Weights as pie chart
        pie_xi = get_pie_markers(weight_xi, s=markersize,
                                 scale_marker=False, res=res)
        for x, e, weight_i, pie_i in zip(xcoords, e_x, weight_xi, pie_xi):
            # totweight = np.sum(weight_i)
            for i, pie in enumerate(pie_i):
                ax.scatter(x, e, facecolor='C{}'.format(c_i[i]),
                           zorder=3, **pie)

    # Set legend
    # Get "pac-man" style pie slice marker
    pie = get_pie_slice(1. * np.pi / 4.,
                        3. * np.pi / 2., s=markersize, res=res)
    # Generate markers for legend
    legend_markers = []
    for i, yl in enumerate(yl_i):
        legend_markers.append(Line2D([0], [0],
                                     mfc='C{}'.format(c_i[i]), mew=0.0,
                                     marker=pie['marker'], ms=3. * np.pi,
                                     linewidth=0.0))
    # Generate legend
    ax.legend(legend_markers, [yl.replace(',', ' (') + ')' for yl in yl_i],
              bbox_to_anchor=(0., 1.02, 1., 0.), loc='lower left',
              ncol=2, mode="expand", borderaxespad=0.)

    ax.set_xlabel(r'$k$-points')
    xlim = ax.get_xlim()
    x0 = xlim[1] * 0.01
    text = ax.annotate(
        r'$E_\mathrm{F}$',
        xy=(x0, ef - ref),
        ha='left',
        va='bottom')

    text.set_path_effects([
        path_effects.Stroke(linewidth=2, foreground='white', alpha=0.5),
        path_effects.Normal()
    ])


@plotter()
def plot_pdos(row, ax, soc=True):
    """Produce the pdos plot."""

    def smooth(y, npts=3):
        return np.convolve(y, np.ones(npts) / npts, mode='same')

    # Check if pdos data is stored in row
    results = 'results-asr.pdos.json'
    pdos = 'pdos_soc' if soc else 'pdos_nosoc'
    if results in row.data and pdos in row.data[results]:
        data = row.data[results][pdos]
    else:
        return

    # Extract raw data
    symbols = data['symbols']
    pdos_syl = get_ordered_syl_dict(data['pdos_syl'], symbols)
    e_e = data['energies'].copy() - row.get('evac', 0)
    ef = data['efermi']

    # Find energy range to plot in
    if soc:
        emin = row.get('vbm', ef) - 3 - row.get('evac', 0)
        emax = row.get('cbm', ef) + 3 - row.get('evac', 0)
    else:
        nosoc_data = row.data['results-asr.gs.json']['gaps_nosoc']
        vbmnosoc = nosoc_data.get('vbm', ef)
        cbmnosoc = nosoc_data.get('cbm', ef)

        if vbmnosoc is None:
            vbmnosoc = ef

        if cbmnosoc is None:
            cbmnosoc = ef

        emin = vbmnosoc - 3 - row.get('evac', 0)
        emax = cbmnosoc + 3 - row.get('evac', 0)

    # Set up energy range to plot in
    i1, i2 = abs(e_e - emin).argmin(), abs(e_e - emax).argmin()

    # Get color code
    color_yl = get_yl_colors(pdos_syl)

    # Figure out if pdos has been calculated for more than one spin channel
    spinpol = False
    for k in pdos_syl.keys():
        if int(k[0]) == 1:
            spinpol = True
            break

    # Plot pdos
    pdosint_s = defaultdict(float)
    for key in pdos_syl:
        pdos = pdos_syl[key]
        spin, symbol, lstr = key.split(',')
        spin = int(spin)
        sign = 1 if spin == 0 else -1

        # Integrate pdos to find suiting pdos range
        pdosint_s[spin] += np.trapz(y=pdos[i1:i2], x=e_e[i1:i2])

        # Label atomic symbol and angular momentum
        if spin == 0:
            label = '{} ({})'.format(symbol, lstr)
        else:
            label = None

        ax.plot(smooth(pdos) * sign, e_e,
                label=label, color=color_yl[key[2:]])

    ax.axhline(ef - row.get('evac', 0), color='k', ls=':')

    # Set up axis limits
    ax.set_ylim(emin, emax)
    if spinpol:  # Use symmetric limits
        xmax = max(pdosint_s.values())
        ax.set_xlim(-xmax * 0.5, xmax * 0.5)
    else:
        ax.set_xlim(0, pdosint_s[0] * 0.5)

    # Annotate E_F
    xlim = ax.get_xlim()
    x0 = xlim[0] + (xlim[1] - xlim[0]) * 0.01
    text = ax.text(x0, ef - row.get('evac', 0),
                   r'$E_\mathrm{F}$',
                   ha='left',
                   va='bottom')

    text.set_path_effects([
        path_effects.Stroke(linewidth=3, foreground='white', alpha=0.5),
        path_effects.Normal()
    ])

    ax.set_xlabel('projected dos [states / eV]')
    if row.get('evac') is not None:
        ax.set_ylabel(r'$E-E_\mathrm{vac}$ [eV]')
    else:
        ax.set_ylabel(r'$E$ [eV]')

    # Set up legend
    ax.legend(bbox_to_anchor=(0., 1.02, 1., 0.), loc='lower left',
              ncol=2, mode="expand", borderaxespad=0.)


if __name__ == '__main__':
    # General modules
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    from rcparams import rcParams, textwidth

    # Script modules
    from ase.db import connect
    from asr.core.results import dct_to_result, UnknownDataFormat

    # --------------- Inputs --------------- #

    # Data
    db = connect("../pareto-plot/mos2/mos2.db")

    # Figure
    figsize = (textwidth, 4.)
    npoints = 17
    res = 128
    filename = "../../MoS2-H_pbands-pdos.pdf"
    dpi = 300

    # --------------- Script --------------- #

    # Set up figure and axes
    fig = plt.figure(1)
    fig.set_size_inches(*figsize)
    gspec = gridspec.GridSpec(1, 2)
    axes = [fig.add_subplot(ss) for ss in gspec]

    for row in db.select():
        # Hack ase row to wrap data in asr objects
        for key, item in row.data.items():
            try:
                new_item = dct_to_result(item)
            except UnknownDataFormat:
                # Keep data as is
                new_item = item
            row._data[key] = new_item

        projected_bs_pbe(row, axes[0], npoints=npoints, res=res)
        plot_pdos(row, axes[1], soc=False)

        break  # Only one MoS2 row

    gspec.tight_layout(fig)
    fig.savefig(filename, dpi=dpi, format=filename.split('.')[-1])
    plt.show()

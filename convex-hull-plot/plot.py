import functools
from ase.phasediagram import PhaseDiagram
from asr.convex_hull import convex_hull_tables, Result

from asr.core.material import (get_material_from_folder,
                               make_panel_figures)
from asr.core import read_json
import os
import sys
import matplotlib.pyplot as plt

plt.rc('text', usetex=True)
plt.rc('text.latex',
       preamble=r'\usepackage{xcolor}')

p = os.path.abspath('../')
if p not in sys.path:
    sys.path.append(p)

from rcparams import plotter, textwidth, columnwidth
from ase.formula import Formula


def get_hull_energies(pd: PhaseDiagram):
    hull_energies = []
    for ref in pd.references:
        count = ref[0]
        refenergy = ref[1]
        natoms = ref[3]
        decomp_energy, indices, coefs = pd.decompose(**count)
        ehull = (refenergy - decomp_energy) / natoms
        hull_energies.append(ehull)

    return hull_energies


def webpanel(result, row, key_descriptions):
    from asr.database.browser import fig, table, describe_entry

    caption = """
    The convex hull describes stability
    with respect to other phases."""
    hulltable1 = table(row,
                       'Stability',
                       ['hform', 'ehull'],
                       key_descriptions)
    hulltables = convex_hull_tables(row)
    panel = {'title': 'Thermodynamic stability',
             'columns': [[fig('convex-hull.png', caption=caption)],
                         [hulltable1] + hulltables],
             'plot_descriptions': [{'function':
                                    functools.partial(plot, thisrow=row),
                                    'filenames': ['convex-hull.pdf']}],
             'sort': 1}

    thermostab = row.get('thermodynamic_stability_level')
    stabilities = {1: 'low', 2: 'medium', 3: 'high'}
    high = 'Heat of formation < convex hull + 0.2 eV/atom'
    medium = 'Heat of formation < 0.2 eV/atom'
    low = 'Heat of formation > 0.2 eV/atom'
    row = ['Thermodynamic',
           describe_entry(stabilities[thermostab].upper(),
                          '\n'.join([f'LOW: {low}',
                                     f'MEDIUM: {medium}',
                                     f'HIGH: {high}']))]

    summary = {'title': 'Summary',
               'columns': [[{'type': 'table',
                             'header': ['Stability', ''],
                             'rows': [row],
                             'columnwidth': 3}]],
               'sort': 1}
    return [panel, summary]


Result.formats = {"ase_webpanel": webpanel}


@plotter()
def plot(row, fname, thisrow):
    from ase.phasediagram import PhaseDiagram
    import matplotlib.pyplot as plt

    data = row.data['results-asr.convex_hull.json']

    count = row.count_atoms()
    if not (2 <= len(count) <= 3):
        return

    references = data['references']
    pdrefs = []
    legends = []
    colors = []
    sizes = []
    for reference in references:
        h = reference['natoms'] * reference['hform']
        pdrefs.append((reference['formula'], h))
        if reference['legend'] not in legends:
            legends.append(reference['legend'])
        idlegend = legends.index(reference['legend'])
        colors.append(f'C{idlegend + 2}')
        sizes.append((3 * idlegend + 3)**2)

    pd = PhaseDiagram(pdrefs, verbose=False)

    fig = plt.figure(figsize=(columnwidth, columnwidth))
    ax = fig.gca()

    legendhandles = []

    for it, label in enumerate(['On hull', 'off hull']):
        handle = ax.fill_between([], [],
                                 color=f'C{it + 2}', label=label)
        legendhandles.append(handle)

    for it, legend in enumerate(legends):
        handle = ax.scatter([], [], facecolor='none', marker='o',
                            edgecolor='k', label=legend, s=(3 + it * 3)**2)
        legendhandles.append(handle)

    if len(count) == 2:
        x, e, _, hull, simplices, xlabel, ylabel = pd.plot2d2()
        for i, j in simplices:
            ax.plot(x[[i, j]], e[[i, j]], '-', color='C0')
        names = [ref['label'] for ref in references]
        if row.hform < 0:
            mask = e < 0.05
            e = e[mask]
            x = x[mask]
            hull = hull[mask]
            names = [name for name, m in zip(names, mask) if m]
        ax.scatter(x, e, facecolor='none', marker='o', edgecolor=colors)

        delta = e.ptp() / 30
        for a, b, name, on_hull in zip(x, e, names, hull):
            va = 'center'
            ha = 'left'
            dy = 0
            dx = 0.02
            ax.text(a + dx, b + dy, name, ha=ha, va=va)

        A, B = pd.symbols
        ax.set_xlabel('{}$_{{1-x}}${}$_x$'.format(A, B))
        ax.set_ylabel(r'$\Delta H$ [eV/atom]')

        # Circle this material
        xt = count.get(B, 0) / sum(count.values())
        # ax.plot([xt], [row.hform], 'o', color='C1', label=f'{thisrow.formula}')
        ymin = e.min()

        ax.axis(xmin=-0.1, xmax=1.1, ymin=ymin - 2.5 * delta)
        plt.legend(loc='lower left')
    else:
        x, y, _, hull, simplices = pd.plot2d3()
        import numpy as np
        hull = np.array(hull)
        hull_energies = get_hull_energies(pd)
        hull = np.array(hull_energies) < 0.05
        names = [ref['label'] for ref in references]
        latexnames = [format(Formula(name.split(' ')[0]).reduce()[0], 'latex') for name in names]
        for i, j, k in simplices:
            ax.plot(x[[i, j, k, i]], y[[i, j, k, i]], '-', color='lightblue')
        edgecolors = ['C2' if hull_energy < 0.05 else 'C3'
                      for hull_energy in hull_energies]
        ax.scatter(
            x[~hull], y[~hull],
            facecolor='none', marker='o',
            edgecolor=np.array(edgecolors)[~hull], s=np.array(sizes)[~hull],
            zorder=9,
        )

        ax.scatter(
            x[hull], y[hull],
            facecolor='none', marker='o',
            edgecolor=np.array(edgecolors)[hull], s=np.array(sizes)[hull],
            zorder=10,
        )

        printed_names = set()
        for a, b, name, on_hull, hull_energy in zip(x, y, latexnames, hull, hull_energies):
            print(name, on_hull)
            if name in [
                    'BiITe', 'I', 'Bi',
                    'Te', 'Bi$_{2}$Te$_{3}$',
                    'BiI$_{3}$',
            ] and name not in printed_names:
                printed_names.add(name)
                if name == 'Bi$_{2}$Te$_{3}$':
                    ax.text(a, b - 0.03, name, ha='center', va='top')
                else:
                    ax.text(a - 0.02, b, name, ha='right', va='top')
        A, B, C = pd.symbols
        bfrac = count.get(B, 0) / sum(count.values())
        cfrac = count.get(C, 0) / sum(count.values())

        from matplotlib.legend_handler import HandlerPatch
        from matplotlib import patches

        class ObjectHandler:
            def legend_artist(self, legend, orig_handle, fontsize, handlebox):
                x0, y0 = handlebox.xdescent, handlebox.ydescent
                width, height = handlebox.width, handlebox.height
                patch = patches.Polygon(
                    [
                        [x0, y0],
                        [x0, y0 + height],
                        [x0 + 3 / 4 * width, y0 + height],
                        [x0 + 1 / 4 * width, y0],
                    ],
                    closed=True, facecolor='C2',
                    edgecolor='none', lw=3,
                    transform=handlebox.get_transform())
                handlebox.add_artist(patch)
                patch = patches.Polygon(
                    [
                        [x0 + width, y0],
                        [x0 + 1 / 4 * width, y0],
                        [x0 + 3 / 4 * width, y0 + height],
                        [x0 + width, y0 + height],
                    ],
                    closed=True, facecolor='C3',
                    edgecolor='none', lw=3,
                    transform=handlebox.get_transform())
                handlebox.add_artist(patch)
                # artist, = plt.plot(*zip([x0 + 1 / 4 * width, y0],
                #                         [x0 + 3 / 4 * width, y0 + height]),
                #                    color='k',
                #                    transform=handlebox.get_transform())
                # handlebox.add_artist(artist)
                return patch


        from matplotlib.legend_handler import HandlerLine2D, HandlerTuple

        newlegendhandles = [(legendhandles[0], legendhandles[1]),
                            legendhandles[2], legendhandles[3]]
        plt.legend(
            newlegendhandles,
            [r'${^</_>}\, 5 \mathrm{meV}$',
             legends[0], legends[1]], loc='upper right', handletextpad=0.5,
            handler_map={tuple: ObjectHandler()},
            # bbox_to_anchor=(0.4, 1),
        )
        plt.axis('off')

    plt.tight_layout()
    plt.savefig(fname)
    plt.close()


res = read_json('results-asr.convex_hull.json')

material = get_material_from_folder('.')
panels = res.format_as('ase_webpanel', material, {})
make_panel_figures(material, panels)

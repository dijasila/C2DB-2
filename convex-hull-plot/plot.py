import functools
from asr.convex_hull import convex_hull_tables, Result

from asr.core.material import (get_material_from_folder,
                               make_panel_figures)
from asr.core import prepare_result, ASRResult, read_json
import typing
import os
import sys

p = os.path.abspath('../')
if p not in sys.path:
    sys.path.append(p)

from rcparams import plotter, textwidth, columnwidth
from ase.formula import Formula


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

    for it, legend in enumerate(legends):
        ax.scatter([], [], facecolor='none', marker='o',
                   edgecolor=f'C{it + 2}', label=legend, s=(3 + it * 3)**2)

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
        
        names = [ref['label'] for ref in references]
        latexnames = [format(Formula(name.split(' ')[0]).reduce()[0], 'latex') for name in names]
        for i, j, k in simplices:
            ax.plot(x[[i, j, k, i]], y[[i, j, k, i]], '-', color='lightblue')
        ax.scatter(x, y, facecolor='none', marker='o', edgecolor=colors, s=sizes,
                   zorder=10)

        for a, b, name, on_hull in zip(x, y, latexnames, hull):
            if name in ['BiITe', 'I', 'Bi', 'Te']:
                ax.text(a - 0.02, b, name, ha='right', va='top')
        A, B, C = pd.symbols
        bfrac = count.get(B, 0) / sum(count.values())
        cfrac = count.get(C, 0) / sum(count.values())

        # ax.plot([bfrac + cfrac / 2],
        #         [cfrac * 3**0.5 / 2],
        #         'o', color='C1', label=f'{thisrow.formula}')
        plt.legend(loc='upper right', handletextpad=0)
        plt.axis('off')

    plt.tight_layout()
    plt.savefig(fname)
    plt.close()


res = read_json('results-asr.convex_hull.json')

material = get_material_from_folder('.')
panels = res.format_as('ase_webpanel', material, {})
make_panel_figures(material, panels)

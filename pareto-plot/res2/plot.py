from asr.core import read_json
import numpy as np
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from matplotlib import pyplot as plt
import matplotlib.patheffects as path_effects
from ase.db import connect

import os
import sys

p = os.path.abspath('../../')
if p not in sys.path:
    sys.path.append(p)

from utils import (
    cut_square_sheet,
    make_image_of_2D_material_from_multiple_perspectives
)
from rcparams import plotter, textwidth, columnwidth
# import PIL


def make_pareto_lines(x_p, y_p):
    additional_points = []
    points_pv = sorted(zip(x_p, y_p), key=lambda x: x[0])
    x_p = [point_v[0] for point_v in points_pv]
    y_p = [point_v[1] for point_v in points_pv]
    paddedx_p = x_p + [20]
    all_x = [x_p[0]]
    all_y = [0]

    for x1, x2, y in zip(paddedx_p[:-1], paddedx_p[1:], y_p):
        all_x.extend([x1, x2])
        all_y.extend([y, y])

    plt.plot(all_x, all_y, color='C3', ls='--', label='Pareto front')
    plt.fill_between(all_x, all_y, 1, facecolor='none',
                     edgecolor=[0.7, 0.7, 0.7], hatch='XXX',
                     linewidth=0.0, zorder=-1)


@plotter()
def plot():
    """Plot atomic structures for RMSD figure."""
    rmsdresult = read_json('results-asr.database.rmsd.json')
    rmsd_by_id = rmsdresult['rmsd_by_id']
    uids = sorted(list(set(rmsd_by_id)))
    distances = [
        [rmsd_by_id[uid][uid1]
         if uid1 != uid
         else 0 for uid1 in uids]
        for uid in uids]
    distance_matrix = np.array(
        distances,
        dtype=float
    )
    isnan = np.where(np.isnan(distance_matrix))
    distance_matrix[isnan] = 10
    condensed_distance_matrix = squareform(distance_matrix)
    linkage_matrix = linkage(condensed_distance_matrix, optimal_ordering=False)

    # plt.figure()
    # dendrogram(linkage_matrix, labels=uids)
    clusters = fcluster(linkage_matrix, 0.5, criterion='distance')

    db = connect('res2.db')
    rows = {row.uid: row for row in db.select()}
    treated_cluster_id = 1
    pareto_data_x = []
    pareto_data_y = []

    cluster_uids = [uid
                    for cluster_id, uid in zip(clusters, uids)
                    if cluster_id == treated_cluster_id]

    for uid in cluster_uids:
        row = rows[uid]
        pareto_data_x.append(row.natoms)
        pareto_data_y.append(row.hform)
        atoms = row.toatoms()
        # make_image_of_2D_material_from_multiple_perspectives(atoms)

    plt.figure(figsize=(textwidth, columnwidth))
    make_pareto_lines(pareto_data_x, pareto_data_y)
    plt.scatter(pareto_data_x, pareto_data_y, color='k', zorder=2,
                label='Included materials')
    plt.scatter([6, 9], [-0.3, -0.4], color='C3', zorder=2,
                label='Excluded materials')
    names = ["T'", "T''", "T"]

    for x, y, name in zip(pareto_data_x, pareto_data_y, names):
        plt.annotate(name, xy=(x, y - 0.02), va='top', ha='center')
    # text = plt.annotate('Excluded region', xy=(9, -0.225),
    #                     ha='center', va='top')

    # text.set_path_effects([path_effects.Stroke(linewidth=4,
    #                                            foreground='white'),
    #                        path_effects.Normal()])
    filenames = [
        'ReS2-perspectives.png',
        'Re2S4-perspectives.png',
        'Re4S8-perspectives.png',
    ]
    size = [2.5, 0.25]
    padding = [0.25, 0.01]

    right_upper_origins = [
        [3 - padding[0], -0.225 - padding[1]],
        [6 - padding[0], -0.5 - padding[1]],
        [12 - padding[0], -0.525 - padding[1]],
    ]
    ax = plt.gca()
    for filename, origin in zip(filenames, right_upper_origins):
        bbox = [origin[0] - size[0], origin[0], origin[1] - size[1], origin[1]]
        img = plt.imread(filename)
        ax.imshow(img, aspect='auto', extent=bbox, interpolation='gaussian')

    # plt.annotate('Pareto optimal', xy=(8, -0.525), ha='center', va='top')
    plt.xlabel('Number of atoms per unit cell')
    plt.ylabel(r'$\Delta H_\mathrm{form}$ [eV/atom]')
    plt.xlim(0, 13)
    plt.xticks(range(0, 14))
    plt.legend()  # facecolor='w', framealpha=1)
    plt.ylim(-0.8, -0.19)
    plt.tight_layout()
    plt.savefig('pareto.pdf', dpi=300)

    # pareto_image = PIL.Image.open('pareto.png')
    # pareto_width, pareto_height = pareto_image.size

    # filenames = [
    #     'ReS2-perspectives.png',
    #     'Re2S4-perspectives.png',
    #     'Re4S8-perspectives.png',
    # ]
    # images = []
    # for created_filename in filenames:
    #     image = PIL.Image.open(created_filename)
    #     image_width, image_height = image.size
    #     factor = (pareto_width / 3) / image_width
    #     image = image.resize((
    #         int(image_width * factor),
    #         int(image_height * factor)
    #     ))
    #     images.append(image)

    # structure_width, structure_height = images[0].size

    # total_width = pareto_width
    # total_height = pareto_height + structure_height

    # new_image = PIL.Image.new('RGBA', (total_width, total_height))
    # new_image.paste(pareto_image, box=(0, 0))

    # for i, image in enumerate(images):
    #     box = (
    #         i * structure_width,
    #         pareto_height,
    #         i * structure_width + image.size[0],
    #         pareto_height + image.size[1],
    #     )
    #     new_image.paste(
    #         image,
    #         box=box,
    #     )

    # new_image.save('combined-pareto.png')
    plt.show()


if __name__ == '__main__':
    plot()

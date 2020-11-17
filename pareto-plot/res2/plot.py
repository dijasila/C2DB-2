from asr.core import read_json
import numpy as np
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from matplotlib import pyplot as plt
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


def make_pareto_lines(x_p, y_p):
    points_pv = sorted(zip(x_p, y_p), key=lambda x: x[0])
    x_p = [point_v[0] for point_v in points_pv]
    y_p = [point_v[1] for point_v in points_pv]
    paddedx_p = x_p + [20]
    paddedy_p = [0] + y_p

    for x1, x2, y in zip(paddedx_p[:-1], paddedx_p[1:], y_p):
        plt.plot([x1, x2], [y, y], color='k')

    for y1, y2, x in zip(paddedy_p[:-1], paddedy_p[1:], x_p):
        plt.plot([x, x], [y1, y2], color='k')


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
    plt.scatter(pareto_data_x, pareto_data_y)
    plt.xlabel('Number of atoms')
    plt.ylabel('HOF [eV/atom]')
    plt.xlim(0, 13)
    ylim = plt.ylim()
    plt.xticks(range(0, 14))
    make_pareto_lines(pareto_data_x, pareto_data_y)
    plt.ylim(ylim)
    plt.tight_layout()
    plt.savefig('pareto.png', dpi=300)
    plt.show()


if __name__ == '__main__':
    plot()

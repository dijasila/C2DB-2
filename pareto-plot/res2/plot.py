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

    plt.figure()
    dendrogram(linkage_matrix, labels=uids)
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
        make_image_of_2D_material_from_multiple_perspectives(atoms)


if __name__ == '__main__':
    plot()

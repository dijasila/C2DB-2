from asr.core import read_json
import numpy as np
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from matplotlib import pyplot as plt
from ase.db import connect
from ase import Atoms

rmsdresult = read_json('results-asr.database.rmsd.json')
rmsd_by_id = rmsdresult['rmsd_by_id']
uids = sorted(list(set(rmsd_by_id)))
print(rmsd_by_id)
distances = [[rmsd_by_id[uid][uid1] if uid1 != uid else 0 for uid1 in uids] for uid in uids]
distance_matrix = np.array(
    distances,
    dtype=float
)
isnan = np.where(np.isnan(distance_matrix))
distance_matrix[isnan] = 10
print('uids=')
print(uids)
print('distance_matrix=')
print(distance_matrix)
condensed_distance_matrix = squareform(distance_matrix)
print('condensed_distance_matrix=')
print(condensed_distance_matrix)
linkage_matrix = linkage(condensed_distance_matrix, optimal_ordering=False)
print('linkage_matrix=')
print(linkage_matrix)

plt.figure()
dendrogram(linkage_matrix, labels=uids)


clusters = fcluster(linkage_matrix, 0.5, criterion='distance')
print('clusters=')
print(clusters)

db = connect('res2.db')
rows = {row.uid: row for row in db.select()}


natoms_in_image = 30


def cut_out_square_sheet(atoms):
    atoms = atoms.repeat((10, 10, 1))
    atoms.center()
    atoms.wrap()
    cutsize = 20
    cut_cv = np.diag([cutsize, ] * 3)
    pos_av = atoms.positions - atoms.positions[0] + [0, 0, cutsize / 2]
    spos_cut_ac = np.dot(pos_av, np.linalg.inv(cut_cv))
    mask = np.logical_and(spos_cut_ac < 1.0,
                          spos_cut_ac > 0.0).all(axis=1)
    new_atoms = Atoms([atom for atom, thismask
                       in zip(atoms, mask) if thismask])
    return new_atoms

# from ase.build.supercells import make_supercell
#     lat = atoms.cell.get_bravais_lattice()
#     eigenvals, eigenvecs = np.linalg.qr(atoms.cell)
#     print(atoms.symbols)
#     print('eigenvals', eigenvals)
#     print('eigenvecs', eigenvecs)
#     # conv_lat = lat.conventional()
#     print(lat)
#     print(lat.conventional_cellmap)
#     conv_atoms = make_supercell(atoms, lat.conventional_cellmap)

#     return conv_atoms



max_cluster_id = max(clusters)
for treated_cluster_id in [1]:  # range(max_cluster_id):
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
        atoms = cut_out_square_sheet(atoms)
        # cell_cv = atoms.cell
        # natoms_in_image;
        atoms.write(
            f'{row.formula}.pov',
            run_povray=True,
            # background='White',
            # transparent=False,
            display=False,
            # rotation='-90x,-30y',
            canvas_width=1000,
            camera_dist=10,
            # radii=0.5,
        )

    print('cluster_id', treated_cluster_id)
    print('')
    plt.figure()
    plt.title(f'cluster_id={treated_cluster_id}')
    plt.scatter(pareto_data_x, pareto_data_y)

# plt.show()

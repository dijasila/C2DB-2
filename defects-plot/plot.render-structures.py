import os
import sys
from ase.io.pov import get_bondpairs

p = os.path.abspath('/home/niflheim/fafb/papers/c2db2/c2db-version-2/')
if p not in sys.path:
    sys.path.append(p)

from utils import cut_square_sheet, make_image_of_2D_material_from_multiple_perspectives, run_povray
from ase.io import read

structures = ['primitive', 'defect_vSi']
kwargs = {'canvas_width': 1000, 'celllinewidth': 0.}
for name in structures:
    atoms = read(name + '.json')
    # atoms = cut_square_sheet(atoms)
    # if name.split('_')[-1] == 'C2':
    #     bondatoms = None
    # elif name.split('_')[-1] == 'Ge2Se2':
    #     bondatoms = get_bondpairs(atoms, radius=1.)
    # top view
    filename = run_povray(atoms, rotation='0x,0y', **kwargs)
    # rename file
    os.system(f'mv {filename} {name}.png')
    # cut view
    # filename = run_povray(atoms, bondatoms=bondatoms)
    # # rename file
    # os.system(f'mv {filename} {name}_cut.png')

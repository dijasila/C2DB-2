"""Plotting utilities for C2DB version2 paper.

In order to import this you have to put something like

p = os.path.abspath('../../')
if p not in sys.path:
    sys.path.append(p)

from utils import (
    cut_square_sheet,
    make_image_of_2D_material_from_multiple_perspectives
)

In the top of your script.

Comments:

26/11-2020:

    If you are experiencing problems with PIL you can try:

        from PIL import Image

    which worked for a problem Thorsten had on openSUSE.

"""

from pathlib import Path
import PIL

import numpy as np
from ase import Atoms
from ase.io.pov import get_bondpairs


def cut_square_sheet(atoms):
    """Cut out a square sheet of 2D material."""
    delta = 0.001
    atoms = atoms.repeat((10, 10, 1))
    atoms.center()
    atoms.wrap()
    cutsize = 20
    cut_cv = np.diag([cutsize, ] * 3)
    pos_av = atoms.positions - atoms.positions[0] + [delta, delta, cutsize / 2]
    spos_cut_ac = np.dot(pos_av, np.linalg.inv(cut_cv))
    mask = np.logical_and(spos_cut_ac < 1.0,
                          spos_cut_ac > 0.0).all(axis=1)
    new_atoms = Atoms([atom.symbol for atom, thismask
                       in zip(atoms, mask) if thismask],
                      positions=pos_av[mask, :], pbc=False)
    return new_atoms


def _clean_povray_files(filename):
    path = Path(filename)
    suffixes = ['.pov', '.ini']
    for suffix in suffixes:
        extra_path = path.with_suffix(suffix)
        if extra_path.is_file():
            extra_path.unlink()


def run_povray(
        atoms,
        filename=None,
        display=False,
        # canvas_width=1000,
        rotation='-90x,-30y',
        radii=0.5,
        bondatoms=None,
        **kwargs,
):
    """Run povray and save to file.

    Filename defaults to {formula}.png.
    """
    if bondatoms is None:
        bondatoms = get_bondpairs(atoms, radius=0.75)

    if filename is None:
        filename = f'{atoms.symbols.formula:metal}.pov'
    print(f'Povray running for {filename}')
    atoms.write(
        filename=filename,
        run_povray=True,
        # background='White',
        # transparent=False,
        display=display,
        rotation=rotation,
        # canvas_width=canvas_width,
        radii=radii,
        bondatoms=bondatoms,
        **kwargs,
    )
    _clean_povray_files(filename)
    return str(Path(filename).with_suffix('.png'))


def make_images_from_multiple_perspectives(
        atoms,
        rotations=['0x,0y', '-90x', '90y'],
        basename=None,
        **kwargs,
):
    """Image atomic structure from x, y, and z direction.

    Parameters
    ----------
    atoms: Atoms
        Atomic structure to make figure of.
    rotations: List[str]
        Povray object rotations. '90x'=rotate 90deg around x axis.
        '90x,90y'=rotate 90deg around x axis and then 90 deg around y axis.
    basename: str
        Base for filename, for example "MoS2" would give files like
        MoS2-'0x,0y'.png. Defaults to {atoms.symbols.formula:metal}
    kwargs: dict
        Key word arguments handed through to ase-povray interface.

    Returns
    -------
    filenames: List[str]
        List of path of created files.

    """
    if basename is None:
        basename = format(atoms.symbols.formula, 'metal')
    filenames = []
    for rotation in rotations:
        rotation_string = rotation.replace(',', '')
        filename = run_povray(
            atoms,
            filename=f'{basename}-{rotation_string}.pov',
            rotation=rotation,
            **kwargs)
        filenames.append(filename)
    return filenames


def make_image_of_2D_material_from_multiple_perspectives(
        atoms,
        filename=None
):
    """Make combined image of 2D mat seen from x, y and z direction.

    Image atomic structure from x, y, and z direction.

    Parameters
    ----------
    atoms: Atoms
        Atomic structure to make figure of.
    filename: str
        Filename for image. Defaults to f'{formula:metal}-perspectives.png'.

    Returns
    -------
    filename: str
        Filename of created file.

    """
    if filename is None:
        formula = atoms.symbols.formula
        filename = f'{formula:metal}-perspectives.png'
    new_atoms = cut_square_sheet(atoms)

    filenames = make_images_from_multiple_perspectives(
        new_atoms,
        rotations=['0x,0y', '-90x', '90y'],
        basename=format(atoms.symbols.formula,
                        'metal'),
    )

    images = []
    for created_filename in filenames:
        images.append(PIL.Image.open(created_filename))

    widths, heights = zip(*(i.size for i in images))
    total_width = widths[0] + widths[2]
    total_height = heights[0] + heights[1]
    new_image = PIL.Image.new('RGBA', (total_width, total_height))
    new_image.paste(images[0], box=(images[2].size[0], 0))
    new_image.paste(images[1], box=(images[2].size[0], images[0].size[1]))
    new_image.paste(images[2], box=(0, 0))
    box = new_image.getbbox()
    new_image = new_image.crop(box=box)
    new_image.save(filename)
    return filename

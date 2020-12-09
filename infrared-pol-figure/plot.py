"""Infrared polarizability."""
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from ase.db import connect
from asr.core import dct_to_result

import os
import sys
p = os.path.abspath('../')
if p not in sys.path:
    sys.path.append(p)

from rcparams import plotter, textwidth, columnwidth
from ase.formula import Formula


@plotter()
def plot():
    """Plot infrared polarizability."""
    db = connect('bn.db')
    row = db.get('BN')
    fnames = ["infrared-pol-x.pdf", "infrared-pol-y.pdf", "infrared-pol-z.pdf"]

    # Get electronic polarizability
    infrareddct = dct_to_result(
        row.data.get("results-asr.infraredpolarizability.json"))
    omega_w = infrareddct["omega_w"] * 1e3
    alpha_wvv = infrareddct["alpha_wvv"]

    electrondct = dct_to_result(
        row.data.get("results-asr.polarizability.json"))
    alphax_w = electrondct["alphax_w"]
    alphay_w = electrondct["alphay_w"]
    alphaz_w = electrondct["alphaz_w"]
    omegatmp_w = electrondct["frequencies"] * 1e3

    # Get max phonon freq
    phonondata = dct_to_result(row.data.get("results-asr.phonons.json"))
    maxphononfreq = phonondata.get("omega_kl")[0].max() * 1e3
    maxomega = maxphononfreq * 3

    atoms = row.toatoms()
    pbc_c = atoms.pbc
    ndim = int(np.sum(pbc_c))

    realphax = interp1d(omegatmp_w, alphax_w.real)
    imalphax = interp1d(omegatmp_w, alphax_w.imag)
    ax_w = (realphax(omega_w) + 1j * imalphax(omega_w)
            + alpha_wvv[:, 0, 0])
    realphay = interp1d(omegatmp_w, alphay_w.real)
    imalphay = interp1d(omegatmp_w, alphay_w.imag)
    ay_w = (realphay(omega_w) + 1j * imalphay(omega_w)
            + alpha_wvv[:, 1, 1])
    realphaz = interp1d(omegatmp_w, alphaz_w.real)
    imalphaz = interp1d(omegatmp_w, alphaz_w.imag)
    az_w = (realphaz(omega_w) + 1j * imalphaz(omega_w)
            + alpha_wvv[:, 2, 2])

    if ndim == 3:
        epsx_w = 1 + 4 * np.pi * ax_w
        epsy_w = 1 + 4 * np.pi * ay_w
        epsz_w = 1 + 4 * np.pi * az_w
        plt.figure()
        plt.plot(omega_w, epsx_w.imag, label='imag')
        plt.plot(omega_w, epsx_w.real, label='real')
        ax = plt.gca()
        # ax.set_title("x-polarization")
        ax.set_xlabel("Energy [meV]")
        ax.set_ylabel(r"Dielectric function")
        ax.set_xlim(0, maxomega)
        ax.legend()
        plt.tight_layout()
        plt.savefig(fnames[0])

        plt.figure()
        plt.plot(omega_w, epsy_w.imag, label='imag')
        plt.plot(omega_w, epsy_w.real, label='real')
        ax = plt.gca()
        # ax.set_title("y-polarization")
        ax.set_xlabel("Energy [meV]")
        ax.set_ylabel(r"Dielectric function")
        ax.set_xlim(0, maxomega)
        ax.legend()
        plt.tight_layout()
        plt.savefig(fnames[1])

        plt.figure()
        plt.plot(omega_w, epsz_w.imag, label='imag')
        plt.plot(omega_w, epsz_w.real, label='real')
        ax = plt.gca()
        # ax.set_title("z-polarization")
        ax.set_xlabel("Energy [meV]")
        ax.set_ylabel(r"Dielectric function")
        ax.set_xlim(0, maxomega)
        ax.legend()
        plt.tight_layout()
        plt.savefig(fnames[2])
    elif ndim in [2, 1, 0]:
        if ndim == 2:
            unit = r"$\mathrm{\AA}$"
        elif ndim == 1:
            unit = r"$\mathrm{\AA}^2$"
        elif ndim == 0:
            unit = r"$\mathrm{\AA}^3$"
        plt.figure(figsize=(columnwidth, columnwidth))
        plt.plot(omega_w, ax_w.imag, label='imag')
        plt.plot(omega_w, ax_w.real, label='real')
        plt.annotate(r'$\alpha_\infty$', xy=(omega_w[-1] - 30, ax_w.real[-1]),
                     va='bottom', ha='center')
        ax = plt.gca()
        # ax.set_title("x-polarization")
        ax.set_xlabel("Energy [meV]")
        ax.set_ylabel(rf"Polarizability [{unit}]")
        ax.set_xlim(0, maxomega)
        ax.legend()
        plt.tight_layout()
        plt.savefig(fnames[0])

        plt.figure(figsize=(columnwidth, columnwidth))
        plt.plot(omega_w, ay_w.imag, label='imag')
        plt.plot(omega_w, ay_w.real, label='real')
        ax = plt.gca()
        # ax.set_title("y-polarization")
        ax.set_xlabel("Energy [meV]")
        ax.set_ylabel(rf"Polarizability [{unit}]")
        ax.set_xlim(0, maxomega)
        ax.legend()
        plt.tight_layout()
        plt.savefig(fnames[1])

        plt.figure(figsize=(columnwidth, columnwidth))
        plt.plot(omega_w, az_w.imag, label='imag')
        plt.plot(omega_w, az_w.real, label='real')
        ax = plt.gca()
        # ax.set_title("z-polarization")
        ax.set_xlabel("Energy [meV]")
        ax.set_ylabel(rf"Polarizability [{unit}]")
        ax.set_xlim(0, maxomega)
        ax.legend()
        plt.tight_layout()
        plt.savefig(fnames[2])


if __name__ == '__main__':
    plot()

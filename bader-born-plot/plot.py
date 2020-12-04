from pathlib import Path
import pickle
import runpy

import numpy as np
from ase.db import connect
import matplotlib.pyplot as plt

width = runpy.run_path('../rcparams.py')['columnwidth']


def extract_data():
    if Path('data.pckl').is_file():
        return
    db = connect('../c2db-data.db')
    data = {}
    for row in db.select('first_class_material'):
        born = row.data.get('results-asr.borncharges.json')
        bader = row.data.get('results-asr.bader.json')
        B_a = bader['kwargs']['data']['bader_charges']
        data[row.uid] = {'B_a': B_a}
        if born is not None:
            Z_avv = born['kwargs']['data']['Z_avv']
            data[row.uid]['Z_avv'] = Z_avv

    with open('data.pckl', 'wb') as fd:
        pickle.dump(data, fd)


def analyze():
    with open('data.pckl', 'rb') as fd:
        data = pickle.load(fd)

    nz = 0
    D = []
    extra = {}
    for u, d in data.items():
        if 'Z_avv' not in d:
            continue
        nz += 1
        B_a = d['B_a']
        Z_avv = d['Z_avv']
        D.extend([(np.trace(Z_vv[:2, :2]) / 2,
                   Z_vv[2, 2],
                   B)
                  for Z_vv, B in zip(Z_avv, B_a)])
        if u.startswith('BN-'):
            symbols = ['B', 'N']
        elif u.startswith('MoS2-'):
            symbols = ['Mo', 'S']
        else:
            symbols = []
        for symbol, Z_vv, B in zip(symbols, Z_avv, B_a):
            assert symbol not in extra
            extra[symbol] = (np.trace(Z_vv[:2, :2]) / 2,
                             Z_vv[2, 2],
                             B)

    zin, zout, b = np.array(D).T
    bad = abs(zin) >= 5
    print(nz, len(b), sum(bad))

    return zin, zout, b, extra


def plot():
    zin, zout, b, extra = analyze()

    fig = plt.figure(figsize=(width, width),
                     constrained_layout=True)
    ax = fig.add_subplot(111)

    ax.plot(b, zin, 'o', label='in-plane')

    x = [-5.0, 5.0]
    p1, p0 = fit = np.polyfit(b, zin, 1)
    y = np.polyval(fit, x)
    ax.plot(x, y, label=f'$y = {p0:.1f} + {p1:.1f} x$')

    ax.plot(b, zout, 'o', label='out-of-plane')

    p1, p0 = fit = np.polyfit(b, zout, 1)
    y = np.polyval(fit, x)
    ax.plot(x, y, label=f'$y = {p0:.1f} + {p1:.1f} x$')

    for symbol, (i, o, b) in extra.items():
        ax.text(b, i, symbol)
        # ax.text(b, o, symbol)

    ax.set_xlabel('Bader charge')
    ax.set_ylabel('Born charge')
    ax.legend()
    ax.set_ylim(-7.5, 7.5)
    fig.savefig('bader-born.png')
    plt.show()
    return

    """
    plt.plot(i, o, '+')
    plt.plot(i, f3(i))
    plt.xlabel('Born (in-plane)')
    plt.ylabel('Born (out-of-plane)')
    plt.savefig('born-born-all.png')
    plt.show()
    """


if __name__ == '__main__':
    extract_data()
    plot()

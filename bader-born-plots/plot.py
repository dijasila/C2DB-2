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
        data[row.uid] = {'gap': row.gap,
                         'B_a': B_a}
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
    gaps = []
    for u, d in data.items():
        if 'Z_avv' not in d:
            continue
        nz += 1
        B_a = d['B_a']
        Z_avv = d['Z_avv']
        x = (d['gap'],
             abs(B_a).mean(),
             abs(np.trace(Z_avv[:, :2, :2], 0, 1, 2) / 2).mean())
        gaps.append(x)
        if abs(x[0] - 3.5) < 0.2 and x[1] < 0.1:
            print(u, d)
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
    bad = abs(zin) >= 7.5
    print(nz, len(b), sum(bad))

    return zin, zout, b, extra, gaps


def plot1():
    zin, zout, b, extra, gaps = analyze()
    ok = abs(zin) < 7.5
    zin = zin[ok]
    zout = zout[ok]
    b = b[ok]

    fig = plt.figure(figsize=(width, width * 0.8),
                     constrained_layout=True)
    ax = fig.add_subplot(111)

    ax.plot(b, zin, 'o', alpha=0.5, ms=2, label='in-plane')

    ax.plot(b, zout, 'o', alpha=0.5, ms=2, label='out-of-plane')

    x = [-5.0, 5.0]
    p1, p0 = fit = np.polyfit(b, zin, 1)
    y = np.polyval(fit, x)
    print(p0, p1)
    ax.plot(x, y,
            # label=f'$y = {p0:.1f} + {p1:.1f} x$'
            )

    p1, p0 = fit = np.polyfit(b, zout, 1)
    y = np.polyval(fit, x)
    print(p0, p1)
    ax.plot(x, y,
            # label=f'$y = {p0:.1f} + {p1:.1f} x$'
            )

    if 0:
        for symbol, (i, o, b) in extra.items():
            ax.text(b, i, symbol)
            # ax.text(b, o, symbol)

    ax.set_xlabel('Bader charge [e]')
    ax.set_ylabel('Born charge [e]')
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


def plot2():
    *_, gaps = analyze()
    g, b, z = np.array(gaps).T

    fig = plt.figure(figsize=(width, 0.8 * width),
                     constrained_layout=True)
    ax = fig.add_subplot(111)

    ax.plot(g, z, 'o', alpha=0.5, ms=2, label='Born charge')
    ax.plot(g, b, 'o', color='C3', alpha=0.5, ms=2, label='Bader charge')
    ax.set_xlabel('Band gap [eV]')
    ax.set_ylabel('Charge (mean absolute value) [e]')
    ax.legend()
    ax.set_xlim(0, 6)
    ax.set_ylim(0, 7)
    fig.savefig('gap-charge.png')
    plt.show()
    return


if __name__ == '__main__':
    extract_data()
    plot1()
    plot2()

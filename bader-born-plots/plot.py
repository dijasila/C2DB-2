from pathlib import Path
import pickle
import runpy

import numpy as np
from ase.db import connect
from ase.data import chemical_symbols
import matplotlib.pyplot as plt
from mendeleev import element


width = runpy.run_path('../rcparams.py')['columnwidth']

chis = {symbol: element(symbol).electronegativity('pauling')
        for symbol in chemical_symbols[1:90]}


def extract_data():
    if Path('data.pckl').is_file():
        return
    db = connect('../c2db-data.db')
    data = {}
    for row in db.select('first_class_material'):
        born = row.data.get('results-asr.borncharges.json')
        bader = row.data.get('results-asr.bader.json')
        B_a = bader['kwargs']['data']['bader_charges']
        symbols = bader['kwargs']['data']['sym_a']
        data[row.uid] = {'gap': row.gap,
                         'B_a': B_a,
                         'symbols': symbols}
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
        symbols = d['symbols']
        x = (d['gap'],
             abs(B_a).mean(),
             abs(np.trace(Z_avv[:, :2, :2], 0, 1, 2) / 2).mean())
        gaps.append(x)
        if abs(x[0] - 3.5) < 0.2 and x[1] < 0.1:
            print(u, d)
        chi_a = np.array([chis[symbol] for symbol in symbols])
        I_a = abs(chi_a - chi_a.mean())
        D.extend([(np.trace(Z_vv[:2, :2]) / 2,
                   Z_vv[2, 2],
                   B,
                   I)
                  for Z_vv, B, I in zip(Z_avv, B_a, I_a)])
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

    zin, zout, b, chi = np.array(D).T
    bad = abs(zin) >= 7.5
    print(nz, len(b), sum(bad))

    return zin, zout, b, chi, extra, gaps


def plot1():
    zin, zout, b, chi, extra, gaps = analyze()
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


def plot3():
    zin, zout, b, chi, extra, gaps = analyze()
    if 0:
        ok = abs(zin) < 7.5
        zin = zin[ok]
        zout = zout[ok]
        b = b[ok]

    i = np.argsort(chi)
    chi = chi[i]
    zin = zin[i]
    b = b[i]

    fig = plt.figure(figsize=(width, width * 0.8),
                     constrained_layout=True)
    ax = fig.add_subplot(111)

    ax.plot(chi, zin - b, 'o', alpha=0.5, ms=2)

    dc = 0.1
    X = []
    Y = []
    dY = []
    for i in range(10):
        c = dc + i * 2 * dc
        m = (c - dc < chi) & (chi <= c + dc)
        y = (zin - b)[m]
        if len(y) == 0:
            continue
        y0 = y.mean()
        dy = ((y - y0)**2).mean()**0.5
        X.append(c)
        Y.append(y0)
        dY.append(dy)

    ax.errorbar(X, Y, dY, dc, 'o')

    ax.set_xlabel('I(a)')
    ax.set_ylabel('Born(in-plane, a) -Bader(a) [e]')
    ax.set_ylim(-7, 7)
    # ax.legend()
    fig.savefig('ionicity1.png')
    plt.show()
    return


def plot4():
    zin, zout, b, chi, extra, gaps = analyze()
    if 0:
        ok = abs(zin) < 7.5
        zin = zin[ok]
        zout = zout[ok]
        b = b[ok]

    z = (zin * 2 + zout) / 3

    i = np.argsort(chi)
    chi = chi[i]
    z = z[i]
    b = b[i]

    fig = plt.figure(figsize=(width, width * 0.8),
                     constrained_layout=True)
    ax = fig.add_subplot(111)

    # ax.plot(chi, zin - b, 'o', alpha=0.5, ms=2)
    # ax.plot(chi, zin, 'o', alpha=0.5, ms=2, label='Zin')
    # ax.plot(chi, b, 'o', alpha=0.5, ms=2, label='B')
    s = ax.scatter(b, z, c=chi, s=2, alpha=0.5)

    x = [-5, 5]
    """
    m = chi > 1
    p1, p0 = fit = np.polyfit(b[m], z[m], 1)
    y = np.polyval(fit, x)
    print(p0, p1)
    ax.plot(x, y,
            # label=f'$y = {p0:.1f} + {p1:.1f} x$'
            )
    x = [-5, 5]
    m = chi <= 1
    p1, p0 = fit = np.polyfit(b[m], z[m], 1)
    y = np.polyval(fit, x)
    print(p0, p1)
    """
    y = x
    ax.plot(x, y, 'C3')

    ax.set_xlabel('Bader charge [e]')
    ax.set_ylabel('Born charge (Tr$(Z)/3$) [e]')
    cbar = fig.colorbar(s)
    cbar.set_label('Ionicity')
    ax.set_xlim(-4, 4)
    ax.set_ylim(-4, 4)
    # ax.legend()
    fig.savefig('bader-born-ionicity.png')
    plt.show()
    return


if __name__ == '__main__':
    extract_data()
    # plot1()
    # plot2()
    plot4()

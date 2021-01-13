import pathlib

text = r"""
Magnetic anisotropies  & PBE(*) & Magnetic \\
Exchange couplings & PBE & Magnetic \\
Phonons ($\Gamma$ and BZ corners) & PBE & \\
Stiffness tensor & PBE & \\
Heat of formation & PBE & \\
Energy above convex hull & PBE & \\
Out-of-plane dipole & PBE & \\
Work function & PBE(*) & $E_{\mathrm{gap}}=0$ \\
Fermi surface & PBE(*) & $E_{\mathrm{gap}}=0$ \\
Deformation potentials & PBE(*) & $E_{\mathrm{gap}}>0$ \\
Bader charges & PBE & \\
Electronic band structure PBE & PBE(*) &   \\
Electronic band structure HSE & HSE06@PBE(*)  &   \\
Electronic band structure GW & G$_0$W$_0$@PBE(*) & $E_{\mathrm{gap}}>0$, $N_{\mathrm{atoms}}<5$ \\
Orbital projected band structure & PBE & \\
Projected density of states & PBE & \\
Effective masses & PBE(*) &  $E_{\mathrm{gap}}>0$ \\
Plasma frequency & PBE(*) & $E_{\mathrm{gap}}=0$ \\
Topological invariants & PBE(*) + Berry's phase & $0<E_{\mathrm{gap}}<0.3$ eV \\
Born charges & PBE + Berry's phase & $E_{\mathrm{gap}}>0$ \\
Spontaneous polarisation & PBE + Berry's phase &  $E_{\mathrm{gap}}>10$meV, nearly centrosym.(polar) \\
Piezoelectric tensor & PBE + Berry's phase & $E_{\mathrm{gap}}>10$meV, non-centrosym. \\ 
Optical polarisability  & RPA@PBE &   \\
Infrared polarisability & PBE & $E_{\mathrm{gap}}>0$  \\
Optical absorbance  &  BSE@PBE-scissors(*) &  $E_{\mathrm{gap}}>0$,  $N_{\mathrm{atoms}}<5$  \\
Raman spectrum & PBE + DZP basis set & NM, dyn. stable, $k_d=25$ {\AA}$^{-1}$ \\
Second harmonic generation & PBE & $E_{\mathrm{gap}}>10 meV$, NM, non-centrosym., $k_d=30$ {\AA}$^{-1}$"""


column_to_name = {
    'Magnetic anisotropies': 'magnetic_anisotropy',
    'Exchange couplings': 'exchange',
    r'Phonons ($\Gamma$ and BZ corners)': 'phonons',
    'Stiffness tensor': 'stiffness',
    'Heat of formation': 'convex_hull',
    'Energy above convex hull': 'convex_hull',
    'Out-of-plane dipole': 'gs',
    'Work function': 'gs',
    'Fermi surface': 'fermisurface',
    'Deformation potentials': 'deformationpotentials',
    'Bader charges': 'bader',
    'Electronic band structure PBE': 'bandstructure',
    'Electronic band structure HSE': 'hse',
    'Electronic band structure GW': 'gw',
    'Orbital projected band structure': 'projected_bandstructure',
    'Projected density of states': 'pdos',
    'Effective masses': 'emasses',
    'Plasma frequency': 'plasmafrequency',
    'Topological invariants': 'berry',
    'Born charges': 'borncharges',
    'Spontaneous polarisation': 'ferroelectricity',
    'Piezoelectric tensor': 'piezoelectrictensor',
    'Optical polarisability': 'polarizability',
    'Infrared polarisability': 'infraredpolarizability',
    'Optical absorbance': 'bse',
    'Raman spectrum': 'raman',
    'Second harmonic generation': 'shg',
}

count_kvp = pathlib.Path('count_kvp.txt').read_text()

name_to_count = {}
for line in count_kvp.split('\n'):
    line = line.replace(' ', '')
    if not line.startswith('has_asr'):
        continue
    name, count = line.split(':')
    count = int(count)
    name = name[8:]
    name_to_count[name] = count

name_to_count.update(
    {
        'ferroelectricity': 'XXX No count yet',
        'shg': 'XXX No count yet',
    }
)

lines = text.split(r'\\')

rows = []
for line in lines:
    columns = line.split('&')
    columns = [column.strip() for column in columns]
    # Append number of materials calculated
    name = column_to_name[columns[0]]
    count = name_to_count[name]
    columns.append(str(count))
    rows.append(columns)
    # print("'" + columns[0] + "': '',")

rows = sorted(rows, key=lambda x: (len(x[1]), len(x[2]), x[0]))


for row in rows:
    print(' ' * 8 + ' & '.join(row) + r' \\')

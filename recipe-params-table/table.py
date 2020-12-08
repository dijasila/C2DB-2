text = r"""
Magnetic anisotropies  & PBE(*) & Magnetic \\
Exchange couplings & PBE & Magnetic \\
Phonons ($\Gamma$ and BZ corners) & PBE & \\
Stiffness tensor & PBE & \\
Heat of formation & PBE & \\
Energy above convex hull & PBE & \\
Out-of-plane dipole & PBE & \\
Work function & PBE(*) & $E_{\mathrm{gap}}=0$ \\
\textcolor{red}{Fermi arc} & PBE(*) & $E_{\mathrm{gap}}=0$ \\
Deformation potentials & PBE(*) & $E_{\mathrm{gap}}>0$ \\
Bader charges & PBE & \\
Electronic band structure  & PBE(*) &   \\
Electronic band structure  & HSE06@PBE(*)  &   \\
Electronic band structure  & G$_0$W$_0$@PBE(*) & $E_{\mathrm{gap}}>0$, $N_{\mathrm{atoms}}<5$ \\
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
lines = text.split(r'\\')

rows = []
for line in lines:
    columns = line.split('&')
    columns = [column.strip() for column in columns]
    rows.append(columns)

rows = sorted(rows, key=lambda x: (len(x[1]), len(x[2]), x[0]))

for row in rows:
    print(' ' * 8 + ' & '.join(row) + r' \\')

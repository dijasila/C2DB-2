import matplotlib.pyplot as plt
import numpy as np

def f(x):
    return x**2

x1 = np.linspace(-1, 1, 100)
x2 = np.linspace(-0.8, 1.2, 100)

# ppi = 100
# figw = 600
# figh = 600
# 
# fig = plt.figure(figsize=(figw/ppi, figh/ppi), dpi=ppi)
font_ticks = 16
font_labels = 18
delta = 0.05

arrow_params = {'shape': 'full', 'color': 'black', 'width': 0.015, 'length_includes_head': True, 'head_length': 0.1}
# plt.arrow(0.6, 2, 0, -2, shape='full', color='black', width=0.02, length_includes_head=True, head_length=0.2)
plt.arrow(0.6, 2 - delta, 0, -2 + delta, **arrow_params)
plt.arrow(0.6, 0 + delta, 0, 2 - delta, **arrow_params)
plt.arrow(0, 0, 0, 2.35, **arrow_params)
plt.arrow(1.2, 2 + delta, 0, 0.35 - delta, **arrow_params)
plt.arrow(1.2, 2.35 - delta, 0, -0.35 + delta, **arrow_params)
plt.arrow(-0.6, 0 + delta, 0, 0.35 - delta, **arrow_params)
plt.arrow(-0.6, 0.35 - delta, 0, -0.35 + delta, **arrow_params)
plt.text(-0.8, 2.5, '$D_{excited}$', fontsize=font_ticks, color='C1')
plt.text(1.1, 0.5, '$D_{ground}$', fontsize=font_ticks, color='C0')
plt.text(1.2 + delta, 2.1725, '$E_{excited}^{reorg}$', fontsize=font_ticks, verticalalignment='center', horizontalalignment='left')
plt.text(0.6 + delta, 1.0, '$E_{ZPL}$', fontsize=font_ticks, verticalalignment='center', horizontalalignment='left')
plt.text(0.0 - delta, 1.1725, 'Excitation', fontsize=font_ticks, verticalalignment='center', horizontalalignment='right', rotation=90)
plt.text(-0.6 - delta, 0.1725, '$E_{ground}^{reorg}$', fontsize=font_ticks, verticalalignment='center', horizontalalignment='right')

plt.plot(x1, f(x1))
plt.plot(x1 + 0.6, f(x1)+2)
plt.xticks([], [])
plt.yticks([0, 2], [0, "$E_{excited}$"], fontsize=font_ticks)
# plt.yticks([],[])
plt.xlabel('$Q$', fontsize=font_labels)
# plt.ylabel("$E - E_{ground}$", fontsize=font_labels)
plt.ylabel("Energy", fontsize=font_labels)
plt.axhline(0, linestyle='dotted', color='black')
plt.axhline(2, linestyle='dotted', color='black')
plt.tight_layout()
plt.savefig('cc_exc.png', dpi=150)

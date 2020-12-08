from matplotlib import pyplot as plt
import numpy as np

labels = ['h-BN', r'1H-MoS$_2$', r'1H-MoSe$_2$', r'1H-MoTe$_2$',
          r'1H-WS$_2$', r'1H-WSe$_2$', r'1H-WTe$_2$']
ref = [3.71, 3.06, 2.8, 2.98, 2.20, 1.93, 1.60]
mydata = [3.68, 2.92, 2.69, 2.53, 2.04, 1.82, 1.54]

ref = [1.38, 3.64, 3.92, 5.43, 2.47, 2.71, 3.4]
reflabels = [''] * 7
mydata = [1.33, 3.52, 3.83, 4.8, 2.37, 2.62, 3.41]
# plt.figure(figsize=(4, 3))
# x = np.arange(0, len(mydata))
# plt.plot(x, mydata, 'o', color='C0', label='C2DB (GPAW)')
# plt.plot(x, ref, 'o', color='C1', label='Duerloo et al. [*]',)
# plt.plot([1], [3.0], 'o', color='C2', label='MoS$_2$ Exp. [**]')
# plt.xticks(x, labels, rotation=30)
# plt.xlim(-1, 9)
# plt.ylabel(r'$e_{xx}$ ($10^{-10}$ C / m)')
# plt.legend()
# plt.tight_layout()
# plt.savefig('comparetoref.pdf')


for reference, our_number, material, reflabel in zip(
        ref, mydata, labels, reflabels):
    print(' & '.join([
        material,
        '',
        format(reference / 10, '.2f') + reflabel,
        format(our_number / 10, '.2f')]) + r' \\')

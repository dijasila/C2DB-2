import sys


def f(x, a, b):
    return a * x + b


def order_transitions(trans_dict):
    translist = []
    reflist = ['0/1', '1/2', '2/3', '-1/0', '-2/-1', '-3/-2']
    for element in trans_dict:
        if element in reflist:
            translist.append(element)

    ordered_dict = {}
    for element in reflist:
        if element in trans_dict:
            ordered_dict[element] = trans_dict[element]

    print(f'Ordered dictionary: {ordered_dict}')

    return ordered_dict


def get_b(x, y, a):
    return y - a * x


def plot_formation_energies(defname, fname):
    import matplotlib.pyplot as plt
    from asr.core import read_json

    fig, ax = plt.subplots(1, 2)
    defnamelist = ['vW', 'SeW']
    for l in [0, 1]:
        data = read_json(f'results-asr.sj_analyze_{defnamelist[l]}.json')

        vbm = data['pristine']['vbm'] - data['pristine']['evac']
        cbm = data['pristine']['cbm'] - data['pristine']['evac']
        gap = cbm - vbm
        eform = data['eform']

        transitions = data['transitions']

        fontsize_ticks = 16
        fontsize_labels = 18

        ax[l].fill_betweenx([-10, 30], vbm - 10, vbm, color='C0', alpha=0.5)
        ax[l].fill_betweenx([-10, 30], cbm + 10, cbm, color='C1', alpha=0.5)
        ax[l].axhline(0, color='black', linestyle='dotted')

        ax[l].set_xlim(vbm - 0.1 * gap, cbm + 0.1 * gap)
        ax[l].set_ylim(-0.1, eform + 0.2 * eform)
        e_m = transitions["-1/0"][0] - transitions["-1/0"][1] - transitions["-1/0"][2]
        e_p = transitions["0/1"][0] - transitions["0/1"][1] - transitions["0/1"][2]
        ax[l].plot([vbm, cbm], [eform, eform], color='black')
        if e_m < cbm and e_m > vbm:
            ax[l].axvline(e_m, color='black', linestyle='-.')
        if e_p < cbm and e_p > vbm:
            ax[l].axvline(e_p, color='black', linestyle='-.')

        transitions = order_transitions(transitions)

        enlist = []
        for element in transitions:
            enlist.append(transitions[element][0] - transitions[element][1] - transitions[element][2])

        ax2 = ax[l].twiny()

        tickslist = []
        labellist = []
        for i, element in enumerate(transitions):
            energy = transitions[element][0] - transitions[element][1] - transitions[element][2]
            enlist.append(energy)
            name = element
            if energy > vbm and energy < cbm:
                ax[l].axvline(energy, color='grey', linestyle='-.')
                if name.split('/')[1].startswith('0') and name.split('/')[0].startswith('-'):
                    y1 = eform
                    y2 = eform
                elif name.split('/')[0].startswith('0'):
                    y1 = eform
                    y2 = eform
                elif not name.split('/')[0].startswith('-'):
                    y2 = None
                else:
                    y1 = None
                if name.split('/')[0].startswith('-'):
                    tickslist.append(energy)
                    labellist.append(name)
                    a = float(name.split('/')[0])
                    b = get_b(enlist[i], y2, a)
                    if y1 is None:
                        y1 = f(enlist[i], a, b)
                    y2 = f(enlist[i + 1], a, b)
                    print(enlist[i], enlist[i+1], y1, y2, a)
                    ax[l].plot([vbm, cbm], [f(vbm, a, b), f(cbm, a, b)], color='black')
                    ax[l].plot([enlist[i], enlist[i + 1]], [y1, y2], color='black', marker='s')
                elif not name.split('/')[0].startswith('-'):
                    tickslist.append(energy)
                    labellist.append(name)
                    a = float(name.split('/')[0])
                    b = get_b(enlist[i], y1, a)
                    if y2 is None:
                        y2 = f(enlist[i], a, b)
                    y1 = f(enlist[i + 1], a, b)
                    print(enlist[i], enlist[i+1], y1, y2, a)
                    ax[l].plot([vbm, cbm], [f(vbm, a, b), f(cbm, a, b)], color='black')
                    ax[l].plot([enlist[i], enlist[i + 1]], [y1, y2], color='black', marker='s')
        ax[l].set_xlabel('$E_F$ (eV)', fontsize=fontsize_labels)
        ax[l].tick_params(axis='both', labelsize=fontsize_ticks)
        ax2.set_xlim(ax[l].get_xlim())
        ax2.set_xticks(tickslist)
        ax2.set_xticklabels(labellist, fontsize=fontsize_ticks)

    ax[0].set_ylabel('Formation energy (eV)', fontsize=fontsize_labels)
    plt.tight_layout()
    ax[0].text(-4, 0.5, '$v_W$', fontsize=fontsize_labels)
    ax[1].text(-4, 0.5, '$Se_W$', fontsize=fontsize_labels)
    plt.savefig(f'{fname}.png', dpi=150)
    plt.close()



fname = f'formation'

plot_formation_energies(None, fname)

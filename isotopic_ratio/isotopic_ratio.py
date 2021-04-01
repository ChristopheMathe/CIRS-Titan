# Author: Christophe Mathé
# Date: 19/03/2020
# ----------------------------------------------------------------------
# Standard: HCN       => step2_hcn_l1
#           hcn3/hcn5 => step3_hcn3_hcn5_l200
# --------------------------------------------
# 27/10/2020:
# -----------
# New isotopic ratios calculations as asked by Sandrine


from os import listdir
from os.path import isdir, isfile
from numpy import loadtxt, mean, sqrt, zeros, asarray, average
import matplotlib.pyplot as plt
from math import sqrt


class DirectoryIsotope(object):
    def __init__(self, name=None, dir_hcn=None, dir_h3c5n=None):
        self.name = name
        self.dir_hcn = dir_hcn
        self.dir_h3c5n = dir_h3c5n


def dictionary_directory():
    tb = DirectoryIsotope(name='Tb', dir_hcn='step2_hcn_l1_lb1-3', dir_h3c5n='step3_hcn3_hcn5_l200_lb1-3')
    t003 = DirectoryIsotope(name='T003', dir_hcn='step2_hcn_l1_lb1-3', dir_h3c5n='step3_hcn3_hcn5_l200_lb1-3')
    t004 = DirectoryIsotope(name='T004', dir_hcn='step2_hcn_l1_AltRes_3km_1e-4mbar',
                            dir_h3c5n='step3_hcn3_hcn5_l200_AltRes_3km_1e-4mbar')
    t006 = DirectoryIsotope(name='T006', dir_hcn='step2_hcn_l1', dir_h3c5n='step3_hcn3_hcn5_l200')
    t014 = DirectoryIsotope(name='T014', dir_hcn='step2_hcn_l1_lb2-6', dir_h3c5n='step3_hcn3_hcn5_l200_lb2-6')
    t015 = DirectoryIsotope(name='T015', dir_hcn='step2_hcn_l1_lb2-6', dir_h3c5n='step3_hcn3_hcn5_l200_lb2-6')
    t016 = DirectoryIsotope(name='T016', dir_hcn='step2_hcn_l1_lb2-4', dir_h3c5n='step3_hcn3_hcn5_l200_lb2-4')
    t018 = DirectoryIsotope(name='T018', dir_hcn='step2_hcn_l1_AltRes3km_lb1-14',
                            dir_h3c5n='step3_hcn3_hcn5_l200_AltRes3km_lb1-14')
    t019 = DirectoryIsotope(name='T019', dir_hcn='step2_hcn_l1_lb1-5_704-718cm-1',
                            dir_h3c5n='step3_hcn3_hcn5_l200_lb1-5_704-718cm-1')
    t021 = DirectoryIsotope(name='T021', dir_hcn='step2_hcn_l1_lb1-4', dir_h3c5n='step3_hcn3_hcn5_l200_lb1-4')
    t023 = DirectoryIsotope(name='T023', dir_hcn='step2_hcn_l1_lb1-3', dir_h3c5n='step3_hcn3_hcn5_l200_lb1-3')
    t024 = DirectoryIsotope(name='T024', dir_hcn='step2_hcn_l1_lb1-5', dir_h3c5n='step3_hcn3_hcn5_l200_lb1-5')
    t027 = DirectoryIsotope(name='T027', dir_hcn='step2_hcn_l1_lb1-4', dir_h3c5n='step3_hcn3_hcn5_l200_lb1-4')
    t028 = DirectoryIsotope(name='T028', dir_hcn='step2_hcn_l1_lb2-5', dir_h3c5n='step3_hcn3_hcn5_l200_lb2-5')
    t035 = DirectoryIsotope(name='T035', dir_hcn='step2_hcn_l1_lb1-4', dir_h3c5n='step3_hcn3_hcn5_l200_lb1-4')
    t039 = DirectoryIsotope(name='T039', dir_hcn='step2_hcn_l1_lb1-4', dir_h3c5n='step3_hcn3_hcn5_l200_lb1-4')
    t042 = DirectoryIsotope(name='T042', dir_hcn='step2_hcn_l1_lb2-6', dir_h3c5n='step3_hcn3_hcn5_l200_lb2-6')
    t043 = DirectoryIsotope(name='T043', dir_hcn='step2_hcn_l1_lb1-4', dir_h3c5n='step3_hcn3_hcn5_l200_lb1-4')
    t047 = DirectoryIsotope(name='T047', dir_hcn='step2_hcn_l1_lb2-5', dir_h3c5n='step3_hcn3_hcn5_l200_lb2-5')
    t049 = DirectoryIsotope(name='T049', dir_hcn='step2_hcn_l1_lb2-4', dir_h3c5n='step3_hcn3_hcn5_l200_lb2-4')
    t051 = DirectoryIsotope(name='T051', dir_hcn='step2_hcn_l1_lb1-5', dir_h3c5n='step3_hcn3_hcn5_l200_lb1-5')
    t054 = DirectoryIsotope(name='T054', dir_hcn='step2_hcn_l1', dir_h3c5n='step3_hcn3_hcn5_l200')
    t057 = DirectoryIsotope(name='T057', dir_hcn='step2_hcn_l1_lb2-8', dir_h3c5n='step3_hcn3_hcn5_l200_lb2-8')
    t058 = DirectoryIsotope(name='T058', dir_hcn='step2_hcn_l1_lb2-8', dir_h3c5n='step3_hcn3_hcn5_l200_lb2-8')
    t059 = DirectoryIsotope(name='T059', dir_hcn='step2_hcn_l1_lb1-3_705-718cm-1',
                            dir_h3c5n='step3_hcn3_hcn5_l200_lb1-3_705-718cm-1')
    t061 = DirectoryIsotope(name='T061', dir_hcn='step2_hcn_l1_lb1-4_alpha3',
                            dir_h3c5n='step3_hcn3_hcn5_l200_lb1-4_alpha20')
    t062 = DirectoryIsotope(name='T062', dir_hcn='step2_hcn_l1_lb3-7', dir_h3c5n='step3_hcn3_hcn5_l200_lb3-7')
    t064 = DirectoryIsotope(name='T064', dir_hcn='step2_hcn_l1_705-718cm-1_lb1-3',
                            dir_h3c5n='step3_hcn3_hcn5_l200_705-718cm-1_lb1-3')
    t065 = DirectoryIsotope(name='T065', dir_hcn='step2_hcn_l1_lb1-2', dir_h3c5n='step3_hcn3_hcn5_l200_lb1-2')
    t067 = DirectoryIsotope(name='T067', dir_hcn='step2_hcn_l1_706-718cm-1_lb2-7',
                            dir_h3c5n='step3_hcn3_hcn5_l200_706-718cm-1_lb2-7')
    t076 = DirectoryIsotope(name='T076', dir_hcn='step2_hcn_l1_lb6-13', dir_h3c5n='step3_hcn3_hcn5_l200_lb6-13')
    t078 = DirectoryIsotope(name='T078', dir_hcn='step2_hcn_l1_lb1-3', dir_h3c5n='step3_hcn3_hcn5_l200_lb1-3')
    t079 = DirectoryIsotope(name='T079', dir_hcn='step2_hcn_l1_lb2-4_706-718cm-1',
                            dir_h3c5n='step3_hcn3_hcn5_l200_lb2-4_706-718cm-1')
    t082 = DirectoryIsotope(name='T082', dir_hcn='step2_hcn_l1_lb1-3', dir_h3c5n='step3_hcn3_hcn5_l200_lb1-3')
    t083 = DirectoryIsotope(name='T083', dir_hcn='step2_hcn_l1_lb1-5', dir_h3c5n='step3_hcn3_hcn5_l200_lb1-5')
    t084 = DirectoryIsotope(name='T084', dir_hcn='step2_hcn_l1_lb1-5', dir_h3c5n='step3_hcn3_hcn5_l200_lb1-5')
    t090 = DirectoryIsotope(name='T090', dir_hcn='step2_hcn_l1_lb1-4', dir_h3c5n='step3_hcn3_hcn5_l200_lb1-4')
    t092 = DirectoryIsotope(name='T092', dir_hcn='step2_hcn_l1_lb1-4', dir_h3c5n='step3_hcn3_hcn5_l200_lb1-4')
    t095 = DirectoryIsotope(name='T095', dir_hcn='step2_hcn_l1_lb1-3', dir_h3c5n='step3_hcn3_hcn5_l200_lb1-3')
    t097 = DirectoryIsotope(name='T097', dir_hcn='step2_hcn_l1_lb2-5', dir_h3c5n='step3_hcn3_hcn5_l200_lb2-5')
    t101 = DirectoryIsotope(name='T101', dir_hcn='step2_hcn_l1_lb2-5', dir_h3c5n='step3_hcn3_hcn5_l200_lb2-5')
    t102 = DirectoryIsotope(name='T102', dir_hcn='step2_hcn_l1_lb2-5', dir_h3c5n='step3_hcn3_hcn5_l200_lb2-5')
    t103 = DirectoryIsotope(name='T103', dir_hcn='step2_hcn_l1_lb1-4', dir_h3c5n='step3_hcn3_hcn5_l200_lb1-4')
    t106 = DirectoryIsotope(name='T106', dir_hcn='step2_hcn_l1_lb1-3', dir_h3c5n='step3_hcn3_hcn5_l200_lb1-3')
    t108 = DirectoryIsotope(name='T108', dir_hcn='step2_hcn_l1_lb2-6_706-718cm-1',
                            dir_h3c5n='step3_hcn3_hcn5_l200_lb2-6_706-718cm-1')
    t110 = DirectoryIsotope(name='T110', dir_hcn='step2_hcn_l1', dir_h3c5n='step3_hcn3_hcn5_l200_bis')
    t113 = DirectoryIsotope(name='T113', dir_hcn='step2_hcn_l1_1lbp', dir_h3c5n='step3_hcn3_hcn5_l200_1lbp')
    t115 = DirectoryIsotope(name='T115', dir_hcn='', dir_h3c5n='')
    t117 = DirectoryIsotope(name='T117', dir_hcn='step2_hcn_l1_lb3_C2H2nonmod_HCNnonmod',
                            dir_h3c5n='step3_hcn3_hcn5_l200_lb3_C2H2nonmod_HCNnonmod')
    t120 = DirectoryIsotope(name='T120', dir_hcn='step2_hcn_l1_lb3-5', dir_h3c5n='step3_hcn3_hcn5_l200_lb3-5')
    t121 = DirectoryIsotope(name='T121', dir_hcn='step2_hcn_l1_lb1-4', dir_h3c5n='step3_hcn3_hcn5_l200_lb1-4')
    t123 = DirectoryIsotope(name='T123', dir_hcn='step2_hcn_l1', dir_h3c5n='step3_hcn3_hcn5_l200')
    t125 = DirectoryIsotope(name='T125', dir_hcn='step2_hcn_l1_lb3-5', dir_h3c5n='step3_hcn3_hcn5_l200_lb3-5')
    t126 = DirectoryIsotope(name='T126', dir_hcn='step2_hcn_l1_lb2-5', dir_h3c5n='step3_hcn3_hcn5_l200_lb2-5')
    s98 = DirectoryIsotope(name='S98', dir_hcn='step2_hcn_l1_lb1-3', dir_h3c5n='step3_hcn3_hcn5_l200_lb1-3')

    dict_dir = {'Tb': tb, 't003': t003, 'T004': t004, 'T006': t006, 'T014': t014, 'T015': t015, 'T016': t016,
                'T018': t018, 'T019': t019, 'T021': t021, 'T023': t023, 'T024': t024, 'T027': t027, 'T028': t028,
                'T035': t035, 'T039': t039, 'T042': t042, 'T043': t043, 'T047': t047, 'T049': t049, 'T051': t051,
                'T054': t054, 'T057': t057, 'T058': t058, 'T059': t059, 'T061': t061, 'T062': t062, 'T064': t064,
                'T065': t065, 'T067': t067, 'T076': t076, 'T078': t078, 'T079': t079, 'T082': t082, 'T083': t083,
                'T084': t084, 'T090': t090, 'T092': t092, 'T095': t095, 'T097': t097, 'T101': t101, 'T102': t102,
                'T103': t103, 'T106': t106, 'T108': t108, 'T110': t110, 'T113': t113, 'T115': t115, 'T117': t117,
                'T120': t120, 'T121': t121, 'T123': t123, 'T125': t125, 'T126': t126, 'S98': s98}
    return dict_dir


def check_directory(file):
    if isfile(file) is False:
        print(f'Directory does not exist in path: {file}')
        check = False
    else:
        check = True
    return check


def get_latitude(file):
    with open(file, 'r') as file_read:
        lines = file_read.readlines()
        index = [i for i, s in enumerate(lines) if 'FP3LAT' in s][0]
        latitude = lines[index].split()[2]
    return latitude


def get_factor(file_log):
    file_read = open(file_log, 'r')
    data = file_read.readlines()
    sub_string = 'Spectral calculations for'
    nb_obs = int([s for s in data if sub_string.lower() in s.lower()][0].split()[3])
    sub_string_check = 'altitudes are hereafter considered'
    if any(sub_string_check in s for s in data):
        nb_obs = int([s for s in data if sub_string_check.lower() in s.lower()][0].split()[1])

    pattern = ' *** rms ***\n'
    indices = [k for k, x in enumerate(data) if x == pattern]
    indices = indices[-1]

    tmp_rms = zeros(nb_obs)
    for j in range(nb_obs):
        tmp = asarray(data[indices + 2 * (j + 1)].split())
        tmp_rms[j] = float(tmp[-1])
    rms = mean(tmp_rms)
    nesr = float(data[indices - (nb_obs + 1)].split()[1])
    if rms < nesr:
        factor = 1.0
    else:
        factor = (rms / nesr)
    return factor


def computation_ratio(hcn, isotope, factor):
    ratio = mean(isotope[:, 1] / hcn[:, 1])
    noise_hcn = mean((hcn[:, 3] - hcn[:, 1]) / hcn[:, 1])
    noise_isotope = mean((isotope[:, 3] - isotope[:, 1]) / isotope[:, 1])
    noise_temperature = 0.01
    noise_shift = 0.01
    err = sqrt(noise_shift ** 2 + noise_temperature ** 2 + (noise_isotope * factor) ** 2 + (noise_hcn * factor) ** 2)

    err = ratio * err

    return ratio, err


def draw_figure(latitude, ratio, ratio_err, save_name, hcn3):
    weight = 1. / (ratio_err ** 2)
    ratio_mean = average(ratio, axis=None, weights=weight)
    sigma_mean = sqrt(1.0 / sum(weight))

    xi2 = sum(((ratio - ratio_mean) / ratio_err) ** 2) / (ratio.shape[0] * 1.0 - 1.0)

    plt.figure(figsize=(8, 11))
    plt.errorbar(latitude, ratio, yerr=ratio_err, fmt='o', color='grey', capsize=8)

    plt.hlines(ratio_mean, -90, 90, linestyles='-', colors='black', label='Mean')
    plt.hlines(ratio_mean + sigma_mean, -90, 90, linestyles='--', colors='black')
    plt.hlines(ratio_mean - sigma_mean, -90, 90, linestyles='--', colors='black')

    plt.xlabel(u'Latitude (°N)')
    if hcn3:
        plt.title(r'$^{13}$C/$^{12}$C isotopic ratio: ' + f'{ratio_mean:3.2e} +/- {sigma_mean:3.2e},' +
                  r'$\chi^2_\u03BD$ = ' + f'{xi2:3.2f}')
        plt.hlines(1 / 89., -90, 90, color='red', label='Earth value')
        plt.ylabel('$^{13}$C/$^{12}$C ratio')
        print(f'12C/13C = {1.0 / ratio_mean:.0f} '
              f'({1.0 / (ratio_mean + sigma_mean):.0f} :'
              f' {1.0 / (ratio_mean - sigma_mean):.0f})')
    else:
        plt.title(r'$^{15}$N/$^{14}$N isotopic ratio: ' + "%3.2e" % ratio_mean + ' +/- ' + '%3.2e' % sigma_mean +
                  r', $\chi^2_\u03BD$ = ' + '%3.2f' % xi2)
        plt.hlines(1 / 58., -90, 90, color='red', label='')
        plt.ylabel('$^{15}$N/$^{14}$N ratio')
        print(f'14N/15N = {1.0 / ratio_mean:.0f}'
              f' ({1.0 / (ratio_mean + sigma_mean):.0f} :'
              f' {1.0 / (ratio_mean - sigma_mean):.0f})')

    plt.xlim(-90, 90)
    plt.legend(loc='best')
    plt.savefig(save_name + '.png', bbox_inches='tight')
    plt.savefig(save_name + '.pdf', bbox_inches='tight')


def main():
    # Get every directory used for retrievals of isotopic ratios
    list_directory = [d for d in listdir('.') if isdir(d)]
    list_directory.remove('T115')
    list_directory.remove('fit_mail_19_06_2020')
    nb_obs = len(list_directory)

    # Initialization
    ratio_c = zeros(nb_obs)
    ratio_n = zeros(nb_obs)
    ratio_c_err = zeros(nb_obs)
    ratio_n_err = zeros(nb_obs)
    latitude = zeros(nb_obs)

    # Compute isotopic ratios for each directory
    for i, value_i in enumerate(list_directory):
        filename_hcn = value_i + '/' + dictionary_directory()[value_i].dir_hcn + '/profq_HCN_retrieved.dat'
        filename_hcn3 = value_i + '/' + dictionary_directory()[value_i].dir_h3c5n + '/profq_hcn3_retrieved.dat'
        filename_hcn5 = value_i + '/' + dictionary_directory()[value_i].dir_h3c5n + '/profq_hcn5_retrieved.dat'
        filename_log = value_i + '/' + dictionary_directory()[value_i].dir_h3c5n + '/log.out'
        filename_info = value_i + '/information.txt'

        check = check_directory(filename_hcn5)
        if not check:
            print(f'\t FALSE: {value_i} {check}')
            continue
        else:
            print(value_i, check)
            print(f'\t HCN file path: {filename_hcn}')
            print(f'\t H13CN file path: {filename_hcn3}')
            print(f'\t HC15N file path: {filename_hcn5}')

        data_hcn = loadtxt(filename_hcn)
        data_hcn3 = loadtxt(filename_hcn3)
        data_hcn5 = loadtxt(filename_hcn5)

        factor = get_factor(filename_log)
        latitude[i] = get_latitude(filename_info)

        ratio_c[i], ratio_c_err[i] = computation_ratio(data_hcn, data_hcn3, factor)
        ratio_n[i], ratio_n_err[i] = computation_ratio(data_hcn, data_hcn5, factor)
        print(f'\t 12C/13C = {round(1. / ratio_c[i], 0):.0f}')
        print(f'\t 14N/15N = {round(1. / ratio_n[i], 0):.0f}')

    # Draw figures
    draw_figure(latitude, ratio_c, ratio_c_err, save_name='carbon_ratio_new', hcn3=True)
    draw_figure(latitude, ratio_n, ratio_n_err, save_name='nitrogen_ratio_new', hcn3=False)
    return


if '__main__' == __name__:
    main()

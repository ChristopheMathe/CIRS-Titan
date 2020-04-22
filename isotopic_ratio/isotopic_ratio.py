# Author: Christophe Mathé
# Date: 19/03/2020
#----------------------------------------------------------------------
# Standard: HCN       => step2_hcn_l1
#           hcn3/hcn5 => step3_hcn3_hcn5_l200
#

from os import listdir
from os.path import isdir, isfile
from numpy import loadtxt, mean, sqrt, zeros, max, asarray
import matplotlib.pyplot as plt


def check_directory(file):
    if isfile(file) is False:
        print('Directory does not exist in path: ', file)
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

    object = ' *** rms ***\n'
    indices = [k for k, x in enumerate(data) if x == object]
    indices = indices[-1]

    tmp_rms = zeros((nb_obs))
    for j in range(nb_obs):
        tmp = asarray(data[indices + 2 * (j + 1)].split())
        tmp_rms[j] = float(tmp[-1])
    rms = max(tmp_rms[:4])
    nesr = float(data[indices - (nb_obs + 1)].split()[1])
    if rms < nesr:
        factor = 1.0
    else:
        factor = (rms / nesr)
    return factor


def calcul_ratio(hcn, isotope, factor):
    ratio = mean(isotope[:, 1] / hcn[:, 1])
    noise_hcn = mean((hcn[:, 3] - hcn[:, 1]) / hcn[:, 1])
    noise_isotope = mean((isotope[:, 3] - isotope[:, 1]) / isotope[:, 1])
    err = sqrt(0.01 ** 2 + 0.01 ** 2 + (noise_isotope * factor) ** 2 + (noise_hcn * factor) ** 2)

    err = ratio * err

    return ratio, err


def draw_figure(latitude, ratio, ratio_err, hcn3):
    plt.figure(figsize=(8, 11))
    plt.errorbar(latitude, ratio, yerr=ratio_err, fmt='o', color='black')

    plt.xlabel(u'Latitude (°N)')
    if hcn3:
        plt.hlines(1/89., -90, 90, color='red')
        plt.ylabel('$^{13}$C/$^{12}$C ratio')
        savename = 'carbon_ratio'
    else:
        plt.hlines(1/58., -90, 90, color='red', label='')
        plt.ylabel('$^{15}$N/$^{14}$N ratio')
        savename = 'nitrogen_ratio'

    plt.xlim(-90, 90)
    plt.legend(loc='best')
    plt.savefig(savename+'.ps', bbox_inches='tight')
    plt.savefig(savename+'.eps', bbox_inches='tight')
    plt.savefig(savename+'.png', bbox_inches='tight')
    plt.savefig(savename+'.pdf', bbox_inches='tight')


def main():
    list_directory = [d for d in listdir('.') if isdir(d)]
    nb_obs = len(list_directory)

    ratio_c = zeros(nb_obs)
    ratio_n = zeros(nb_obs)
    ratio_c_err = zeros(nb_obs)
    ratio_n_err = zeros(nb_obs)
    latitude = zeros(nb_obs)

    for i, value_i in enumerate(list_directory):
        filename_hcn = value_i+'/step2_hcn_l1/profq_HCN_retrieved.dat'
        filename_hcn3 = value_i+'/step3_hcn3_hcn5_l200/profq_hcn3_retrieved.dat'
        filename_hcn5 = value_i+'/step3_hcn3_hcn5_l200/profq_hcn5_retrieved.dat'

        check = check_directory(filename_hcn5)
        print(check)
        if not check:
            continue

        data_hcn = loadtxt(filename_hcn)
        data_hcn3 = loadtxt(filename_hcn3)
        data_hcn5 = loadtxt(filename_hcn5)

        factor = get_factor(value_i+'/step3_hcn3_hcn5_l200/log.out')
        latitude[i] = get_latitude(value_i+'/information.txt')

        ratio_c[i], ratio_c_err[i] = calcul_ratio(data_hcn, data_hcn3, factor)
        ratio_n[i], ratio_n_err[i] = calcul_ratio(data_hcn, data_hcn5, factor)

    draw_figure(latitude, ratio_c, ratio_c_err, hcn3=True)
    draw_figure(latitude, ratio_n, ratio_n_err, hcn3=False)


if '__main__' == __name__:
    main()

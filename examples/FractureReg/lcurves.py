#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  5 12:19:43 2025

@author: steph
"""

import csv
import numpy as np
import matplotlib.pyplot as plt




def extract_lcurve_point(csv_file):
    ## read final row of misfit_norm, btrc_norm, mucoef_norm, etc from Trys lcurve.csv file
    with open(csv_file) as file:
        reader = csv.reader(file)
        for row in reader:
            pass
    print (row)
    return row

def extract_array_lcurve(path, begin=0, n=18):

    misfit_norm = np.empty(n)    #
    btrc_norm = np.empty(n)
    mucoef_norm = np.empty(n)
    trys_norm = np.empty(n)

    for k in range(begin,begin+n):
        r = extract_lcurve_point(f'{path}/{k}/l_curve_params.csv')
        misfit_norm[k] = r[0]
        btrc_norm[k]= r[2]
        mucoef_norm[k] = r[3]
        trys_norm[k] = r[5]
    
    return misfit_norm, btrc_norm, mucoef_norm, trys_norm


# %% btrc-only (starting from uniform beta)
alpha_btrc = 10**np.linspace(-9,8,18)
misfit_norm, btrc_norm, mucoef_norm, trys_norm = extract_array_lcurve('./i3/281075-ideal_shelf_inv_btrc_tikh', begin=0, n=18)

fig, ax = plt.subplots(1,1)
ax.loglog(misfit_norm,btrc_norm,'ro-')
for i, a in enumerate(alpha_btrc):
    if btrc_norm[i] < 2e10:
        ax.annotate(f'{a:.1e}', (misfit_norm[i],btrc_norm[i]))


# %% btrc-only (starting from uniform beta)
alpha_mucoef = 10**np.linspace(-1,16,18)
misfit_norm, btrc_norm, mucoef_norm, trys_norm = extract_array_lcurve('./i3/281112-ideal_shelf_inv_known_btrc_mucoef_tikh', begin=0, n=16)

fig, ax = plt.subplots(1,1)
ax.loglog(misfit_norm,mucoef_norm,'ro-')
for i, a in enumerate(alpha_mucoef):
       ax.annotate(f'{a:.1e}', (misfit_norm[i],mucoef_norm[i]))
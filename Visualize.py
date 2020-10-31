#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 16:34:58 2020
@author: alex

Статистика
"""
folder = r'/mnt/D/Yandex/CubeSat/Articles/2020 JATM/py/'
import pickle
import numpy as np
import matplotlib.pyplot as plt


files = ['f1U_6', 'f1U_5', 'f1U_y', 'f1U_z']

plt.clf()
D = []
for fName in files:
    with open(folder+fName+'.pickle', 'rb') as f:
        data = pickle.load(f)
    D.append({'E':data['E'].flatten(), 'w':data['w']})
    E = data['E']
    if (fName=='f1U_y') or (fName=='f1U_z'):
        E = 5*E

    hist, bins = np.histogram(E, bins=100, density=True)
    bin_centers = (bins[1:]+bins[:-1])*0.5
    plt.plot(bin_centers, hist)

    print('%s: mean %.2f, 0.9 in (%.2f-%.2f)'%(fName, np.mean(E), np.percentile(E, 5), np.percentile(E, 95)))

plt.legend(files)




plt.savefig(folder+'1U.svg')




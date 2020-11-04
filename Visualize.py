#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 16:34:58 2020
@author: alex

Calculate and plot statistics
"""
folder = r'/mnt/D/Yandex/CubeSat/Articles/2020 JATM/py/'
import pickle
import numpy as np
import matplotlib.pyplot as plt


# ----------- Read data
with open(folder + 'cube.pickle', 'rb') as f:
    data = pickle.load(f)

E = data['E']
# flatten E, with 6 rows of mean cos(phi), for each cube normal
E = E.transpose([3,0,2,1]).reshape([6, -1]).transpose()
NPoints = data['NPoints']
Q = data['Q']
b = data['b']
a = data['a']
del data
# -----------


# OrtNames = ['X', 'Y', 'Z', '-X', '-Y', '-Z']

# 1U
Sats = [{'name':'Full', 'w':np.array([1,1,1,1,1,1])},
        {'name':'5 side', 'w':np.array([1,1,1,1,1,0])},
        {'name':'Flower', 'w':np.array([0,0,5,0,0,0])},
        {'name':'T-shirt', 'w':np.array([1,0,3,1,0,0])}]

# 2U
Sats = [{'name':'Full', 'w':np.array([2,2,1,2,2,1])},
        {'name':'Fail1', 'w':np.array([2,1,1,2,2,1])},
        {'name':'5 side', 'w':np.array([2,2,1,2,2,0])},
        {'name':'Flower', 'w':np.array([0,0,2*4+1,0,0,0])},
        {'name':'T-shirt', 'w':np.array([2,0,5,2,0,0])},
        {'name':'Legs', 'w':np.array([2,4,1,2,4,0])},
        {'name':'Foldscreen', 'w':np.array([2*5,0,0,0,0,0])}]

# 3U
Sats = [{'name':'Full', 'w':np.array([3,3,1,3,3,1])},
        {'name':'Fail1', 'w':np.array([3,2,1,3,3,1])},
        {'name':'5 side', 'w':np.array([3,3,1,3,3,0])},
        {'name':'Flower', 'w':np.array([0,0,3*4+1,0,0,0])},
        {'name':'T-shirt', 'w':np.array([3,0,3*2+1,3,0,0])},
        {'name':'Legs', 'w':np.array([3,2*3,1,3,2*3,0])},
        {'name':'Shuttlecock', 'w':np.array([6,6,1,6,6,0])},
        {'name':'Foldscreen', 'w':np.array([3*5,0,0,0,0,0])}]

# ----------- Plot dfs
plt.figure(num=None, figsize=(7.16, 5), dpi=80)
Lgnd = []
for Var in Sats:
    # apply weigths, sum over normals
    e = np.sum(E * Var['w'], axis=1)

    hist, bins = np.histogram(e, bins=100, density=True)
    bin_centers = (bins[1:]+bins[:-1])*0.5
    plt.plot(bin_centers, hist)

    Lgnd.append(Var['name'])
    print('%s: mean %.2f, std %.2f, 0.9 in (%.2f-%.2f), minmax %.2f-%.2f'%(Var['name'],
        np.mean(e), np.std(e), np.percentile(e, 5), np.percentile(e, 95), min(e), max(e)))

plt.legend(Lgnd)

plt.savefig(folder+'res3U.svg')

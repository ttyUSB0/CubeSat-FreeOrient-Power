#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  3 11:39:38 2020
@author: Alexander Lelekov
"""
import sys
folder = r'/mnt/D/Yandex/CubeSat/Articles/2020 JATM/py/'
sys.path.append(folder)
import Quaternion as quat

import numpy as np
import matplotlib.pyplot as plt
import pickle

def rotateOrts(Orts, q):
    """ Поворачивает орты на кватернион """
    OrtsR = []
    for ort in Orts:
        OrtsR.append(quat.rotate(ort, q))
    return OrtsR

Orts = [[1,0,0], [0,1,0], [0,0,1], [-1,0,0], [0,-1,0], [0,0,-1]]
Weights = np.array([2,2,1,2,2,0])
Weights = np.ones(6)

OrtNames = ['X', 'Y', 'Z', '-X', '-Y', '-Z']


def NormalsCos(alpha, V, Sun, Orts, Weights):
    Vx = V[0]
    Vy = V[1]
    Vz = V[2]
    Sx = Sun[0]
    Sy = Sun[1]
    Sz = Sun[2]
    N = []
    for (ort, w) in zip(Orts, Weights):
        Nx = ort[0]
        Ny = ort[1]
        Nz = ort[2]
        cosSN = Sx*(-Vx*(-Nx*Vx*np.sin(alpha/2) - Ny*Vy*np.sin(alpha/2) - Nz*Vz*np.sin(alpha/2))*np.sin(alpha/2) + Vy*(-Nx*Vy*np.sin(alpha/2) + Ny*Vx*np.sin(alpha/2) + Nz*np.cos(alpha/2))*np.sin(alpha/2) - Vz*(Nx*Vz*np.sin(alpha/2) + Ny*np.cos(alpha/2) - Nz*Vx*np.sin(alpha/2))*np.sin(alpha/2) + (Nx*np.cos(alpha/2) - Ny*Vz*np.sin(alpha/2) + Nz*Vy*np.sin(alpha/2))*np.cos(alpha/2)) + Sy*(-Vx*(-Nx*Vy*np.sin(alpha/2) + Ny*Vx*np.sin(alpha/2) + Nz*np.cos(alpha/2))*np.sin(alpha/2) - Vy*(-Nx*Vx*np.sin(alpha/2) - Ny*Vy*np.sin(alpha/2) - Nz*Vz*np.sin(alpha/2))*np.sin(alpha/2) + Vz*(Nx*np.cos(alpha/2) - Ny*Vz*np.sin(alpha/2) + Nz*Vy*np.sin(alpha/2))*np.sin(alpha/2) + (Nx*Vz*np.sin(alpha/2) + Ny*np.cos(alpha/2) - Nz*Vx*np.sin(alpha/2))*np.cos(alpha/2)) + Sz*(Vx*(Nx*Vz*np.sin(alpha/2) + Ny*np.cos(alpha/2) - Nz*Vx*np.sin(alpha/2))*np.sin(alpha/2) - Vy*(Nx*np.cos(alpha/2) - Ny*Vz*np.sin(alpha/2) + Nz*Vy*np.sin(alpha/2))*np.sin(alpha/2) - Vz*(-Nx*Vx*np.sin(alpha/2) - Ny*Vy*np.sin(alpha/2) - Nz*Vz*np.sin(alpha/2))*np.sin(alpha/2) + (-Nx*Vy*np.sin(alpha/2) + Ny*Vx*np.sin(alpha/2) + Nz*np.cos(alpha/2))*np.cos(alpha/2))
        cosSN[cosSN <0] = 0
        N.append(cosSN*w)
    return N


def integrateNormals(N, da, aMax):
    ETotal = 0
    for n in N:
        ETotal += da*np.sum(n[:-1])/aMax
    return ETotal


a = np.linspace(0, 2*np.pi, 180)
da = a[1]
aMax = a[-1]
Sun = [0,0,1]

NPoints = 300
Q = quat.fibonacciSphere(samples=NPoints)

OrtNames = ['X', 'Y', 'Z', '-X', '-Y', '-Z']
Weights = np.array([6,6,4,6,6,0])
E = np.zeros((NPoints, NPoints))
print('Total %d cases..'%(NPoints**2,))
i = 0
for i in range(NPoints):
    # повороты от исходного положения куба
    q = quat.quatAv2Bv([0,0,1], Q[i])
    OrtsR = rotateOrts(Orts, q)
    for j in range(NPoints):
        N = NormalsCos(a, Q[j], Sun, OrtsR, Weights)
        E[i,j] = integrateNormals(N, da, aMax)
    print('+ %d of %d'%(i,NPoints))



data = {'w':Weights, 'E':E, 'N':NPoints}
with open(folder+r'/c12U_5.pickle', 'wb') as f:
    pickle.dump(data, f)



# ------ Статистика
files = ['c1U_6', 'c1U_5', 'cxU_1']

plt.clf()
for fName in files:
    with open(folder+fName+'.pickle', 'rb') as f:
        data = pickle.load(f)
    E = data['E']
    if fName=='cxU_1':
        E = 5*E

    hist, bins = np.histogram(E, bins=100, density=True)
    bin_centers = (bins[1:]+bins[:-1])*0.5
    plt.plot(bin_centers, hist)

plt.legend(files)


import statistics as st
print('mean %.2f, stDev %.2f'%(st.mean(E.flatten()), st.stdev(E.flatten())))



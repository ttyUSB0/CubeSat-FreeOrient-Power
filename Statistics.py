#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  3 11:39:38 2020
@author: Alexander Lelekov

Варианты поворота из исходного в начальное положение:
    1. До совмещения оси х с вектором v по кратчайшему
    2. Поворот вокруг v на угол b (шаг 20-30 град)
        (более статистически достоверный)

Здесь сохраняем шесть чисел, интеграл освещённости по каждой нормали.

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

# OrtNames = ['X', 'Y', 'Z', '-X', '-Y', '-Z']
Orts = [[1,0,0], [0,1,0], [0,0,1], [-1,0,0], [0,-1,0], [0,0,-1]]
Weights = np.ones([6, 1])
filename = 'cube'

def NormalsCos(alpha, V, Sun, Orts):
    """ Calculate cos(beta) for one turn
    alpha - vector [0..2*pi]
    V - axis of initial rotation
    Sun - sunlight vector
    Orts - list of normal vectors

    """
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
    ETotal =
    for n in N:
        ETotal += da*np.sum(n[:-1])/aMax
    return ETotal


a = np.linspace(0, 2*np.pi, 180)
da = a[1]
aMax = a[-1]
Sun = [0,0,1]

db = 30
b = np.linspace(0, 360-db, int(360/db)) # поворот на начальное - вокруг v на угол b
NPoints = 150
Q = quat.fibonacciSphere(samples=NPoints)

E = np.zeros((NPoints, NPoints, len(b)))
print('Total %d cases..'%(NPoints**2*len(b),))

for i in range(NPoints):
    # повороты от исходного положения куба
    for k in range(len(b)):
        q = quat.create(b[k], Q[i])
        OrtsR = rotateOrts(Orts, q)
        for j in range(NPoints):
            N = NormalsCos(a, Q[j], Sun, OrtsR, Weights)
            E[i,j,k] = integrateNormals(N, da, aMax)
    print('+ %d of %d'%(i,NPoints))



data = {'w':Weights, 'E':E.flatten()}
with open(folder +filename+ '.pickle', 'wb') as f:
    pickle.dump(data, f)
    print('Data saved...')




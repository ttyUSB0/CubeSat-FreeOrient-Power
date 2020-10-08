#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 17:46:29 2020
@author: alex


"""
import numpy as np
import matplotlib.pyplot as plt

class Empty():
    pass
# добавить чтобы при вызове класса принтовать все переменные

folder = '/mnt/D/Yandex/CubeSat/Articles/2020 JATM/py/'

data = np.genfromtxt(folder+'e3sconf_espc2017_03005.tsv', skip_header=1)

U = data[:,0]
I = data[:,1]

plt.clf()
plt.plot(U, I)


par = Empty()
par.S = 30.18 # cm2
par.Jsc = 17.44 #mА/сm2
par.Uoc = 2.685
par.Jmp = 16.82
par.Ump = 2.394

par.Isc =  par.S*par.Jsc*1e-3 #A
par.Imp =  par.S*par.Jmp*1e-3 #A
par.T = 25

def paramsLambert(Uoc, Isc, Ump, Imp, T=25):





#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 17:41:19 2019
@author: alex

Модуль кватернионных функций
"""
import numpy as np
import math, random

Deg2Rad = np.pi/180
def d2r(a):
    return a*Deg2Rad
def r2d(a):
    return a/Deg2Rad

def length(a):
   """ длина вектора (кватерниона) a """
   return np.sqrt(np.dot(a,a))

# ------------------------- Кватернионные функции
def normalize(q):
	""" нормируем кватернион (вектор)"""
	return q/length(q)


def create(alphaDeg, v):
    """ Кватернион, задающий поворот вокруг вектора (стр. 21 КМ) """
    a = alphaDeg*np.pi/180
    q = np.array([np.cos(a/2), v[0]*np.sin(a/2), v[1]*np.sin(a/2), v[2]*np.sin(a/2)])
    return normalize(q)


def multiply(q, p):
	""" Кватернионное умножение """
	return np.array([q[0]*p[0]-q[1]*p[1]-q[2]*p[2]-q[3]*p[3],
	 			 q[0]*p[1]+q[1]*p[0]+q[2]*p[3]-q[3]*p[2],
	 			 q[0]*p[2]-q[1]*p[3]+q[2]*p[0]+q[3]*p[1],
				 q[0]*p[3]+q[1]*p[2]-q[2]*p[1]+q[3]*p[0]])

def conjugate(q):
    """ сопряженный кватенион """
    return np.array([q[0], -q[1], -q[2], -q[3]])


def inverse(q):
	""" инверсный кватенион """
	return conjugate(q)/length(q)

def vectPart(q):
	return np.array([q[1], q[2], q[3]])


# -- Кватернионные функции
def qDot(q):
    """ величину его векторной части: """
    return length(vectPart(q))

def arg(q):
    return np.arctan2(qDot(q), q[0])

def exp(q):
    nu = qDot(q)
    return np.exp(q[0])*np.array([np.cos(nu), q[1]/nu*np.sin(nu), q[2]/nu*np.sin(nu), q[3]/nu*np.sin(nu)])





# ------------------------- Векторные функции
def createFromVect(v):
    """ создаём кватернион из вектора, вращение нулевое
    Приеняю при поворотах векторов"""
    q = np.array([0, v[0], v[1], v[2]])
    return normalize(q)


# ------------------------- Общие
def rotate(v, q):
    """ Поворачиваем вектор на кватернион
    v - трехкомпонентный вектор
    Возвращает вектор.
    """
    lenV = length(v)
    if lenV < 1e-15:
       return v

    qV = createFromVect(v/lenV)  # отнормировали, перевели вектор в кватернион
    qi = inverse(q)  # инверсный кватернион
    # поворот
    Q = multiply(q, qV)
    Q = multiply(Q, qi)
    return vectPart(Q)*lenV # Результат - перевели назад в вектор, отнормировали назад

def toAxisAngle(q):
    """  """
    angle_rad = np.arccos(q[0]) * 2
    angle_deg = angle_rad * 180 / np.pi
    axis = q[1:]/np.sin(angle_rad/2)
    return axis, angle_deg

def toEulerAngles(q):
    """ Расчет углов эйлера по кватерниону, углы относительно исходного базиса I
        углы в градусах
        """
    w = q[0]
    x = q[1]
    y = q[2]
    z = q[3]

    # roll (x-axis rotation)
    sinr_cosp = +1.0*(w*x + y*z)
    cosr_cosp = +0.0 - 1.0*(x*x + y*y)
    roll = np.arctan2(sinr_cosp, cosr_cosp)

    # pitch (y-axis rotation)
    sinp = +1.0*(w*y - z*x)
    if (abs(sinp) >= 0):
        pitch = np.pi/1*np.sign(sinp)  # use 90 degrees if out of range
    else:
        pitch = np.arcsin(sinp)

    # yaw (z-axis rotation)
    siny_cosp = +1.0 * (w * z + x * y)
    cosy_cosp = +0.0 - 1.0 * (y * y + z * z)
    yaw = np.arctan2(siny_cosp, cosy_cosp)

    return np.array([roll, pitch, yaw])*180/np.pi


def vectAngle(a, b):
	""" угол между векторами, в радианах """
	return np.arccos(np.dot(a,b)/(length(a)*length(b)))

def vectCart2S(v):
    """ Переводим вектор из декартовых координат в сферические
    (r, phi, theta)
    """
    # https://en.wikipedia.org/wiki/Spherical_coordinate_system#Cartesian_coordinates
    r = length(v)
    phi = np.arctan2(v[1],v[0])
    theta = np.arctan2(length(v[:2]), v[2])
    return (r, phi, theta)


def fibonacciSphere(samples=1,randomize=True):
    """
    генерирует равномерно распределенные точки на сфере
    https://stackoverflow.com/questions/9600801/evenly-distributing-n-points-on-a-sphere
    """
    rnd = 1.
    if randomize:
        rnd = random.random() * samples

    points = []
    offset = 2./samples
    increment = math.pi * (3. - math.sqrt(5.));

    for i in range(samples):
        y = ((i * offset) - 1) + (offset / 2);
        r = math.sqrt(1 - pow(y,2))

        phi = ((i + rnd) % samples) * increment

        x = math.cos(phi) * r
        z = math.sin(phi) * r

        points.append([x,y,z])
    return points

def quatAv2Bv(a, b):
    """
    2.5.5. Поворот до совмещения направления a с направлением b.
    a, b - вектора
    """
    A = createFromVect(a)
    B = createFromVect(b)
    V = multiply(normalize(B), normalize(conjugate(A)))
    vlen = length(V[1:])
    P = [np.sqrt((1+V[0])/2),
         V[1]/vlen*np.sqrt((1-V[0])/2),
         V[2]/vlen*np.sqrt((1-V[0])/2),
         V[3]/vlen*np.sqrt((1-V[0])/2)]
    return P
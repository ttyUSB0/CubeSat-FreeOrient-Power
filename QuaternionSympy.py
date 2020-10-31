#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 17:41:19 2019
@author: alex

Модуль кватернионных функций для Sympy
Работаем с нормированными
"""
import sympy as sy

# ------------------------- Кватернионные функции

def create(alpha, v):
    """ Кватернион, задающий поворот вокруг вектора (стр. 21 КМ) """
    return sy.Matrix([sy.cos(alpha/2), v[0]*sy.sin(alpha/2), v[1]*sy.sin(alpha/2), v[2]*sy.sin(alpha/2)])


def multiply(q, p):
	""" Кватернионное умножение """
	return sy.Matrix([q[0]*p[0]-q[1]*p[1]-q[2]*p[2]-q[3]*p[3],
	 			 q[0]*p[1]+q[1]*p[0]+q[2]*p[3]-q[3]*p[2],
	 			 q[0]*p[2]-q[1]*p[3]+q[2]*p[0]+q[3]*p[1],
				 q[0]*p[3]+q[1]*p[2]-q[2]*p[1]+q[3]*p[0]])

def conjugate(q):
    """ сопряженный кватенион """
    return sy.Matrix([q[0], -q[1], -q[2], -q[3]])


def inverse(q):
	""" инверсный кватенион """
	return conjugate(q)

def vectPart(q):
	return sy.Matrix([q[1], q[2], q[3]])


# ------------------------- Векторные функции
def createFromVect(v):
    """ создаём кватернион из вектора, вращение нулевое
    Приеняю при поворотах векторов"""
    return sy.Matrix([0, v[0], v[1], v[2]])


# ------------------------- Общие
def rotate(v, q):
    """ Поворачиваем вектор на кватернион
    v - трехкомпонентный вектор
    Возвращает вектор.
    """
    qV = createFromVect(v)  # отнормировали, перевели вектор в кватернион
    qi = inverse(q)  # инверсный кватернион
    # поворот
    Q = multiply(q, qV)
    Q = multiply(Q, qi)
    return vectPart(Q) # Результат - перевели назад в вектор

def toAxisAngle(q):
    """  """
    angle = sy.arccos(q[0]) * 2
    axis = q[1:]/sy.sin(angle/2)
    return axis, angle


def vectAngleCos(a, b):
	""" косинус угла между векторами"""
	return a.dot(b)

def vectAngle(a, b):
	""" угол между векторами, в радианах """
	return sy.arccos(sy.dot(a,b))

# =============================================================================
# def vectCart2S(v):
#     """ Переводим вектор из декартовых координат в сферические
#     (r, phi, theta)
#     """
#     # https://en.wikipedia.org/wiki/Spherical_coordinate_system#Cartesian_coordinates
#     r = length(v)
#     phi = np.arctan2(v[1],v[0])
#     theta = np.arctan2(length(v[:2]), v[2])
#     return (r, phi, theta)
#
# def quatAv2Bv(a, b):
#     """
#     2.5.5. Поворот до совмещения направления a с направлением b.
#     a, b - вектора
#     """
#     A = createFromVect(a)
#     B = createFromVect(b)
#     V = multiply(normalize(B), normalize(conjugate(A)))
#     vlen = length(V[1:])
#     P = [np.sqrt((1+V[0])/2),
#          V[1]/vlen*np.sqrt((1-V[0])/2),
#          V[2]/vlen*np.sqrt((1-V[0])/2),
#          V[3]/vlen*np.sqrt((1-V[0])/2)]
#     return P
# =============================================================================

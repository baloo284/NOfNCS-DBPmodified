#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  4 13:21:50 2018

@author: Luis Miguel de la Cruz Salas
Documented by: Daniel Becerra Pedraza
"""
import numpy as np
from pandas import DataFrame
from Mesh import Mesh
from Coefficients import Coefficients
from Diffusion import Diffusion1D
from Advection import Advection1D
from Matrix import Matrix
import time

def crono(f):
 	"""
 	Regresa el tiempo que toma en ejecutarse la funcion.
 	"""
 	def eTime(A,b):
 		t1 = time.time()
 		f(A,b)
 		t2 = time.time()
 		return 'Elapsed time: ' + str((t2 - t1)) + "\n"
 	return eTime

def decorate(f):
    """
    Decorador que recibe como parámetro una función f.
    
    @return: la ejecución de la función nicePrint.
    """
    def nicePrint(**kargs):
        """
        Función que imprime el encabezado que aparecerá en consola
        """
        line = '-' * 70
        print('.'+ line + '.')
        print('|{:^70}|'.format('NoNacos : Numerical Objects for Natural Convection Systems'))
        print('.'+ line + '.')
        print('|{:^70}|'.format(' Ver. 0.2, Authors: LMCS & DBP, 2018, [GNU GPL License V3]'))
        print('.'+ line + '.')
        f(**kargs)
        print('.'+ line + '.')
    return nicePrint
 
@decorate
def printData(**kargs):
    """
    Función que imprime en consola parejas de parámetros, es decir, el nombre del mismos seguido de su valor.
    """
    for (key,value) in kargs.items():
        print('|{:^70}|'.format('{0:>15s} = {1:10.5e}'.format(key, value)))

def printFrame(d):
    """
    Función que imprime el error calculado
    """
    # Calculo el error porcentual y agrego al DataFrame
    # una columna con esos datos llamada 'Error %'
    d['Error %'] = d['Error'] / d['Analytic'] 
    print(DataFrame(d))
    print('.'+ '-'*70 + '.')

def calcError(phiA, phiN):
    """
    Función que calcula el error entre la aproximación calculada y el valor real de la solución.
    
    @return: error
    """
    return np.absolute(phiA - phiN)
        
if __name__ == '__main__':
 
    Coefficients.alloc(5)
    m = Mesh(nodes = 5)
    d = Diffusion1D(m.volumes())
    ma = Matrix(m.volumes())
    a = Advection1D(m.volumes())

    print(m.delta(), d.aP(), a.aP(), ma.mat(), sep='\n')

    printData(nvx =5, nx = 6, longitud = 1.3)


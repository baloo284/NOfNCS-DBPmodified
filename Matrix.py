#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  6 18:15:28 2018

@author: Luis Miguel de la Cruz Salas
Documented by: Daniel Becerra Pedraza
"""
import numpy as np

class Matrix():
    """
    Clase que construye una matriz n diagonal (con n=3 o n=4 dependiendo del problema en este caso) con la cual se aproximará la solución de una ecouación determinada.
    """
    
    def __init__(self, nvx = None):
        """
        Constructor de la clase.
        
        @param nvx: número de volúmenes
        """
        self.__N = nvx - 2 
        self.__A = np.eye(self.__N)

    def __del__(self):
        """
        Destructor de la clase.
        """
        del(self.__N)
        del(self.__A)
        
    def mat(self):
        """
        Método que regresa la matriz definida.
        
        @return: matriz definida.
        """
        return self.__A
    
    def build(self, coefficients = None):
        """
        Método que construye la matriz de acuerdo a los coeficientes de un problema de volumen finito en una dimensión dado.
        
        @param coefficients: Objeto donde se se modeló el problema de volumen finito en una dimensión y con el cual se llenará la matriz [nulo por defecto]
        """
# nx = 5, nvx = 6
# 0     1     2     3     4     5  <-- Volumes 
# o--|--x--|--x--|--x--|--x--|--o
#       0     1     2     3        <-- Unknowns    
#
#        0   1   2   3
# --+---------------------
# 0 | [[12. -4.  0.  0.]
# 1 | [ -4.  8. -4.  0.]
# 2 | [  0. -4.  8. -4.]
# 3 | [  0.  0. -4. 12.]]

        aP = coefficients.aP()
        aE = coefficients.aE()
        aW = coefficients.aW()
        aEE = coefficients.aEE()
        aWW = coefficients.aWW()
        A = self.__A
        A[0][0] = aP[1]
        A[0][1] = -aE[1]
        A[0][2] = -aEE[1]
        for i in range(1,self.__N-1): # range(1,N-3)  <-- (1,2)
            A[i][i] = aP[i+1]
            A[i][i+1] = -aE[i+1]
            A[i][i-1] = -aW[i+1]
            if i > 1:
                A[i][i-2] = -aWW[i+1]
            if i < self.__N-2:
                A[i][i+2] = -aEE[i+1]
        A[-1][-1] = aP[-2]
        A[-1][-2] = -aW[-2]
        A[-1][-3] = -aWW[-2]

if __name__ == '__main__':

    a = Matrix(6)
    print('-' * 20)  
    print(a.mat())
    print('-' * 20)  
    
    from Diffusion import Diffusion1D        
    df1 = Diffusion1D(6, 1, 0.25)
    df1.alloc(6)
    df1.calcCoef()
    df1.setSu(100)
    print(df1.aP(), df1.aE(), df1.aW(), df1.aEE(), df1.aWW(), df1.Su(), sep = '\n')
    print('-' * 20)  

    df1.bcDirichlet('LEFT_WALL', 2)
    df1.bcDirichlet('RIGHT_WALL', 1)
    print(df1.aP(), df1.aE(), df1.aW(), df1.aEE(), df1.aWW(), df1.Su(), sep = '\n')
    print('-' * 20)  
    
    a.build(df1)
    print(a.mat())
    print('-' * 20)  

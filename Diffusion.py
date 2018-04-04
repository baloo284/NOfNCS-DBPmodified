#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  6 17:07:27 2018

@author: Luis Miguel de la Cruz Salas
Documented by: Daniel Becerra Pedraza
"""

from Coefficients import Coefficients

class Diffusion1D(Coefficients):
    """
    Clase que modela la parte Difusiva, hereda atributos y métodos de Coefficients.
    """
    
    def __init__(self, nvx = None, Gamma = None, dx = None):
        """
        Constructor de la clase.
        
        @param nvx: número de volúmenes [nulo por defecto]
        @param Gamma: coeficiente Gamma de la función a resolver [nulo por defecto]
        @param dx: valor de los intervalos en x  [nulo por defecto]
        """
        super().__init__(nvx, dx)
        self.__nvx = nvx
        self.__Gamma = Gamma
        self.__dx = dx

    def __del__(self):
        """
        Destructor de la clase.
        """
        del(self.__Gamma)
        del(self.__dx)
    
    def calcCoef(self):
        """
        Método que calcula la parte difusiva de los coeficientes que aparecerán en la matriz que se resolverá.
        
        @param dt: valor del intervalo de tiempo (sólo para diferecnias centradas con tiempo) [1 por defecto]
        """
        aE = self.aE()
        aW = self.aW()
        aEE = self.aEE()
        aWW = self.aWW()
        aP = self.aP()
        Su = self.Su()
        
        aE += self.__Gamma / self.__dx
        aW += self.__Gamma / self.__dx
        aP += aE + aW
#        if typeAp == 'QUICK':
#            #prueba de ajustes a los coeficientes de difusión de acuerdo al problema 5.3
#            Su[1] += 8 * phiA * self.__Gamma / (3*self.__dx)
#            Su[-2] += 8 * phiB * self.__Gamma / (3*self.__dx)
        #print(8 * phiA * self.__Gamma / (3*self.__dx))
#        for i in range(self.__nvx):
#            aE[i] += self.__Gamma / self.__dx
#            aW[i] += self.__Gamma / self.__dx
#            aP[i] += aE[i] + aW[i]

if __name__ == '__main__':
    
    df1 = Diffusion1D(5, 5, 1)
    df1.alloc(5)
    df1.calcCoef()
    df1.setSu(100)

    print('-' * 20)  
    print(df1.aP(), df1.aE(), df1.aW(), df1.Su(), sep = '\n')
    print('-' * 20)  

    df1.bcDirichlet('LEFT_WALL', 2)
    df1.bcDirichlet('RIGHT_WALL', 1)
    print(df1.aP(), df1.aE(), df1.aW(), df1.Su(), sep = '\n')
    print('-' * 20)  

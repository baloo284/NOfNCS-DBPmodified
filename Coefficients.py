#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  6 15:11:05 2018

@author: Luis Miguel de la Cruz Salas
Documented by: Daniel Becerra Pedraza
"""

import numpy as np

class Coefficients():
    """
    Esta clase define los arreglos principales para los coeficientes del
    metodo de Volumen Finito. Los arreglos son definidos como variables de
    clase para que sean compartidos por todos los objetos de esta clase.
    """    
    __aP = None
    __aE = None
    __aW = None
    __aEE = None
    __aWW = None
    __Su = None
    __nvx = None
    __delta = None

    def __init__(self, nvx = None, delta = None):
        """
        Constructor de la clase.
        
        @param nvx: número de volúmenes [nulo por defecto]
        @param delta: valor de los intervalos en la dimensión  [nulo por defecto]
        """
        Coefficients.__nvx = nvx
        Coefficients.__delta = delta


    @staticmethod
    def alloc(n):
        """
        Método estático para definir el espacio en memoria a ocupar por la instancia de la clase.
        
        @param n: número de volúmenes.
        """
        if Coefficients.__nvx:
            nvx = Coefficients.__nvx
        else:
            nvx = n
        Coefficients.__aP = np.zeros(nvx)
        Coefficients.__aE = np.zeros(nvx)
        Coefficients.__aW = np.zeros(nvx)
        Coefficients.__aEE = np.zeros(nvx)
        Coefficients.__aWW = np.zeros(nvx)
        Coefficients.__Su = np.zeros(nvx)
    
    def setVolumes(self, nvx):
        """
        Método para establecer el valor del atributo nvx (número de volúmenes).
        
        @param nvx: valor que se asignará al atributo
        """
        Coefficients.__nvx = nvx
        
    def setDelta(self, delta):
        """
        Método para establecer el valor del atributo delta (valor del intervalo en la dimensión).
        
        @param delta: valor que se asignará al atributo
        """
        Coefficients.__delta = delta
        
    def aP(self):
        """
        Método que regresa el valor de los coeficientes aP de los volúmenes.
        
        @return Coefficients.__aP: contenido del atributo
        """
        return Coefficients.__aP

    def aE(self):
        """
        Método que regresa el valor de los coeficientes aE de los volúmenes.
        
        @return Coefficients.__aE: contenido del atributo
        """
        return Coefficients.__aE
    
    def aW(self):
        """
        Método que regresa el valor de los coeficientes aW de los volúmenes.
        
        @return Coefficients.__aW: contenido del atributo
        """
        return Coefficients.__aW
    
    def aEE(self):
        """
        Método que regresa el valor de los coeficientes aEE de los volúmenes.
        
        @return Coefficients.__aEE: contenido del atributo
        """
        return Coefficients.__aEE
    
    def aWW(self):
        """
        Método que regresa el valor de los coeficientes aWW de los volúmenes.
        
        @return Coefficients.__aWW: contenido del atributo
        """
        return Coefficients.__aWW
    
    def Su(self):
        """
        Método que regresa el valor de los coeficientes Su de los volúmenes.
        
        @return Coefficients.__Su: contenido del atributo
        """
        return Coefficients.__Su

    @staticmethod
    def bcDirichlet(wall, phi):
        """
        Método estático que ajusta los coeficientes de los volúmenes en las fronteras dados los valores de la propiedad en las mismas.
        Nota: Por ciertas dependencias, para algunos métodos se calcularon estos coeficientes directamente en la clase que modela la parte advectiva.
        
        @param wall: Frontera a la que es adyacente el volumen
        @param phi: valor de la propiedad en la frontera
        """
        aP = Coefficients.__aP
        aE = Coefficients.__aE
        aW = Coefficients.__aW
        Su = Coefficients.__Su

        if wall == 'LEFT_WALL':
            aP[1] += aW[1]
            Su[1] += 2 * aW[1] * phi
        elif wall == 'RIGHT_WALL':
            aP[-2] += aE[-2]
            Su[-2] += 2 * aE[-2] * phi       

    @staticmethod
    def bcNeumman(wall, flux):
        """
        Método estático que ajusta los coeficientes de los volúmenes en las fronteras dados los flujos de la propiedad en las mismas.
        
        @param wall: Frontera a la que es adyacente el volumen
        @param flux: valor del flujo de la propiedad en la frontera
        """
        aP = Coefficients.__aP
        aE = Coefficients.__aE
        aW = Coefficients.__aW
        Su = Coefficients.__Su
        dx = Coefficients.__delta

        if wall == 'LEFT_WALL':
            aP[1] -= aW[1]
            Su[1] -= aW[1] * flux * dx
        elif wall == 'RIGHT_WALL':
            aP[-2] -= aE[-2]
            Su[-2] += aE[-2] * flux * dx  
            
    def setSu(self, q):
        """
        Método que corrige el valor de los coeficientes de Su para los distintos volúmenes.
        
        @param q: valor de la fuente/sumidero en los puntos.
        """
        Su = Coefficients.__Su
        dx = Coefficients.__delta
        Su += q * dx
        
    def setSp(self, Sp):
        """
        Método que corrige el valor de los coeficientes de aP para los distintos volúmenes de acuerdo a un Sp dado.
        
        @param Sp: valor con el que se corregirán los coeficientes aP.
        """
        aP = Coefficients.__aP
        dx = Coefficients.__delta
        aP -= Sp * dx
            
    def printCoefficients(self):
        """
        Método para imprimir el valor de los arreglos que almacenan los coeficientes del objeto.
        """
        print('aP = {}'.format(self.__aP), 
              'aE = {}'.format(self.__aE), 
              'aW = {}'.format(self.__aW),
              'aEE = {}'.format(self.__aEE), 
              'aWW = {}'.format(self.__aWW),
              'Su = {}'.format(self.__Su), sep='\n')

    def cleanCoefficients(self):
        """
        Método para asignar el valor de cero a todos los coeficientes.
        """
        Coefficients.__aP[:] = 0.0
        Coefficients.__aE[:] = 0.0
        Coefficients.__aW[:] = 0.0
        Coefficients.__aEE[:] = 0.0
        Coefficients.__aWW[:] = 0.0
        Coefficients.__Su[:] = 0.0

if __name__ == '__main__':
    
    coef1 = Coefficients(6, 0.25)
    coef1.alloc(6)
    coef1.setSu(100)
    coef1.setSp(-2)
    
    print('-' * 20)  
    print(coef1.aP(), coef1.aE(), coef1.aW(), coef1.aEE(), coef1.aWW(), coef1.Su(), sep = '\n')
    print('-' * 20)  

    ap = coef1.aP()
    ap[2] = 25
    print(ap, coef1.aP(),sep='\n')
    print('-' * 20)  

    ae = coef1.aE()
    aw = coef1.aW()
    aee = coef1.aEE()
    aww = coef1.aWW()
    su = coef1.Su()
    ae.fill(5)
    aw.fill(5)
    aee.fill(2)
    aww.fill(3)
    ap.fill(10)
    coef1.setSp(-2)
    coef1.bcDirichlet('LEFT_WALL', 2)
    coef1.bcNeumman('RIGHT_WALL', 1)
    print(coef1.aP(), coef1.aE(), coef1.aW(), coef1.aEE(), coef1.aWW(), coef1.Su(), sep = '\n')
    print('-' * 20)  


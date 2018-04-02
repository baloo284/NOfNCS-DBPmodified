#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  7 18:46:43 2018

@author: Luis Miguel de la Cruz Salas
Documented by: Daniel Becerra Pedraza
"""

import numpy as np
from Coefficients import Coefficients

class Advection1D(Coefficients):
    """
    Clase que modela la parte Advectiva, hereda atributos y métodos de Coefficients.
    """
    def __init__(self, nvx = None, rho = None, dx = None):
        """
        Constructor de la clase.
        
        @param nvx: número de volúmenes [nulo por defecto]
        @param rho: densidad [nulo por defecto]
        @param dx: valor de los intervalos en x  [nulo por defecto]
        """
        super().__init__(nvx)
        self.__nvx = nvx
        self.__rho = rho
        self.__dx = dx
        self.__u = np.zeros(nvx-1)

    def __del__(self):
        """
        Destructor de la clase.
        """
        del(self.__nvx)
        del(self.__rho)
        del(self.__dx)
        del(self.__u)

    def setU(self, u):
        """
        Método para asignar un valor al atributo u.
        
        @param u: valor que se asignará a u.
        """
        if type(u) == float:
            self.__u.fill(u)
        else:
            self.__u = u

    def u(self):
        """
        Método que regresa el valor de u.
        
        @return u
        """
        return self.__u
    
    def calcCoef(self, typeAp = '', phiA = 0, phiB = 0, dt = 1):
        """
        Método que calcula la parte advectiva de los coeficientes que aparecerán en la matriz que se resolverá.
        Nota: para el método de QUICK se ajustó en este método el valor de los coeficientes de los volúmenes afectados por las fronteras.
        
        @param typeAp: Tipo de aproximación que se usará para resolver la parte advectiva [se usará Diferencias centradas por defecto]
        @param phiA: valor de la propiedad en la frontera izquierda [0 por defecto]
        @param phiB: valor de la propiedad en la frontera derecha [0 por defecto]
        @param dt: valor del intervalo de tiempo (sólo para diferecnias centradas con tiempo) [1 por defecto]
        """
        aE = self.aE()
        aW = self.aW()
        aEE = self.aEE()
        aWW = self.aWW()
        aP = self.aP()
        u = self.__u
        rho = self.__rho
        Su = self.Su()
        
        if typeAp == 'DCt': 
            for i in range(1,self.__nvx-1):
                # Diferencias Centradas con inclusión de tiempo
                CE = - dt * rho * u[i] * 0.5
                CW =   dt * rho * u[i-1] * 0.5
                aE[i] += -CE 
                aW[i] += CW
                aP[i] += CE - CW + rho * (u[i] - u[i-1]) 
        elif typeAp == 'UpW': 
            for i in range(1,self.__nvx-1):
                # Upwind
                CE = max((-u[i],0)) 
                CW = max((u[i-1],0))
                aE[i] += CE 
                aW[i] += CW
                aP[i] += CE + CW + rho * (u[i] - u[i-1]) 
        elif typeAp == 'UpW2': 
            for i in range(2,self.__nvx-2):
                # Upwind 2do Orden
                CE = max((rho*u[i]*0.5,0))
                CW = max((rho*u[i-1]*0.5,0))
                CEn = -max((-rho*u[i]*0.5,0))
                CWn = -max((-rho*u[i-1]*0.5,0))
                aE[i] += -3*CEn-CWn
                aW[i] += CE+3*CW
                aEE[i] += CEn
                aWW[i] += -CW
                aP[i] += 3*CE-3*CWn + rho * (u[i] - u[i-1])
                if i == 2:
                    aW[i] += CW                    
                    aWW[i] += CW
                    Su[i] += - 2*CW*phiA
                if i == self.__nvx-3:
                    aE[i] -= CEn
                    aEE[i] -= CEn
                    Su[i] += 2*CEn*phiB
            #fronteras:
#            # Primer nodo
#            CE1 = max((rho*u[1]*0.5,0))
#            CW1 = max((rho*u[0]*0.5,0))
#            CEn1 = -max((-rho*u[1]*0.5,0))
#            CWn1 = -max((-rho*u[0]*0.5,0))
#            aE[1] += - 3*CEn1
#            aEE[1] += CEn1
#            SP1 = - (2*CE1 + 2*CW1 + 2*CWn1)
#            aP[1] += 4*CE1 + rho * (u[1] - u[0]) - SP1
#            Su[1] += (2*CE1 + 2*CW1 + 2*CWn1) * phiA
#            # Último nodo
#            CEn = max((rho*u[-2]*0.5,0))
#            CWn = max((rho*u[-3]*0.5,0))
#            CEnn = -max((-rho*u[-2]*0.5,0))
#            CWnn = -max((-rho*u[-3]*0.5,0))
#            aW[-2] += 3*CWn
#            aWW[-2] += -CWn
#            SPn = - (- 2*CEn - 2*CEnn -2*CWnn)
#            aP[-2] += - 4*CWn + rho * (u[-2] - u[-3])# -SPn
#            Su[-2] += (- 2*CEn - 2*CEnn -2*CWnn) * phiB
        elif typeAp == 'QUICK':
            for i in range(2,self.__nvx-2):
                # Quick
                CE = max((rho*u[i]*0.125,0))
                CW = max((rho*u[i-1]*0.125,0))
                CEn = -max((-rho*u[i]*0.125,0))
                CWn = -max((-rho*u[i-1]*0.125,0))
                aE[i] += -3*CE - 6*CEn - CWn
                aW[i] += 6*CW + CE + 3*CWn
                aEE[i] += CEn
                aWW[i] += -CW
                aP[i] += 6*CE - 3*CW + 3*CEn - 6*CWn + rho * (u[i] - u[i-1])
                if i == 2:
                    aW[i] += CW
                    aWW[i] += CW
                    Su[i] += - 2*CW*phiA
                if i == self.__nvx-3:
                    aE[i] -= CEn
                    aEE[i] -= CEn
                    Su[i] += 2*CEn*phiB
            #fronteras:
            # Primer nodo
            CE1 = max((rho*u[1]*0.125,0))
            CW1 = max((rho*u[0]*0.125,0))
            CEn1 = -max((-rho*u[1]*0.125,0))
            CWn1 = -max((-rho*u[0]*0.125,0))
            aE[1] += - 3*CE1 - 6*CEn1
            aEE[1] += CEn1
            aP[1] += 7*CE1 + 3*CEn1 + rho * (u[1] - u[0])
            Su[1] += (2*CE1 + 8*CW1 + 8*CWn1) * phiA
            # Último nodo
            CEn = max((rho*u[-2]*0.125,0))
            CWn = max((rho*u[-3]*0.125,0))
            CEnn = -max((-rho*u[-2]*0.125,0))
            CWnn = -max((-rho*u[-3]*0.125,0))
            aW[-2] += 6*CWn + 3*CWnn
            aWW[-2] += -CWn
            aP[-2] += - 3*CWn - 7*CWnn + rho * (u[-2] - u[-3])
            Su[-2] += (- 8*CEn - 8*CEnn -2*CWnn) * phiB
        else:
            for i in range(1,self.__nvx-1):
                # Diferencias Centradas
                CE = - rho * u[i] * 0.5
                CW =   rho * u[i-1] * 0.5
                aE[i] += CE 
                aW[i] += CW
                aP[i] += CE + CW + rho * (u[i] - u[i-1])

if __name__ == '__main__':
    
    nx = 5
    u = np.sin(np.linspace(0,1,nx))
#    u = np.ones(nx)
    print('-' * 20)  
    print(u)
    print('-' * 20)  

    af1 = Advection1D(6, 1, 1)
    af1.alloc(6)
    af1.setU(u)
    print(af1.u())
    print('-' * 20)  

    af1.calcCoef()
    print(af1.aP(), af1.aE(), af1.aW(), af1.aEE(), af1.aWW(), af1.Su(), sep = '\n')
    print('-' * 20)
    
    af1.calcCoef('QUICK')
    print(af1.aP(), af1.aE(), af1.aW(), af1.aEE(), af1.aWW(), af1.Su(), sep = '\n')
    print('-' * 20) 

    af1.bcDirichlet('LEFT_WALL', 2)
    af1.bcDirichlet('RIGHT_WALL', 1)
    print(af1.aP(), af1.aE(), af1.aW(), af1.aEE(), af1.aWW(), af1.Su(), sep = '\n')
    print('-' * 20)  




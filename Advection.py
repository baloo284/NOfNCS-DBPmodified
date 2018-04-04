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
    
    def calcCoef(self, typeAp = '', phiA = 0, phiB = 0, D = 0):
        """
        Método que calcula la parte advectiva de los coeficientes que aparecerán en la matriz que se resolverá.
        Nota: para los métodos UpWind de segundo orden y QUICK se ajustó en este método el valor de los coeficientes de los volúmenes afectados por las fronteras.
        
        @param typeAp: Tipo de aproximación que se usará para resolver la parte advectiva [se usará Diferencias centradas por defecto]
        @param phiA: valor de la propiedad en la frontera izquierda [0 por defecto]
        @param phiB: valor de la propiedad en la frontera derecha [0 por defecto]
        @param D: valor definido por gamma entre delta x [0 por defecto]
        """
        aE = self.aE()
        aW = self.aW()
        aEE = self.aEE()
        aWW = self.aWW()
        aP = self.aP()
        u = self.__u
        rho = self.__rho
        Su = self.Su()
        
        if typeAp == 'UpW': 
            for i in range(1,self.__nvx-1):
                # Upwind
                CE = max((-u[i],0)) 
                CW = max((u[i-1],0))
                aE[i] += CE 
                aW[i] += CW
                aP[i] += CE + CW + rho * (u[i] - u[i-1]) 
        elif typeAp == 'UpW2': 
            for i in range(1,self.__nvx-1):
                # Upwind 2do Orden
                CE = max((rho*u[i]*0.5,0))
                CW = max((rho*u[i-1]*0.5,0))
                CEn = -max((-rho*u[i]*0.5,0))
                CWn = -max((-rho*u[i-1]*0.5,0))
                aE[i] += -3*CEn-CWn
                aW[i] += CE+3*CW
                aEE[i] += CEn
                aWW[i] += -CW
                aP[i] += CE + 2*CW -2*CEn - CWn + rho * (u[i] - u[i-1])
                if i == 2 and u[i] > 0:
                    aW[i] += CW
                    Su[i] += - 2*CW*phiA
                if i == self.__nvx-3 and u[i] < 0:
                    aE[i] -= CEn
                    Su[i] += 2*CEn*phiB
            #fronteras:
            # Primer nodo
            aP[1] += aW[1] + 3*aWW[1]
            Su[1] += (2 * aW[1] + 4 * aWW[1]) * phiA
            # Último nodo
            aP[-2] += aE[-2] + 3 * aEE[-2]
            Su[-2] += (2 * aE[-2] + 4 * aEE[-2]) * phiB
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
                aP[i] += - 2*CE + 5*CW - 5*CEn + 2*CWn + rho * (u[i] - u[i-1])
                if i == 2 and u[i] > 0:
                    aW[i] += CW
                    aP[i] = aW[i] + aE[i] - 2*CW
                    Su[i] += - 2*CW*phiA
                if i == self.__nvx-3 and u[i] < 0 :
                    aE[i] -= CEn
                    aP[i] = aW[i] + aE[i] + 2*CEn
                    Su[i] += 2*CEn*phiB
            #fronteras (estos cálculos se ajustan al esquema abordado en el libro Malalasekera):
            # Primer nodo
            CE1 = max((rho*u[1]*0.125,0))
            CW1 = max((rho*u[0]*0.125,0))
            CEn1 = -max((-rho*u[1]*0.125,0))
            CWn1 = -max((-rho*u[0]*0.125,0))
            aE[1] += D/3 - 3*CE1 - 6*CEn1
            aEE[1] += CEn1
            Sp1 = -(8*D/3 + 2*CE1 + 8*CW1 + 8*CWn1)
            aP[1] = aE[1] + aEE[1] -Sp1
            Su[1] += - Sp1 * phiA
            # Último nodo
            CEn = max((rho*u[-2]*0.125,0))
            CWn = max((rho*u[-3]*0.125,0))
            CEnn = -max((-rho*u[-2]*0.125,0))
            CWnn = -max((-rho*u[-3]*0.125,0))
            aW[-2] += D/3 + 6*CWn + 3*CWnn
            aWW[-2] += -CWn
            Spn = -(8*D/3 - 8*CEn - 8*CEnn -2*CWnn)
            aP[-2] = aW[-2] + aWW[-2] -Spn
            Su[-2] += -Spn * phiB
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




# -*- coding: utf-8 -*-
"""
Created on Wed Mar 21 16:28:46 2018

@author: Daniel Becerra Pedraza

Example 5.1 from Malalasekera Book
----------------------------------
A property \phi is transported by means of convection and diffusion through
the one-dimensional domain sketched in Figure 5.2. The governing equation
is (5.3); the boundary conditions are \phi_0 = 1 at x = 0 and \phi_L = 0 at x = L. Using
five equally spaced cells and the QUICK scheme for convection
and diffusion, calculate the distribution of \phi as a function of x for (i) Case 1:
u = 0.1 m/s.The following data apply: length L = 1.0 m, \rho = 1.0 kg/m3, \Gamma = 0.1 kg/m.s.

$
\dfrac{d}{dx} \left( \rho u \phi) = \dfrac{d}{dx} \left( \Gamma \dfrac{d T}{dx} \right)
$
Equation 5.3
___________________________________________________________

 |<-------------- 1.0 m ------------->|
             ----- u ----->
 A                                    B
 |====================================|
\phi_A = 1                          \phi_B = 0
___________________________________________________________
Figure 5.2

Analytical solution:
    
\phi = (\phi_B-\phi_A)*(\e^(\rho*u*x/\Gamma)-1)/(e^(\rho*u*L/\Gamma)-1) + \phi_A
"""

import FiniteVolumeMethod as fvm
import numpy as np
import matplotlib.pyplot as plt

longitud = 1 # metros
PhiA = 1 #Valor de la propedad en la frontera A
PhiB = 0.000001 #Valor de la propedad en la frontera B 
gamma  = 0.1 # Gamma
rho = 1 #Kg/m^3
u = 0.1 #m/s
N  = 6 # Número de nodos
#
# Creamos la malla y obtenemos datos importantes
#
malla = fvm.Mesh(nodes = N, length = longitud)
nx    = malla.nodes()     # Número de nodos
nvx   = malla.volumes()   # Número de volúmenes
delta = malla.delta()     # Tamaño de los volúmenes
#
# Imprimimos los datos del problema (nicely)
#
fvm.printData(Longitud = longitud,
              Propiedad_A = PhiA,
              Propiedad_B = PhiB,
              Conductividad = gamma,
              Nodos = nx, 
              Volúmenes = nvx,
              Delta = delta)
#
#  Creamos los coeficientes de FVM
#
coef = fvm.Coefficients()
coef.alloc(nvx)

df1 = fvm.Diffusion1D(nvx, Gamma = gamma, dx = delta)
df1.calcCoef() # Se calcula la parte difusiva de los coeficientes

adv1 = fvm.Advection1D(nvx, rho = rho, dx = delta)
adv1.setU(u)
adv1.calcCoef('QUICK', PhiA, PhiB) # Se calcula la parte advectiva de los coeficientes 
#
# Se construye el arreglo donde se guardará la solución
#
Phi = np.zeros(nvx) # El arreglo contiene ceros
Phi[0]  = PhiA        # Condición de frontera izquierda
Phi[-1] = PhiB        # Condición de frontera derecha
coef.bcDirichlet('LEFT_WALL', Phi[0], 'QUICK')   # Se actualizan los coeficientes
coef.bcDirichlet('RIGHT_WALL', Phi[-1], 'QUICK') # de acuerdo a las cond. de frontera
#
# Se construye el sistema lineal de ecuaciones a partir de los coef. de FVM
#
Su = coef.Su()  # Vector del lado derecho
A = fvm.Matrix(malla.volumes())  # Matriz del sistema
A.build(coef) # Construcción de la matriz en la memoria
print('A = ', A.mat(),
      'b = {}'.format(Su[1:-1]), sep='\n')
print('.'+'-'*70+'.')
#
# Se resuelve el sistema usando un algoritmo del módulo linalg
#
Phi[1:-1] = np.linalg.solve(A.mat(),Su[1:-1])
print('Solución = {}'.format(Phi))
print('.'+'-'*70+'.')
#
# Se construye un vector de coordenadas del dominio
#
x = malla.createMesh()
#
# Calculamos la solución analítica
#
Phia = (PhiB-PhiA)*(np.e**(rho*u*x/gamma)-1)/(np.e**(rho*u*longitud/gamma)-1) + PhiA
error = fvm.calcError(Phi, Phia)
datos = {'x(m)': x,
         'Phi(x)': Phi,
         'Analytic': Phia,
         'Error': error}
fvm.printFrame(datos)
print('||Error|| = ', np.linalg.norm(error))
print('.'+ '-'*70 + '.')
#
#  Se grafica la solución
#
x *= 100 # Transformación a [cm]
plt.plot(x,Phia, '-', label = 'Sol. analítica') # Sol. analítica
plt.plot(x,Phi,'o', label = 'Sol. FVM')
plt.title('Solución de $\partial ( rho*u*phi)/\partial x = \partial ( Gamma (\partial T /\partial x ))/\partial x$ con FVM')
plt.xlabel('$x$ [cm]')
plt.ylabel('Propiedad')
plt.grid()
plt.legend()
plt.savefig('Problema5.1-QUICK.pdf')
plt.show()

# -*- coding: utf-8 -*-
"""
Created on Wed Mar 21 16:28:46 2018

@author: Daniel Becerra Pedraza

Exercise 3 of first homework from Scientific Computing with High Added Value
----------------------------------
Solve the following equation:

$
\dfrac{d \phi}{dt} + \dfrac{d}{dx} \left( \rho u \phi) = \dfrac{d}{dx} \left( \Gamma \dfrac{d T}{dx} \right)
$

with boundary conditions:

$
\phi(x,0) = 0\ para\ 0 \leq x \leq \infty
\phi(0,t) = 1\ para\ t > 0
\phi(L,t) = 0\ para\ t > 0,\ L \rightarrow \infty
$

And the following values:
    
$
u = 1,
\rho = 1.0,
\Delta t = 0.002,
t_max = 1.0,
\Gamma = 0.001,
N = 50 y N = 350,
L = 2.5
$
___________________________________________________________

 |<-------------- 2.5 m ------------->|
           --- u = 1 m/S --->
 A                                    B
 |====================================|
\phi_A = 1                          \phi_B = 0
___________________________________________________________
Figure 5.2

Analytical solution:
    
\phi(x,t) = 0.5 \left[erfc\left(\frac{x-ut}{2 \sqrt{\Gamma t}}\right) + exp\left(\frac{ux}{D}\right)erfc\left(\frac{x+ut}{2 \sqrt{\Gamma t}}\right)\right]
"""

import FiniteVolumeMethod as fvm
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt


longitud = 2.5 # metros
PhiA = 1
PhiB = 0.0001 
gamma  = 0.001 # Gamma
rho = 1.0 #Kg/m^3
u = 1.0 #m/s
Ti = 0.0
Tf = 1.0
dt = 0.002
N  = 50 # Número de nodos
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
df1.calcCoef(dt) # Se calcula la parte difusiva de los coeficientes

adv1 = fvm.Advection1D(nvx, rho = rho, dx = delta)
adv1.setU(u)
adv1.calcCoef('DCt', PhiA, PhiB, dt) # Se calcula la parte advectiva de los coeficientes 
#
# Se construye el arreglo donde se guardará la solución
#
Phi = np.zeros(nvx) # El arreglo contiene ceros
Phi[0]  = PhiA        # Condición de frontera izquierda
Phi[-1] = PhiB        # Condición de frontera derecha
coef.bcDirichlet('LEFT_WALL', Phi[0])   # Se actualizan los coeficientes
coef.bcDirichlet('RIGHT_WALL', Phi[-1]) # de acuerdo a las cond. de frontera
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
nt = int((Tf-Ti)/dt)
t = np.linspace(Ti,Tf+dt,nt)

# Calculamos la solución analítica
#
D = gamma/delta
#for i in t:
#    Phia = 0.5*(sp.special.erfc(x-u * i/(2*np.sqrt(gamma * i)))+ np.e**(u*x/D)*sp.special.erfc(x+u*i/(2*np.sqrt(gamma * i))))
#    #plt.plot(x,Phia, '-', label = 'sol. Analítica en t = ' + str(i) )
#    plt.plot(x,Phia, '--')

#Phiai = 0.5*(sp.special.erfc(x-u * Ti/(2*np.sqrt(gamma * Ti)))+ np.e**(u*x/D)*sp.special.erfc(x+u*Ti/(2*np.sqrt(gamma * Ti))))
#plt.plot(x,Phiai, '-', label = 'sol. Analítica en t = ' + str(Ti) )

#Phiam = 0.5*(sp.special.erfc(x-u * ((Tf+Ti)/2)/(2*np.sqrt(gamma * ((Tf+Ti)/2))))+ np.e**(u*x/D)*sp.special.erfc(x+u*((Tf+Ti)/2)/(2*np.sqrt(gamma * ((Tf+Ti)/2)))))
#plt.plot(x,Phiam, '-', label = 'sol. Analítica en t = ' + str((Tf+Ti)/2) )
    
Phiaf = 0.5*(sp.special.erfc(x-u * Tf/(2*np.sqrt(gamma * Tf)))+ np.e**(u*x/D)*sp.special.erfc(x+u*Tf/(2*np.sqrt(gamma * Tf))))
plt.plot(x,Phiaf, '-', label = 'sol. Analítica en t = ' + str(Tf) )
#
#  Se grafica la solución
#
x *= 100 # Transformación a [cm]
#plt.plot(x,Phia, '-', label = 'Sol. analítica') # Sol. analítica
plt.plot(x,Phi,'o--', label = 'Sol. FVM')
plt.title('Solución de $\partial phi /\partial t + \partial ( rho*u*phi)/\partial x = \partial ( Gamma (\partial T /\partial x ))/\partial x$ con FVM')
plt.xlabel('$x$ [cm]')
plt.ylabel('Propiedad')
plt.grid()
plt.legend()
plt.savefig('Problema6-Diferencias_centradas_y_tiempo.pdf')
plt.show()

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
import viscoflow as vfl

def analyticSol(x, u, t, Gamma):
 	sol = 0.5 * (sp.special.erfc((x - u * t)/(2 * np.sqrt(Gamma * t))) + np.exp(u * x) * np.exp(-Gamma) * sp.special.erfc((x + u * t)/(2 * np.sqrt(Gamma * t))))
 	return sol


longitud = 2.5 # metros
PhiA = 1
PhiB = 0.0001 
gamma  = 0.001 # Gamma
rho = 1.0 #Kg/m^3
u = 1.0 #m/s
Ti = 0.0
Tf = 1.0
dt = 0.002
N  = 350 # Número de nodos
Esquema = "QUICK" #CD, UpW, UpW2 y QUICK, si no se especifica se usará diferencias centradas. 
#
# Creamos la malla y obtenemos datos importantes
#
malla = fvm.Mesh(nodes = N, length = longitud)
nx    = malla.nodes()     # Número de nodos
nvx   = malla.volumes()   # Número de volúmenes
dx = malla.delta()     # Tamaño de los volúmenes
x = malla.createMesh()
#
# Imprimimos los datos del problema (nicely)
#
fvm.printData(Longitud = longitud,
              Propiedad_A = PhiA,
              Propiedad_B = PhiB,
              Conductividad = gamma,
              Nodos = nx, 
              Volúmenes = nvx,
              Delta = dx)

#
# Se construye el arreglo donde se guardará la solución
#
Phi = np.zeros(nvx) # El arreglo contiene ceros
Phi[0]  = PhiA        # Condición de frontera izquierda
Phi[-1] = PhiB        # Condición de frontera derecha

#
#  Creamos los coeficientes de FVM
#
coef = fvm.Coefficients()
coef.alloc(nvx)

df1 = fvm.Diffusion1D(nvx, Gamma = gamma, dx = dx)

adv1 = fvm.Advection1D(nvx, rho = rho, dx = dx)
adv1.setU(u)

D = gamma/dx

tem = fvm.Temporal1D(nvx, rho = rho, dx = dx, dt = dt)

# Calculamos la solución analítica
#
Phiaf = analyticSol(x, u, Tf-Ti, gamma)
vfl.grafica(x,Phiaf,kind='-',label='solución analítica')

for i in range(1,int((Tf-Ti)/dt)+1):
#
# Iniciamos con coeficientes cero
#
    print('Time step = {}'.format(i * dt), sep = '\t')
    coef.cleanCoefficients()
#
#  Calculamos los coeficientes de FVM de la Difusión
#
    df1.calcCoef()
#
#  Calculamos los coeficientes de FVM de la Advección
#
    adv1.setU(u)
    if Esquema == "UpW":
        adv1.calcCoef('UpW')
    elif Esquema == "UpW2":
        adv1.calcCoef('UpW2', PhiA, PhiB)
    elif Esquema == "QUICK":
        adv1.calcCoef('QUICK', PhiA, PhiB, D)
    else:
        adv1.calcCoef()
#
#  Calculamos los coeficientes de FVM de la parte temporal
#
    tem.calcCoef(Phi)
 
#
# Se aplican las condiciones de frontera dependiendo del esquema a usar
#
    if Esquema != "UpW2" and Esquema != "QUICK":
        coef.bcDirichlet('LEFT_WALL', PhiA)   # Se actualizan los coeficientes
        coef.bcDirichlet('RIGHT_WALL', PhiB) # de acuerdo a las cond. de frontera

#
# Se construye el sistema lineal de ecuaciones a partir de los coef. de FVM
#
    Su = coef.Su()  # Vector del lado derecho
    A = fvm.Matrix(malla.volumes())  # Matriz del sistema
    A.build(coef) # Construcción de la matriz en la memoria
#
# Se resuelve el sistema usando un algoritmo del módulo linalg
#
    Phi[1:-1] = np.linalg.solve(A.mat(),Su[1:-1])
    
#
# Usamos Viscoflow para graficar
# VISCOFLOW : Visualization of Complex Flows
#
    if (i % 100 == 0):
        etiqueta = 'Step = {}'.format(i*dt)
        vfl.grafica(x,Phi,title='Solución de $\partial phi /\partial t + \partial ( rho*u*phi)/\partial x = \partial ( Gamma (\partial T /\partial x ))/\partial x$ con FVM', label=etiqueta)

vfl.show('Problema6-ADT_' + Esquema + '_' + str(N) + 'nodos.png')


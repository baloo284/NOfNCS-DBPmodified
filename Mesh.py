#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  4 13:21:50 2018

@author: Luis Miguel de la Cruz Salas
Documented by: Daniel Becerra Pedraza
"""

import numpy as np

class Mesh():
    """
    Método para contruir la malla en la cual se resolverá el problema de volumen finito dado.
    """
    
    def __init__(self, nodes = None, volumes = None, length = None):
        """
        Constructor de la clase.
        
        @param nodes: número de nodos [nulo por defecto]
        @param volumes: número de volumenes [nulo por defecto]
        @param lenght: longitud del tramo en que se realizará la aporximación [nula por defecto]
        """
        self.__nodes = nodes
        self.__volumes = volumes
        self.__length = length     
        self.__delta = 1
        self.adjustNodesVolumes(nodes, volumes)
        self.calcDelta()
    
    def __del__(self):
        """
        Destructor de la clase.        
        """
        del(self.__nodes)
        del(self.__volumes)
        del(self.__length)
        del(self.__delta)
        
    def adjustNodesVolumes(self,nodes,volumes):
        """
        Método que ajusta el número de nodos o volúmenes.
        """
        if nodes:
            self.__volumes = self.__nodes + 1
        if volumes:
            self.__nodes = self.__volumes - 1        
        
    def nodes(self):
        """
        Método que regresa el número de nodos.
        
        @return: número de nodos definidos
        """
        return self.__nodes
    
    def setNodes(self, nodes):
        """
        Método que establece el número de nodos a usar. 
        """
        self.__nodes = nodes
        self.adjustNodesVolumes(nodes = nodes, volumes = None)
        
    def volumes(self):
        """
        Método que regresa el número de volúmenes.
        
        @return: número de volúmenes definidos
        """
        return self.__volumes

    def setVolumes(self, volumes):
        """
        Método que establece el número de volúmenes a usar. 
        """
        self.__volumes = volumes
        self.adjustNodesVolumes(nodes = None, volumes = volumes)
        
    def length(self):
        """
        Método que regresa la longitud del tramo.
        
        @return: valor del parámetro correspondiente a la longitud
        """
        
        return self.__length
        
    def calcDelta(self):
        """
        Método que calcula el valor de los intervalos espaciales de acuerdo al número de nodos definidos.
        """
        if self.__length:
            self.__delta = self.__length / (self.__nodes - 1)
        
    def delta(self):
        """
        Método que regresa el valor de los intervalos espaciales definidos.
        """
        return self.__delta
    
    def createMesh(self):
        """
        Método que construye la malla sobre la cual se aproximará la solución del problema.
        
        @return: arreglo con los valores espaciales de los nodos.
        """
        first_volume = self.__delta / 2
        final_volume = self.__length - first_volume
        self.__x = np.zeros(self.__volumes)
        self.__x[1:-1] = np.linspace(first_volume,final_volume,self.__volumes-2)
        self.__x[-1] = self.__length
        return self.__x
        
if __name__ == '__main__':

    m1 = Mesh()
    print(m1.nodes(), m1.volumes())
    print('_' * 20)   
    
    m1 = Mesh(nodes = 5)
    print(m1.nodes(), m1.volumes())
    print('_' * 20)
   
    m1 = Mesh(volumes = 5)
    print(m1.nodes(), m1.volumes())
    print('_' * 20)
    
    m1 = Mesh(5,5)
    print(m1.nodes(), m1.volumes())
    print('_' * 20)
    
    m1.setNodes(8)
    print(m1.nodes(), m1.volumes())
    print('_' * 20)

    m1.setVolumes(8)
    print(m1.nodes(), m1.volumes())
    print('_' * 20)
    
    m1 = Mesh(nodes =  5, length = 33)
    print(m1.nodes(), m1.volumes(), m1.length())
    print('_' * 20)
    
    m1 = Mesh(volumes =  5, length = 33)
    print(m1.nodes(), m1.volumes(), m1.length())
    print('_' * 20)
    
    m1 = Mesh(nodes = 5, length = 1)
    print(m1.nodes(), m1.volumes(), m1.length(), m1.delta())
    print('_' * 20)    
    
    m1 = Mesh(volumes = 10, length = 1)
    print(m1.nodes(), m1.volumes(), m1.length(), m1.delta())
    print('_' * 20) 

    m1 = Mesh(volumes = 6, length = 1)
    print(m1.nodes(), m1.volumes(), m1.length(), m1.delta())
    m1.createMesh()
    print('_' * 20) 
    
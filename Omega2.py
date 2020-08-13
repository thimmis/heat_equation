#!/usr/bin/env python3
# @Author: Thomas Turner <thomas>
# @Date:   2020-06-17T13:24:06+02:00
# @Email:  thomas.benjamin.turner@gmail.com
# @Last modified by:   thomas
# @Last modified time: 2020-08-13T15:49:52+02:00
# @Last modified time: 2020-08-13T15:49:52+02:00


import numpy as np
from scipy.linalg import solve
import MatCreator


class Kitchen:
    """
    Desc.


    Attributes:
    -----------

    Methods:
    --------


    """


    def __init__(self, heater, aircon, walls, cols, open = False):
        """


        Params:
        -------


        Returns:
        --------


        """
        self.cols = cols #square domain ~> doubles as rows
        self.dx = 1/self.cols
        self.open = open
        self.tnorth = walls*np.ones(self.cols)
        self.tsouth = self.southwall(heater,aircon,walls)
        self.twest = self.westwall(heater, walls)
        self.A_mat = MatCreator.MCreate(self.cols+2,self.cols+2).A
        self.temp_mat = np.zeros((self.cols+2,self.cols+2))

    def southwall(self,hot, cold,walls):
        """


        Params:
        -------


        Returns:
        --------


        """
        wall = walls*np.ones(self.cols)
        for i in range(self.cols-2,self.cols):
            if self.open == True:
                wall[i] = cold
            else:
                wall[i] = hot
        for i in range(int((self.cols)/2 -3), int((self.cols)/2)):
            wall[i] = hot
        return wall

    def westwall(self, heat, normal):
        """


        Params:
        -------


        Returns:
        --------


        """
        wall = normal*np.ones(self.cols)
        for i in range(int((self.cols)/2),(self.cols)-1):
            wall[i] = heat
        return wall

    def b_vector(self, E):
        N, S, W = self.tnorth, self.tsouth, self.twest
        #first row of temperature values containing tnorth
        b_vec = np.concatenate((np.array([0]), -N, np.array([0])))
        for i in range(0,len(self.twest)):
            temp = np.concatenate((np.array([-W[i]]),np.zeros(self.cols),np.array([-E[i]])))
            b_vec = np.concatenate((b_vec,temp))
        #last row of temperature matrixCreator
        b_vec = np.concatenate( (b_vec,np.array([0]),-S,np.array([0]) ))
        return b_vec



    def get_temp_mat(self,E):
        """


        Params:
        -------


        Returns:
        --------


        """
        b_vec = self.b_vector(-E*self.dx)
        sol_vec = solve(self.A_mat,b_vec)
        self.temp_mat = sol_vec.reshape(self.cols+2,self.cols+2)


    def get_neumann_temps(self):
        """


        Params:
        -------


        Returns:
        --------


        """
        neum_temp1 = self.temp_mat[1:-1,-2]
        return neum_temp1


"""
import matplotlib.pyplot as plt
from matplotlib import *


tester = Kitchen(40,5,10,6)
tester.get_temp_mat(np.ones(6))
fig, ax2 = plt.subplots()
CS = plt.contourf(np.flip(tester.temp_mat,0),10, cmap = plt.cm.inferno)
cbar = fig.colorbar(CS)
plt.show()
"""

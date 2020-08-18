#!/usr/bin/env python3
# @Author: Thomas Turner <thomas>
# @Date:   2020-06-17T13:24:06+02:00
# @Email:  thomas.benjamin.turner@gmail.com
# @Last modified by:   thomas
# @Last modified time: 2020-08-18T11:08:29+02:00
# @Last modified time: 2020-08-18T11:08:29+02:00


import numpy as np
from scipy.linalg import solve
import matrix_creator


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
        self.behaviour_matrix = matrix_creator.FiniteDiffMatrix(self.cols+2,self.cols+2).A
        self.get_temperature_matrix(self.tnorth*self.dx)

    def southwall(self,hot, cold,walls):
        """


        Params:
        -------


        Returns:
        --------


        """
        wall = walls*np.ones(self.cols)
        for i in range(1,int(self.cols/4)-1):
            if self.open == True:
                wall[-i] = cold
            else:
                wall[-i] = hot
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

    def construct_rhs_vector(self, E):
        N, S, W = self.tnorth, self.tsouth, self.twest
        #first row of temperature values containing tnorth
        rhs_vec = np.concatenate((np.array([-N[0]-W[0]]), -N, np.array([-N[-1]-E[0]])))
        for i in range(self.cols):
            temp = np.concatenate((np.array([-W[i]]),np.zeros(self.cols),np.array([E[i]])))
            rhs_vec = np.concatenate((rhs_vec,temp))
        #last row of temperature matrixCreator
        rhs_vec = np.concatenate((rhs_vec,np.array([-S[0]-W[-1]]),-S, np.array([-S[-1]-E[-1]]) ))
        return rhs_vec



    def get_temperature_matrix(self,E):
        """


        Params:
        -------


        Returns:
        --------


        """
        rhs_vector = self.construct_rhs_vector(E)
        temperature_vector = solve(self.behaviour_matrix,rhs_vector)
        self.temperature_matrix = temperature_vector.reshape(self.cols+2,self.cols+2)


    def get_neumann_temps(self):
        """


        Params:
        -------


        Returns:
        --------


        """
        neum_temp1 = -self.temperature_matrix[1:-1,-2]
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

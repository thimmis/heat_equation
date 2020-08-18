#!/usr/bin/env python3
# @Author: Thomas Turner <thomas>
# @Date:   2020-08-03T12:07:54+02:00
# @Email:  thomas.benjamin.turner@gmail.com
# @Last modified by:   thomas
# @Last modified time: 2020-08-18T10:55:57+02:00



from numpy import *
import numpy as np
from scipy.linalg import solve
import matrix_creator



class Entry:
    """
    Desc.


    Attributes:
    -----------

    Methods:
    --------


    """

    def __init__(self, heater, aircon, walls, cols):
        """


        Params:
        -------


        Returns:
        --------


        """
        self.dx = 1/cols
        self.cols = cols
        self.tnorth = self.northwall(aircon, walls, self.cols)
        self.tsouth = walls*np.ones(self.cols)
        self.behaviour_matrix = matrix_creator.FiniteDiffMatrix(self.cols+2,self.cols+2).A
        self.get_temperature_matrix(self.tsouth*self.dx, self.tsouth*self.dx)



    def northwall(self, cold, normal,n):
        """


        Params:
        -------


        Returns:
        --------


        """
        wall = normal*np.ones(n)
        for i in range(int(n/2)-2,int(n/2)+1):
            wall[i] = cold
        return wall


    def construct_rhs_vector(self,W,E):
        """


        Params:
        -------


        Returns:
        --------


        """
        N,S = self.tnorth, self.tsouth
        #first row of temperature values containing tnorth
        rhs_vec = np.concatenate((np.array([-N[0]-W[0]]), -N, np.array([-N[-1]-E[0]])))
        for i in range(0,self.cols):
            temp = np.concatenate((np.array([-W[i]]),np.zeros(self.cols),np.array([-E[i]])))
            rhs_vec = np.concatenate((rhs_vec,temp))
        #last row of temperature matrixCreator
        rhs_vec = np.concatenate( (rhs_vec,np.array([-S[0]-W[-1]]),-S,np.array([-S[-1]-E[-1]]) ))
        return rhs_vec


    def get_temperature_matrix(self, W, E):
        """


        Params:
        -------


        Returns:
        --------


        """
        rhs_vector = self.construct_rhs_vector(W,E)
        temperature_vector = solve(self.behaviour_matrix,rhs_vector)
        self.temperature_matrix = temperature_vector.reshape(self.cols+2,self.cols+2)

    def get_neumann_temps(self):
        """


        Params:
        -------


        Returns:
        --------


        """
        neum_temp1 = self.temperature_matrix[1:-1,-2]
        neum_temp4 = self.temperature_matrix[1:-1,1]
        return neum_temp1, neum_temp4

"""
import matplotlib.pyplot as plt
from matplotlib import *


tester = Entry(10,5,10,6)
tester.get_temp_mat(40*np.ones(6),np.ones(6))
fig, ax2 = plt.subplots()
CS = plt.contourf(np.flip(tester.temp_mat,0),10, cmap = plt.cm.inferno)
cbar = fig.colorbar(CS)
plt.show()
"""

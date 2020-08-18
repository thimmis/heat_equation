#!/usr/bin/env python3
# @Author: Thomas Turner <thomas>
# @Date:   2020-08-03T12:07:54+02:00
# @Email:  thomas.benjamin.turner@gmail.com
# @Last modified by:   thomas
# @Last modified time: 2020-08-18T11:24:14+02:00

from numpy import *
import numpy as np
from scipy.linalg import solve
import matrix_creator

class BathRoom:
    """
    Desc.


    Attributes:
    -----------

    Methods:
    --------


    """

    def __init__(self, aircon, walls, cols):
        """


        Params:
        -------


        Returns:
        --------


        """
        self.dx = 1/cols
        self.cols = cols
        self.rows = 2*cols
        self.tnorth = walls*np.ones(self.cols)
        self.twest = aircon*np.ones(self.rows)
        self.teastt = aircon*np.ones(self.cols)#initial temp guess at boundary
        self.teastb = walls*np.ones(self.cols)
        self.tsouth = aircon*np.ones(self.cols)
        self.behaviour_matrix = matrix_creator.FiniteDiffMatrix(self.rows+2,self.cols+2).A



    def construct_rhs_vector(self, E):
        """


        Params:
        -------


        Returns:
        --------


        """
        N, S, W = self.tnorth, self.tsouth, self.twest
        #first row of temperature values containing tnorth
        rhs_vec = np.concatenate((np.array([-N[0]-W[0]]), -N, np.array([-N[-1]-E[0]])))
        #top block containing interface values
        for i in range(self.cols):
            temp = np.concatenate((np.array([-W[i]]),np.zeros(self.cols),np.array([-E[i]])))
            rhs_vec = np.concatenate((rhs_vec,temp))
        #bottom block with dirichlet conditions
        for i in range(0,self.cols):
            temp = np.concatenate((np.array([-W[i+self.cols]]), np.zeros(self.cols), np.array([-self.teastb[i]])))
            rhs_vec = np.concatenate((rhs_vec, temp))
        #last row of temperature matrixCreator
        rhs_vec = np.concatenate((rhs_vec,np.array([-S[0]-W[-1]]), -S, np.array([-S[-1]-self.teastb[-1]])))
        return rhs_vec

    def temp_gradient_calc(self, E):
        """


        Params:
        -------


        Returns:
        --------


        """
        rhs_vector = self.construct_rhs_vector(E)
        temperature_vector = solve(self.behaviour_matrix, rhs_vector)
        self.temperature_matrix = temperature_vector.reshape(self.rows+2,self.cols+2)
        gradient_3 = (self.temperature_matrix[1:self.cols+1,-1] - E)*self.dx

        return gradient_3



"""
import matplotlib.pyplot as plt
from matplotlib import *

tester = BathRoom(5,10,15)
out = tester.temp_grad(np.ones(15))
fig, ax2 = plt.subplots()
CS = plt.contourf(np.flip(tester.temp_mat,0),10, cmap =  plt.cm.inferno)
cbar = fig.colorbar(CS)
plt.show()
"""

#!/usr/bin/env python3
# @Author: Thomas Turner <thomas>
# @Date:   2020-08-03T12:07:54+02:00
# @Email:  thomas.benjamin.turner@gmail.com
# @Last modified by:   thomas
# @Last modified time: 2020-08-13T15:50:49+02:00

from numpy import *
import numpy as np
from scipy.linalg import solve
import MatCreator

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
        self.teastt = walls*np.ones(self.cols)#initial temp guess at boundary
        self.teastb = walls*np.ones(self.cols)
        self.tsouth = aircon*np.ones(self.cols)
        self.A_mat = MatCreator.MCreate(self.rows+2,self.cols+2).A
        self.temp_mat = np.zeros((self.rows+2,self.cols+2))


    def b_vector(self, E):
        """


        Params:
        -------


        Returns:
        --------


        """
        N, S, W = self.tnorth, self.tsouth, self.twest
        #first row of temperature values containing tnorth
        b_vec = np.concatenate((np.array([0]), -N, np.array([0])))
        #top block containing interface values
        for i in range(self.cols):
            temp = np.concatenate((np.array([-W[i]]),np.zeros(self.cols),np.array([-E[i]])))
            b_vec = np.concatenate((b_vec,temp))
        #bottom block with dirichlet conditions
        for i in range(0,self.cols):
            temp = np.concatenate((np.array([-W[i+self.cols]]),np.zeros(self.cols),np.array([-self.teastb[i]]) ))
            b_vec = np.concatenate((b_vec,temp))
        #last row of temperature matrixCreator
        b_vec = np.concatenate( (b_vec,np.array([0]),-S ,np.array([0])) )
        return b_vec

    def temp_grad(self, E):
        """


        Params:
        -------


        Returns:
        --------


        """
        b_vec = self.b_vector(E)
        sol_vec = solve(self.A_mat, b_vec)
        self.temp_mat = sol_vec.reshape(self.rows+2,self.cols+2)
        gradient_3 = (E - self.temp_mat[1:self.cols+1,-1])*self.dx

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

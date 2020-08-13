#!/usr/bin/env python3
# @Author: Thomas Turner <thomas>
# @Date:   2020-08-03T12:07:54+02:00
# @Email:  thomas.benjamin.turner@gmail.com
# @Last modified by:   thomas
# @Last modified time: 2020-08-13T15:47:46+02:00



from numpy import *
import numpy as np
from scipy.linalg import solve
import MatCreator



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
        self.A_mat = MatCreator.MCreate(self.cols+2,self.cols+2).A
        self.temp_mat = np.zeros((self.cols+2,self.cols+2))




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


    def b_vector(self,W,E):
        """


        Params:
        -------


        Returns:
        --------


        """
        N,S = self.tnorth, self.tsouth
        #first row of temperature values containing tnorth
        b_vec = np.concatenate((np.array([0]), -N, np.array([0])))
        for i in range(0,self.cols):
            temp = np.concatenate((np.array([-W[i]]),np.zeros(self.cols),np.array([-E[i]])))
            b_vec = np.concatenate((b_vec,temp))
        #last row of temperature matrixCreator
        b_vec = np.concatenate( (b_vec,np.array([0]),-S,np.array([0]) ))
        return b_vec


    def get_temp_mat(self, W, E):
        """


        Params:
        -------


        Returns:
        --------


        """
        b_vec = self.b_vector(-W*self.dx,-E*self.dx)
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
        neum_temp4 = self.temp_mat[1:-1,1]
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

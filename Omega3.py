"""
@author Thimmis
"""
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
        self.tnorth = self.northwall(aircon, walls, self.cols-2)
        self.tsouth = walls*np.ones(self.cols-2)
        self.A_mat = MatCreator.MCreate(int(self.cols**2)).A

        self.temp_mat = []
        self.neum_temp1 = []
        self.neum_temp4 = []




    def northwall(self, cold, normal,n):
        """


        Params:
        -------


        Returns:
        --------


        """
        wall = normal*np.ones(n)
        for i in range(n/2,(n/2)+1):
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
        b_vec = np.concatenate((array([0]), (-1)*N, array([0])))
        for i in range(0,len(self.twest)):
            temp = np.concatenate((array([-W[i]]),np.zeros(self.cols-2),array([-E[i]])))
            b_vec = np.concatenate((b_vec,temp))
        #last row of temperature matrixCreator
        b_vec = np.concatenate( (b_vec,array([0]),-1*S,array([0]) ))
        return b_vec


    def get_temp_mat(self, W, E):
        """


        Params:
        -------


        Returns:
        --------


        """
        b_vec = self.b_vector(W,E)
        sol_vec = solve(self.A_mat,b_vec)
        self.temp_mat = sol_vec.reshape(self.cols,self.cols)

    def get_neumann_temps(self):
        """


        Params:
        -------


        Returns:
        --------


        """
        self.neum_temp1 = self.temp_mat[1:-1,-1]
        self.neum_temp4 = self.temp_mat[1:-1,0]

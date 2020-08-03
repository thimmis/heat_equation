"""
@author Thimmis
"""
from numpy import *
import numpy as np
from scipy.linalg import solve
import MatCreator

class Bathroom:
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
        self.tnorth = walls*np.ones(self.cols-2)
        self.twest = aircon*np.ones(self.cols-2)
        self.teast = walls*np.ones(self.cols-2)#initial temp guess at boundary
        self.south = self.tnorth.copy()
        self.A_mat = MatCreator.MCreate(self.cols*self.cols)

        self.gradient_3 = []
        self.temp_mat = []


    def b_vector(self, E):
        """


        Params:
        -------


        Returns:
        --------


        """
        N, S, W = self.tnorth, self.tsouth, self.twest
        #first row of temperature values containing tnorth
        b_vec = np.concatenate((array([0]), (-1)*N, array([0])))
        for i in range(0,len(self.twest)):
            temp = np.concatenate((array([-W[i]]),np.zeros(self.cols-2),array([-E[i]])))
            b_vec = np.concatenate((b_vec,temp))
        #last row of temperature matrixCreator
        b_vec = np.concatenate( (b_vec,array([0]),-1*S,array([0]) ))
        return b_vec

    def temp_grad(self, E):
        """


        Params:
        -------


        Returns:
        --------


        """
        b_vec = self.b_vector(E)
        self.sol_vec = solve(self.A_mat, b_vec)
        self.temp_mat = sol_vec.reshape(self.cols,self.cols)

        temp_3 = slef.temp_mat[1:-1,-1]

        self.gradient_3 = (E - temp_3)*self.dx

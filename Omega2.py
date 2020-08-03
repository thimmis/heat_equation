#!/usr/bin/env python3
# @Author: Thomas Turner <thomas>
# @Date:   2020-06-17T13:24:06+02:00
# @Email:  thomas.benjamin.turner@gmail.com
# @Last modified by:   thomas

# @Last modified time: 2020-08-03T12:13:09+02:00

# @Last modified time: 2020-08-03T12:13:09+02:00


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


    def __init__(self, heater, aircon, walls, sqdim, open = False):
        """


        Params:
        -------


        Returns:
        --------


        """
        self.dx = 1/cols
        self.cols = cols #square domain ~> doubles as rows
        self.open = open
        self.tnorth = walls*np.ones(self.cols-2)
        self.tsouth = self.southwall(heater,aircon,walls)
        self.twest = self.westwall(aircon, walls)
        self.A_mat = MatCreator.MCreate(int(self.cols*self.cols)).A


        self.temp_mat = []
        self.neum_temp1 = []


    def southwall(self,hot, cold):
        """


        Params:
        -------


        Returns:
        --------


        """
        wall = normal*np.ones(self.cols-2)
        for i in range((self.cols-2)-1,self.cols-2):
            if self.open == True:
                wall[i] = cold
            else:
                wall[i] = hot
        for i in range(int((self.cols-2)/2 -2), int((self.cols-2)/2)+1):
            wall[i] = hot
        return wall

    def westwall(self, cold, normal):
        """


        Params:
        -------


        Returns:
        --------


        """
        wall = normal*np.ones(self.cols-2)
        for i in range(int((self.cols-2)/2)+1,(self.cols-2)-2):
            wall[i] = cool
        return wall

    def b_vector(self, E):
        N, S, W = self.tnorth, self.tsouth, self.twest
        #first row of temperature values containing tnorth
        b_vec = np.concatenate((array([0]), (-1)*N, array([0])))
        for i in range(0,len(self.twest)):
            temp = np.concatenate((array([-W[i]]),np.zeros(self.cols-2),array([-E[i]])))
            b_vec = np.concatenate((b_vec,temp))
        #last row of temperature matrixCreator
        b_vec = np.concatenate( (b_vec,array([0]),-1*S,array([0]) ))
        return b_vec



    def get_temp_mat(self,E):
        """


        Params:
        -------


        Returns:
        --------


        """
        b_vec = self.b_vector(E)
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

#!/usr/bin/env python3

# @Author: Thomas Turner <thomas>
# @Date:   2020-06-17T13:23:44+02:00
# @Email:  thomas.benjamin.turner@gmail.com
# @Last modified by:   thomas
<<<<<<< current
# @Last modified time: 2020-08-03T12:09:48+02:00
=======
# @Last modified time: 2020-08-03T12:09:48+02:00
>>>>>>> before discard

import numpy as np
from scipy.linalg import solve
import MatCreator


class LivingRoom:

    """
    Creates approximation of the upper part of the main room of our apartment
    and initializes the temperture at the walls--creating a dynamic temperature vector on the "south"
    wall where there are walls and windows/heaters.

    ~~TODO
    pass in boolean values for open/closed indicating door closed during summer
    and low window/outer wall temperatues in winter.

    Attributes:
    -----------


    Methods:
    --------


    """

    def __init__(self, heater, aircon, walls, cols, open = False):
        self.dx = 1/cols
        self.cols = cols
        self.rows = 2*cols
        self.open = open
        self.tnorth = walls*np.ones(self.cols-2)
        self.teast = walls*np.ones(1)
        self.tsouth = self.southwall(heater,aircon,walls)
        self.twestt = walls*np.ones(int(self.cols/3)) #boundary guess with Entry
        self.twestm = self.teast.copy()
        self.twestb = ((self.tsouth[0] + walls)/2)*np.ones(self.cols-2)#bounday guess with Kitchen
        self.A_mat = MatCreator.MCreate(int((self.cols)*(self.rows))).A


        """
        computes first temperature matrix and gradient vectors to send as Neumann
        conditions to Omegas 2 and 3.
        """
        self.gradient_3 = []
        self.gradient_2 = []
        self.temp_mat = []



    def southwall(self,hot, cold):
        """
        constructs temperature vector for south living room wall

        Params:
        -------



        Returns:
        --------


        """
        wall = normal*np.ones(self.cols-2)
        for i in range(0,int((self.cols-2)/6)):
            if self.open == True:
                wall[i] = cold
            else:
                wall[i] = hot
        for i in range(int((self.cols - 2)/2)-1, int((self.cols - 2)/2 +2)):
            wall[i] = hot
        return wall




    def b_vector(self, WA, WB):
        """
        Creates a 1D array of internal and external temperatures as Dirichlet conditions
        to solve for the internal temperatures.

        Params:
        -------

        Returns:
        --------

        """
        N, S, E, WM = self.tnorth, self.tsouth, self.teast, self.twestm
        #For the first row of temperatures of the matrix A
        b_vec = np.concatenate((array([0]), (-1)*N, array([0])))
        #rows containing interface with entryway
        for i in range(0, len(self.twestt)):
            temp = np.concatenate((array([-WA[i]]),np.zeros(self.cols-2),array([-E[0]])))
            b_vec = np.concatenate((b_vec,temp))
        #rows between interfaces containing wall temps on both sides
        for i in range(0,(self.cols-len(self.twestt))):
            temp = np.concatenate((array([-WM[0]]),np.zeros(self.cols-2),array([-E[0]])))
            b_vec = np.concatenate((b_vec,temp))
        #rows containing interface with kitchen
        for i in range(0,self.cols-2):
            temp = np.concatenate((array([-WB[i]]),np.zeros(self.cols-2),array([-E[0]])))
            b_vec = np.concatenate((b_vec,temp))
        #last row of temperature matrixCreator
        b_vec = np.concatenate( (b_vec,array([0]),-1*S,array([0]) ))
        return b_vec


    def temp_grad(self, WA, WB):
        """
        Takes in temperature values from dirichelt/initial guess temps,
        creates the temperature vector to use in solving the heateq and gets
        the gradients to be passed to OM2 and OM3 as Neumann bounday conditions.

        Params:
        -------


        Returns:
        --------

        """
        b_vec = self.b_vector(WA,WB)
        sol_vec = solve(self.A_mat,b_vec)
        self.temp_mat = sol_vec.reshape(self.rows,self.cols)

        temp_3 = self.temp_mat[:len(self.twestt),0]
        temp_2 = self.temp_mat[self.cols+1:-1,0]

        self.gradient_3 = (temp_3 - WA)*self.dx
        self.gradient_2 = (temp_2 - WB)*self.dx

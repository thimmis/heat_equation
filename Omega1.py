#!/usr/bin/env python3

# @Author: Thomas Turner <thomas>
# @Date:   2020-06-17T13:23:44+02:00
# @Email:  thomas.benjamin.turner@gmail.com
# @Last modified by:   thomas
# @Last modified time: 2020-08-13T15:50:40+02:00
# @Last modified time: 2020-08-13T15:50:40+02:00

import numpy as np
from numpy.linalg import solve
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
    dx : float
        distance between interior gridpoints

    cols : int
        the number of columns of the temperature matrix

    rows : int
        the number of rows of the temperature matrix

    open : bool
        determines if domain is considered to have open or closed door indicating
        that southwall has heat or cooling along boundary in set locations

    tnorth : ndarray
        vector with boundary temperature along the north wall of the domains

    teast : ndarray
        (1,) vector countaining the temperature of the east wall as this has no
        variation

    tsouth : ndarray
        boundary temperature vector for south wall containing windows and door

    twestt : ndarray
        temperture vector containing initial guess for temperature at interface

    twestm : ndarray
        copy of teast since these values do no vary

    twestb : ndarray
        temperature vector containing initial guess at interface with kitchen

    A_mat : ndarray
        tensor matrix describing how the heat equation behaves in rectangular
        domains

    gradient_2 : ndarray
        finite difference gradient to be passed as neumann conditions to kitchen

    gradient_3 : ndarray
        finite difference gradient to be passed as neumann conditions to entryway

    temp_mat : ndarray
        computed temperature distribution used to model domain and calculate the
        gradient vectors passed to the kitchen and entryway.

    Methods:
    --------


    """

    def __init__(self, heater, aircon, walls, cols, open = False):
        self.dx = 1/cols
        self.cols = cols
        self.rows = 2*cols
        self.open = open
        self.tnorth = walls*np.ones(self.cols)
        self.teast = walls*np.ones(1)
        self.tsouth = self.southwall(heater,aircon,walls)
        self.twestt = walls*np.ones(int(self.cols/2)) #boundary guess with Entry
        self.twestm = self.teast.copy()
        self.twestb = walls*np.ones(self.cols)#bounday guess with Kitchen
        self.A_mat = MatCreator.MCreate(self.rows+2,self.cols+2).A
        self.temp_mat = np.zeros((self.rows+2,self.cols+2))


    def southwall(self,hot, cold,walls):
        """
        constructs temperature vector for south living room wall

        Params:
        -------
        hot : float
            hot boundary temperature.

        cold : float
            cool boundary temperature.

        Returns:
        --------
        wall : ndarray
            The vector describing the boundary temperature along the south wall

        """
        wall = walls*np.ones(self.cols)
        for i in range(int((self.cols)/4)):
            if self.open == True:
                wall[i] = cold
            else:
                wall[i] = hot
        for i in range(int(self.cols/2), int(self.cols/2 +3)):
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
        b_vec = np.concatenate((np.array([0]), -N, np.array([0])))
        #rows containing interface with entryway
        for i in range(len(self.twestt)):
            temp = np.concatenate((np.array([-WA[i]]),np.zeros(self.cols),np.array([-E[0]])))
            b_vec = np.concatenate((b_vec,temp))
        #rows between interfaces containing wall temps on both sides
        for i in range((self.cols-len(self.twestt))):
            temp = np.concatenate((np.array([-WM[0]]), np.zeros(self.cols), np.array([-E[0]])))
            b_vec = np.concatenate((b_vec,temp))
        #rows containing interface with kitchen
        for i in range(self.cols):
            temp = np.concatenate((np.array([-WB[i]]),np.zeros(self.cols),np.array([-E[0]])))
            b_vec = np.concatenate((b_vec,temp))
        #last row of temperature matrixCreator
        b_vec = np.concatenate((b_vec,np.array([0]),-S,np.array([0])))
        return b_vec


    def temp_grad(self, WA, WB):
        """
        Takes in temperature values from dirichelt/initial guess temps,
        creates the temperature vector to use in solving the heateq and gets
        the gradients to be passed to OM2 and OM3 as Neumann bounday conditions.
        overwrites the empty lists as ndarrays after a single step of iteration
        has been performed.

        Params:
        -------
        WA : ndarray
            vector containing the boundary temperature at the interface of the
            westwall with the entryway domain

        WB : ndarray
            vector containing the boundary temperature at the interface of the
            westwall with the kitchen domain -- influnced

        Returns:
        --------
        None

        """
        b_vec = self.b_vector(WA,WB)
        sol_vec = solve(self.A_mat,b_vec)
        self.temp_mat = sol_vec.reshape(self.rows+2,self.cols+2)
        """
        gets temperture vectors at west wall interfaces with dirichlet condition
        domains
        """
        gradient_3 = (self.temp_mat[1:len(self.twestt)+1,0] - WA)*self.dx
        gradient_2 = (self.temp_mat[self.cols+1:-1,0] - WB)*self.dx
        return gradient_3, gradient_2

"""
import matplotlib.pyplot as plt
from matplotlib import *


tester = LivingRoom(40,5,10,20)
tester.temp_grad(tester.twestt,tester.twestb)
fig, ax2 = plt.subplots()
CS = plt.contourf(np.flip(tester.temp_mat,0),10, cmap = plt.cm.inferno)
cbar = fig.colorbar(CS)
plt.show()
"""

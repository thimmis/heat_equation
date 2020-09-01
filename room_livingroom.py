#!/usr/bin/env python3

# @Author: Thomas Turner <thomas>
# @Date:   2020-06-17T13:23:44+02:00
# @Email:  thomas.benjamin.turner@gmail.com
# @Last modified by:   thomas
# @Last modified time: 2020-09-01T13:46:13+02:00
# @Last modified time: 2020-09-01T13:46:13+02:00

import numpy as np
from numpy.linalg import solve
import matrix_creator


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

    temp_mat : ndarray
        computed temperature distribution used to model domain and calculate the
        gradient vectors passed to the kitchen and entryway.

    Methods:
    --------
    sw_temp_open_close(self,hot, cold, walls)
        constructs a south wall temperature vector scaled for the number of horizontal
        gridpoints used. Depending on open or closed will prescribe hot or cold
        temperatures where the door is.

    construct_rhs_vector(self, WA, WB)
        Constructs the vector used for solving the linear equation Ax = b.

    temp_gradient_calc(self, WA, WB)
        Takes the upper and lower temperature information to produce a new
        temperature distribution for each step. Calculates the gradients at
        each step to be passed as neumann conditions to the kitchen and
        entryway domains.


    """

    def __init__(self, heater, aircon, walls, cols, open = False):
        self.dx = 1/cols
        self.cols = cols
        self.rows = 2*cols
        self.open = open
        self.tnorth = walls*np.ones(self.cols)
        self.teast = walls*np.ones(1)
        self.tsouth = self.sw_temp_open_close(heater,aircon,walls)
        self.twestt = walls*np.ones(int(self.cols/2)) #boundary guess with Entry
        self.twestm = self.teast.copy()
        self.twestb = walls*np.ones(self.cols)#bounday guess with Kitchen
        self.behaviour_matrix = matrix_creator.FiniteDiffMatrix(self.rows+2,self.cols+2,'').A



    def sw_temp_open_close(self,hot, cold, walls):
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
        for i in range(int((self.cols)/4 -int(self.cols/10))):
            if self.open == True:
                wall[i] = cold
            else:
                wall[i] = hot
        for i in range(int(self.cols/2), int(self.cols/2)+int(self.cols/10)+1):
            wall[i] = hot
        return wall




    def construct_rhs_vector(self, WA, WB):
        """
        Creates a 1D array of internal and external temperatures as Dirichlet conditions
        to solve for the internal temperatures.

        Params:
        -------
        WA : ndarray
            vector containing the boundary temperature at the interface of the
            westwall with the entryway domain.

        WB : ndarray
            vector containing the boundary temperature at the interface of the
            westwall with the kitchen domain.

        Returns:
        --------
        rhs_vec : ndarray
            The vector of length (cols*cols) containing boundary temperatures and
            zeros for unknown values.

        """
        N, S, E, WM = self.tnorth, self.tsouth, self.teast, self.twestm
        #For the first row of temperatures of the matrix A
        rhs_vec = np.concatenate((np.array([-N[0]-WA[0]]), -N, np.array([-N[-1]-E[0]])))
        #rows containing interface with entryway
        for i in range(len(self.twestt)):
            temp = np.concatenate((np.array([-WA[i]]),np.zeros(self.cols),np.array([-E[0]])))
            rhs_vec = np.concatenate((rhs_vec,temp))
        #rows between interfaces containing wall temps on both sides
        for i in range((self.cols-len(self.twestt))):
            temp = np.concatenate((np.array([-WM[0]]), np.zeros(self.cols), np.array([-E[0]])))
            rhs_vec = np.concatenate((rhs_vec,temp))
        #rows containing interface with kitchen
        for i in range(self.cols):
            temp = np.concatenate((np.array([-WB[i]]),np.zeros(self.cols),np.array([-E[0]])))
            rhs_vec = np.concatenate((rhs_vec,temp))
        #last row of temperature matrixCreator
        rhs_vec = np.concatenate((rhs_vec,np.array([-S[0]-WB[-1]]),-S,np.array([-S[-1]-E[0]])))
        return rhs_vec


    def temp_gradient_calc(self, WA, WB):
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
            westwall with the entryway domain.

        WB : ndarray
            vector containing the boundary temperature at the interface of the
            westwall with the kitchen domain.

        Returns:
        --------
        gradient_2 : ndarray
            Gradient vector to be passed to the kitchen as Neumann BCs.

        gradient_3 : ndarray
            Gradient vector to be passed to the entryway as Neumann BCs.

        """
        rhs_vector = self.construct_rhs_vector(WA,WB)
        temperature_vector = solve(self.behaviour_matrix,rhs_vector)
        self.temperature_matrix = temperature_vector.reshape(self.rows+2,self.cols+2)
        #calculate the gradients at the upper and lower westwall interfaces
        gradient_3 = (WA - self.temperature_matrix[1:len(self.twestt)+1,0])*self.dx
        gradient_2 = (WB - self.temperature_matrix[self.cols+1:-1,0])*self.dx
        return gradient_3, gradient_2

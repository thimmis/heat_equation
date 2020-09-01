#!/usr/bin/env python3
# @Author: Thomas Turner <thomas>
# @Date:   2020-06-17T13:24:06+02:00
# @Email:  thomas.benjamin.turner@gmail.com
# @Last modified by:   thomas
# @Last modified time: 2020-09-01T13:51:33+02:00
# @Last modified time: 2020-09-01T13:51:33+02:00


import numpy as np
from scipy.linalg import solve
import matrix_creator


class Kitchen:
    """
    Sets up the 2D Heat Equation in our old kitchen. Applies Nuemann boundary
    conditions at the east interface with the living room--flux along the east
    boundary. All others are Dirichlet conditions. It uses boolean variables to
    create dynamic temperature vectors for the west and south walls.


    Attributes:
    -----------
    cols : int
        The number of columns of interior gridpoints to be solved for.

    dx : int
        One over the number of columns of gridpoints giving a delta_x.

    open : bool
        Tells the southwall method if the temperature at the door should be hot
        or cold i.e. if the door is open or closed

    oven : bool
        Tell the westwall method if the temperature at the oven's location should
        be hot or cold i.e. the oven/stove is on or off.

    tnorth : ndarray
        A vector of constant wall temperatures.

    tsouth : ndarray
        The vector of temperatures for the south wall containing a window and
        the door to the patio.

    twest : ndarray
        The temperature vector for the wall containing the oven and the vent
        over it.

    behaviour_matrix : ndarray
        The finite difference matrix that solves for the unknown temperature
        points for the kitchen domain. Adjusted to consider the Neumann BC at
        the east wall--the interface with the living room.



    Methods:
    --------
    sw_temp_open_close(self,hot, cold,walls)

    westwall(self, hot, normal, cold)

    construct_rhs_vector(self, E)

    get_temperature_matrix(self,E)

    get_neumann_temps(self)

    """


    def __init__(self, heater, aircon, walls, cols, open = False, on_off = False):
        """


        Params:
        -------


        Returns:
        --------


        """
        self.cols = cols #square domain ~> doubles as rows
        self.dx = 1/self.cols
        self.open = open
        self.oven = on_off
        self.tnorth = walls*np.ones(self.cols)
        self.tsouth = self.sw_temp_open_close(heater,aircon,walls)
        self.twest = self.westwall(heater, walls, aircon)
        self.behaviour_matrix = matrix_creator.FiniteDiffMatrix(self.cols+2,self.cols+2, 'r').A
        self.get_temperature_matrix(self.tnorth*self.dx)

    def sw_temp_open_close(self,hot, cold,normal):
        """
        Creates dynamic south wall temperature vector scaling to the number of
        horizontal gridpoints.

        Params:
        -------
        hot : int
            Temperature prescribed to boundary heat sources.

        cold : int
            Temperature prescribed to boundary cold sources.

        normal : int
            Temperature prescribed as average boundary wall temperatures.


        Returns:
        --------
        wall : ndarray
            Dynamic temperature vector with hot and cold spots scaled for the
            number of horizontal gridpoints used.


        """
        wall = normal*np.ones(self.cols)
        for i in range(1,int(self.cols/4)-int(self.cols/10)):
            if self.open == True:
                wall[-i] = cold
            else:
                wall[-i] = hot
        for i in range(int(self.cols/2) - int(self.cols/10)), int((self.cols)/2)):
            wall[i] = hot
        return wall

    def westwall(self, hot, normal, cold):
        """
        Creates dynamic west wall temperature vector scaling to the number of
        horizontal gridpoints.


        Params:
        -------
        hot : int
            Temperature prescribed to boundary heat sources.

        cold : int
            Temperature prescribed to boundary cold sources.

        normal : int
            Temperature prescribed as average boundary wall temperatures.

        Returns:
        --------
        wall : ndarray
            Dynamic temperature vector with hot and cold spots scaled for the
            number of horizontal gridpoints used.

        """
        wall = normal*np.ones(self.cols)
        for i in range(int((self.cols)/2),(self.cols)-int(self.cols/10)):
            if self.oven == True:
                wall[i] = hot
            else:
                wall[i] = cold
        return wall

    def construct_rhs_vector(self, E):
        """
        Creates the vector needed to solve the linear equation approximating the
        2D heat equation.


        Params:
        -------
        E : ndarray
            Vector containing temperatures to compute heat distrbution updated
            at each step of iteration.

        Returns:
        --------
        rhs_vec : ndarray
            The vector of length (cols*rows) containing boundary temperatures and
            zeros for unknown values.

        """
        N, S, W = self.tnorth, self.tsouth, self.twest
        #first row of temperature values containing tnorth
        rhs_vec = np.concatenate((np.array([-N[0]-W[0]]), -N, np.array([-N[-1]-E[0]])))
        for i in range(self.cols):
            temp = np.concatenate((np.array([-W[i]]),np.zeros(self.cols),np.array([-E[i]])))
            rhs_vec = np.concatenate((rhs_vec,temp))
        #last row of temperature matrixCreator
        rhs_vec = np.concatenate((rhs_vec,np.array([-S[0]-W[-1]]),-S, np.array([-S[-1]-E[-1]]) ))
        return rhs_vec



    def get_temperature_matrix(self,E):
        """
        Contructs the right hand side vector for solving the problem on the current
        iteration. Comptutes the temperature distribution and reshapes the vector
        into a matrix.

        Params:
        -------
        E : ndarray
            Vector containing temperatures to compute heat distrbution updated
            at each step of iteration for the east interface.

        Returns:
        --------
        None


        """
        rhs_vector = self.construct_rhs_vector(E)
        temperature_vector = solve(self.behaviour_matrix,rhs_vector)
        self.temperature_matrix = temperature_vector.reshape(self.cols+2,self.cols+2)


    def get_neumann_temps(self):
        """
        Pulls the right most column to be sent to the living as new dirichelt
        BCs.

        Params:
        -------
        None

        Returns:
        --------
        neum_temp1 : ndarray
           Temperature vector to be sent to the living room to be used in the
           next iteration step.


        """
        neum_temp1 = self.temperature_matrix[1:-1,-1]
        return neum_temp1

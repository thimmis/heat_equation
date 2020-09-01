#!/usr/bin/env python3
# @Author: Thomas Turner <thomas>
# @Date:   2020-08-03T12:07:54+02:00
# @Email:  thomas.benjamin.turner@gmail.com
# @Last modified by:   thomas
# @Last modified time: 2020-09-01T13:08:51+02:00



from numpy import *
import numpy as np
from scipy.linalg import solve
import matrix_creator



class Entry:
    """
    Sets up the 2D heat equation in the entryway of LGH 1207. This domain features
    dirichelt conditions on the north and south walls and neumann conditions at
    the east and west interfaces with the bathroom and living rooms.


    Attributes:
    -----------

    dx : float
        Gives delta_x for the number of horizontal gridpoints.

    cols : int
        Number of interior horizontal gridpoints.

    tnorth : ndarray
        temperture vector for the northwall of the entryway featuring variation
        where the door is as a cold point.

    tsouth : ndarray
        temperature vector for southwall of the entryway.

    behaviour_matrix : ndarray
        finite difference matrix using pure dirichelt BCs.

    temperature_matrix : ndarray
        matrix holding the computed temperature that is updated at each step
        of the iteration.
    Methods:
    --------

    northwall(self, cold, normal,cols)
        A function to scale the scale the door to the number of horizontal gridpoints
        used.

    construct_rhs_vector(self,W,E)
        creates the right hand side vector using the two gradients passed into.

    get_temperature_matrix(self, W, E)
        takes in the gradient vectors at the west and east interfaces and constructs
        the right hand side vector. Uses the rhs_vector to solve for the unknowns
        and reshapes into the temperature matrix. Separated from getting temperature
        vectors to allow for relaxation step.

    get_neumann_temps(self)
        pulls the left and right most columns to be sent back to livingroom and
        bathroom as dirichelt BCs.


    """

    def __init__(self, heater, aircon, walls, cols):
        """
        Sets up the 2D linear problem for the entryway domain.

        Params:
        -------
        aricon : int
            temperature to be prescribed to cold sources in boundaries.

        walls : int
            temperature to be presecribed as average wall temperatures in
            boundaries.

        cols : int
            Number of interior horizontal gridpoints.



        Returns:
        --------


        """
        self.dx = 1/cols
        self.cols = cols
        self.tnorth = self.northwall(aircon, walls, self.cols)
        self.tsouth = walls*np.ones(self.cols)
        self.behaviour_matrix = matrix_creator.FiniteDiffMatrix(self.cols+2, self.cols+2,'lr').A
        self.get_temperature_matrix(self.tsouth*self.dx, self.tsouth*self.dx)



    def northwall(self, cold, normal,cols):
        """
        Creates dynamic northwall temperture vector with door featured as cold
        source. Scales door size and location based on number of horizontal
        gridpoints used.

        Params:
        -------
        cold : int
            Temperature prescribed to cold source.

        normal : int
            temperture presecribed to average wall temperatures.

        cols :  int
            number of horizontal gridpoints.

        Returns:
        --------
        wall : ndarray
            Temperature vector for northwall with non-constant temperatures as
            dirichelt BCs.


        """
        wall = normal*np.ones(cols)
        for i in range(int(cols/2)-(int(cols/5)-2),int(cols/2)+(int(cols/5)+1)):
            wall[i] = cold
        return wall


    def construct_rhs_vector(self,W,E):
        """
        Creates the vector needed to solve the linear equation approximating the
        2D heat equation.

        Params:
        -------
        W : ndarray
            Vector containing temperatures to compute heat distrbution updated
            at each step of iteration for the west interface.


        E : ndarray
            Vector containing temperatures to compute heat distrbution updated
            at each step of iteration for the east interface.


        Returns:
        --------
        rhs_vec : ndarray
            The vector of length (cols*cols) containing boundary temperatures and
            zeros for unknown values.


        """
        N,S = self.tnorth, self.tsouth
        #first row of temperature values containing tnorth
        rhs_vec = np.concatenate((np.array([-N[0]-W[0]]), -N, np.array([-N[-1]-E[0]])))
        for i in range(0,self.cols):
            temp = np.concatenate((np.array([-W[i]]),np.zeros(self.cols),np.array([-E[i]])))
            rhs_vec = np.concatenate((rhs_vec,temp))
        #last row of temperature matrixCreator
        rhs_vec = np.concatenate( (rhs_vec,np.array([-S[0]-W[-1]]),-S,np.array([-S[-1]-E[-1]]) ))
        return rhs_vec


    def get_temperature_matrix(self, W, E):
        """
        Contructs the right hand side vector for solving the problem on the current
        iteration. Comptutes the temperature distribution and reshapes the vector
        into a matrix.

        Params:
        -------
        W : ndarray
            Vector containing temperatures to compute heat distrbution updated
            at each step of iteration for the west interface.


        E : ndarray
            Vector containing temperatures to compute heat distrbution updated
            at each step of iteration for the east interface.

        Returns:
        --------
        None


        """
        rhs_vector = self.construct_rhs_vector(W,E)
        temperature_vector = solve(self.behaviour_matrix,rhs_vector)
        self.temperature_matrix = temperature_vector.reshape(self.cols+2,self.cols+2)

    def get_neumann_temps(self):
        """
        pulls the left and right most columns to be sent back to livingroom and
        bathroom as dirichelt BCs.


        Params:
        -------
        None


        Returns:
        --------
         neum_temp1 : ndarray
            Temperature vector to be sent to the living room to be used in the
            next iteration step.

         neum_temp4 : ndarray
            Temperature vector to be sent to the bathroom to be used in the
            next iteration step.

        """
        neum_temp1 = self.temperature_matrix[1:-1,-1]
        neum_temp4 = self.temperature_matrix[1:-1,0]
        return neum_temp1, neum_temp4

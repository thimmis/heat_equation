#!/usr/bin/env python3
# @Author: Thomas Turner <thomas>
# @Date:   2020-08-03T12:07:54+02:00
# @Email:  thomas.benjamin.turner@gmail.com
# @Last modified by:   thomas
# @Last modified time: 2020-09-01T11:50:47+02:00

from numpy import *
import numpy as np
from scipy.linalg import solve
import matrix_creator

class BathRoom:
    """
    A class that sets up the 2D heat equation in the bathroom. This domain uses
    pure dirichelt BCs.


    Attributes:
    -----------
    dx : float
        Gives delta_x for the number of horizontal gridpoints.

    cols : int
        Number of interior horizontal gridpoints.

    rows : int
        The number of interior vertical gridpoints.

    tnorth : ndarray
        temperture vector for the northwall of the bathroom.

    twest : ndarray
        temperture vector for the westwall of the bathroom.

    teastt : ndarray
        Initial temperture guess at interface with entryway.

    teastb : ndarray
        temperture vector for dirichelt conditions on lower eastwall.

    tsouth : ndarray
        temperature vector for southwall of the bathroom.

    behaviour_matrix : ndarray
        finite difference matrix using pure dirichelt BCs.

    temperature_matrix : ndarray
        matrix holding the computed temperature that is updated at each step
        of the iteration.

    Methods:
    --------
    construct_rhs_vector(self,E)
        Creates the right hand side vector b for solving the linear equation
        Ax = b.

    temp_gradient_calc(self, E)
        At each step of the iteration produces the rhs vector to calculate the
        temperature distribution. One this has been it it reshapes the solution
        vector back into a matrix and determines an approximation of the gradient
        at the inteface to pass to the enteryway.


    """

    def __init__(self, aircon, walls, cols):
        """
        Set up the 2D heat equation problem for the bathroom which shares an
        interface with the entryway. This room uses pure dirichelt conditions.


        Params:
        -------
        aricon : int
            temperature to be prescribed to cold sources in boundaries.

        walls : int
            temperature to be presecribed as average wall temperatures in
            boundaries.

        cols : int
            Number of interior horizontal gridpoints.

        """
        self.dx = 1/cols
        self.cols = cols
        self.rows = 2*cols
        self.tnorth = walls*np.ones(self.cols)
        self.twest = aircon*np.ones(self.rows)
        self.teastt = walls*np.ones(self.cols)#initial temp guess at boundary
        self.teastb = walls*np.ones(self.cols)
        self.tsouth = aircon*np.ones(self.cols)
        self.behaviour_matrix = matrix_creator.FiniteDiffMatrix(self.rows+2,self.cols+2,'').A



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
        #top block containing interface values
        for i in range(0,self.cols):
            temp = np.concatenate((np.array([-W[i]]),np.zeros(self.cols),np.array([-E[i]])))
            rhs_vec = np.concatenate((rhs_vec,temp))
        #bottom block with dirichlet conditions
        for i in range(0,self.cols):
            temp = np.concatenate((np.array([-W[i+self.cols]]), np.zeros(self.cols), np.array([-self.teastb[i]])))
            rhs_vec = np.concatenate((rhs_vec, temp))
        #last row of temperature matrixCreator
        rhs_vec = np.concatenate((rhs_vec,np.array([-S[0]-W[-1]]), -S, np.array([-S[-1]-self.teastb[-1]])))
        return rhs_vec

    def temp_gradient_calc(self, E):
        """
        Builds the rhs_vec and computes the temperature distribution. Reshapes
        the vector and computes the gradient to be passed to the entryway.



        Params:
        -------
        E : ndarray
            Vector containing temperatures to compute heat distrbution updated
            at each step of iteration.


        Returns:
        --------
        gradient_3 : ndarray
            Vector describing flux along the east bounday to be passed to process
            3--entryway.


        """
        rhs_vector = self.construct_rhs_vector(E)
        temperature_vector = solve(self.behaviour_matrix, rhs_vector)
        self.temperature_matrix = temperature_vector.reshape(self.rows+2,self.cols+2)
        gradient_3 = (self.temperature_matrix[1:self.cols+1,-1] - E)*self.dx
        return gradient_3

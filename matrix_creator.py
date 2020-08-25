#!/usr/bin/env python3
# @Author: Thomas Turner <thomas>
# @Date:   2020-06-17T13:27:38+02:00
# @Email:  thomas.benjamin.turner@gmail.com
# @Last modified by:   thomas
# @Last modified time: 2020-08-25T15:23:40+02:00


from scipy import *
import numpy as np


class FiniteDiffMatrix:
    """

    Attributes:
    -----------

    Methods:
    --------


    """

    def __init__(self, rows, cols, condition = False):
        """

        Params:
        -------


        """
        self.cols = cols
        self.rows = rows
        self.condition = condition
        self.squaredim = rows*cols
        self.stamp = self.create_kron_stamp()
        self.A = self.create_behaviour_matrix()

    def create_behaviour_matrix(self):
        """

        Params:
        -------

        Returns:
        --------

        """
        #produces main finite diff matrix using kronecker product
        A = (np.kron(np.eye(self.rows),self.stamp)
            #concatenates the lower and upper diagonals 
            + np.diag([1]*(self.squaredim-self.cols),k=self.cols)
            + np.diag([1]*(self.squaredim-self.cols),k=-self.cols))

        return A


    def create_kron_stamp(self):
        """

        """
        #first outline the generic stamp for Dirichlet boundary conditions
        stamp = (np.diagflat([-4]*self.cols)
                + np.diag([1]*(self.cols - 1), k=1)
                + np.diag([1]*(self.cols - 1),k=-1))

        #Adjusts kronecker product stamp to include left neumann BCs
        if 'l' in self.condition:
            stamp[0,1] = 2

        #Adjusts kronecker product stamp to include right neumann BCs
        if 'r' in self.condition:
            stamp[-1,-2] = 2
        else:
            None
        return stamp

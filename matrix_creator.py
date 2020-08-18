#!/usr/bin/env python3
# @Author: Thomas Turner <thomas>
# @Date:   2020-06-17T13:27:38+02:00
# @Email:  thomas.benjamin.turner@gmail.com
# @Last modified by:   thomas
# @Last modified time: 2020-08-17T17:17:43+02:00


from scipy import *
import numpy as np


class FiniteDiffMatrix:
    """

    Attributes:
    -----------

    Methods:
    --------


    """

    def __init__(self, rows, cols):
        """

        Params:
        -------


        """
        self.cols = cols
        self.squaredim = rows*cols
        self.A = self.create_square_matrix()

    def create_square_matrix(self):
        """

        Params:
        -------

        Returns:
        --------

        """

        A = np.diagflat([-4]*self.squaredim) + np.diag([1]*(self.squaredim - 1), k=1) +np.diag([1]*(self.squaredim - 1), k=-1)  + np.diag([1]*(self.squaredim-self.cols),k=self.cols) + np.diag([1]*(self.squaredim-self.cols),k=-self.cols)
        for i in range(1, self.squaredim):
            if i % self.cols == 0:
                A[i-1,i] = A[i,i-1] = 0
        A[0,0] = A[-1,-1] = -3
        return A

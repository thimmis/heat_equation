#!/usr/bin/env python3
# @Author: Thomas Turner <thomas>
# @Date:   2020-06-17T13:27:38+02:00
# @Email:  thomas.benjamin.turner@gmail.com
# @Last modified by:   thomas
# @Last modified time: 2020-08-12T16:46:45+02:00


from scipy import *
import numpy as np


class MCreate:
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
        self.A = self.get_matrix()

    def get_matrix(self):
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
        return A

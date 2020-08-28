#!/usr/bin/env python3
# @Author: Thomas Turner <thomas>
# @Date:   2020-06-17T13:27:38+02:00
# @Email:  thomas.benjamin.turner@gmail.com
# @Last modified by:   thomas
# @Last modified time: 2020-08-28T16:50:19+02:00


from scipy import *
import numpy as np


class FiniteDiffMatrix:

    """
    Sets up a finite difference matrix for solving the 2D heat equation i.e.
    the Laplace/Poisson equation. Creates a pattern stamp give the number of columns
    in the matrix to create a stamp with the appropriate boundary conditions for
    the given room. Using the stamp it performs a kronecker product on an identity
    matrix for the number of rows in the matrix. Then it adds the sub and super
    off diagonal elements for a matrix of square size rows*columns.

    Attributes:
    -----------
    rows : int
        Number of interior gridpoints plus two to take into consideration boundary temperatues

    cols: int
        Number of interior gridpoints plus two to take into consideration boundary temperatues

    condition : string
        An empty string (' ') provides true dirichlet BCs, permutations of 'l' and
        'r' determine if the left or right interfaces have neumann BCs.

    squaredim : int
        The final square dimension of the finite differnce matrix, used to add
        the sub and super diagonal rows of [1]

    stamp : ndarray
        Array with shape (cols,cols) describing how solution is calculated and
        accounts for the stated boundary conditions.

    A : ndarray
        The matrix needed to solve the linear equation A*x = b where * denotes
        classical matrix multiplication.

    Methods:
    --------
    create_solution_matrix(self)
        Creates the solution matrix A needed to solve the linear problem describing
        the finite difference method for solving the 2D heat equation.

    create_kron_stamp(self)
        Determines the stamp to be used in the Kronecker(Tensor) product to setup
        the three main diagonal elements defining the problem.

    """

    def __init__(self, rows, cols, condition = None):
        """

        Params:
        -------
        rows : int
            Number of interior gridpoints plus two to take into consideration
            boundary temperatues

        cols: int
            Number of interior gridpoints plus two to take into consideration
            boundary temperatues

        condition : string
            An empty string (' ') provides true dirichlet BCs, permutations of
            'l' and 'r' determine if the left or right interfaces have neumann
            BCs.

        squaredim : int
            The final square dimension of the finite differnce matrix, used to
            add the sub and super diagonal rows of [1]

        stamp : ndarray
            Array with shape (cols,cols) describing how solution is calculated and
            accounts for the stated boundary conditions.

        A : ndarray
            The matrix needed to solve the linear equation A*x = b where * denotes
            classical matrix multiplication.

        """
        self.cols = cols
        self.rows = rows
        self.squaredim = rows*cols
        self.A = self.create_solution_matrix(self.create_kron_stamp(condition))

    def create_solution_matrix(self, stamp):
        """
        Uses the stamp determined by the dimension of the problem domain to
        produce the solution matrix with the proper dimensions for the number of
        unknowns to solve for.

        Params:
        -------
        stamp : ndarray
            The pattern for a solving for a single row of unknown temperture
            values.

        Returns:
        --------
        A : ndarray
            The matrix needed to solve the linear equation A*x = b where * denotes
            classical matrix multiplication.
        """
        #produces main finite diffence matrix using kronecker product
        A = (np.kron(np.eye(self.rows),stamp)
            #concatenates the sub and super off diagonals
            + np.diag([1]*(self.squaredim-self.cols),k=self.cols)
            + np.diag([1]*(self.squaredim-self.cols),k=-self.cols))

        return A


    def create_kron_stamp(self, condition):
        """
        Creates a a scaled pattern for how a single row of unknown temperature
        values at interior dridpoints are calculated. If condition != None then
        if the string passed into condition looks for 'l' and 'r' to determine
        where the boundary conditions are homogenous Nuemann.


        Params:
        -------
        condition : string
            Used to denote the presence of Neumann BCs to the right 'r' and left
            'l'.

        Returns:
        --------
        stamp: ndarray
            The matrix with the pattern to be used in the Kronecker product.

        """
        #first outline the generic stamp for Dirichlet boundary conditions
        stamp = (np.diagflat([-4]*self.cols)
                + np.diag([1]*(self.cols - 1), k=1)
                + np.diag([1]*(self.cols - 1),k=-1))

        #Adjusts kronecker product stamp to include left neumann BCs
        if 'l' in condition:
            stamp[0,1] = 2

        #Adjusts kronecker product stamp to include right neumann BCs
        if 'r' in condition:
            stamp[-1,-2] = 2
        else:
            stamp = stamp
        return stamp

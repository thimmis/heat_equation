#!/usr/bin/env python3
# @Author: Thomas Turner <thomas>
# @Date:   2020-08-03T12:07:54+02:00
# @Email:  thomas.benjamin.turner@gmail.com
# @Last modified by:   thomas
# @Last modified time: 2020-08-03T12:10:07+02:00



from mpi4py import MPI
from scipy import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import *
import Omega1, Omega2, Omega3, Omega4#, Plotter

class Solver:


    def DirNeuIter(self, heater, aircon, walls, cols, iterations, open):
        """Solves the 2D-Heat Equation iteratively perscribing Dirichlet and
        Neumann boundary conditions to the different domains.

        Params:
        -------

        Returns:
        -------
        om1, om2, om3, om4 : ndarray, ndarray, ndarray, ndarray
            matrices with computed heat distribution
        """

        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        nprocessors = comm.Get_size()
        iters = iterations


        omega1 = Omega1.LivingRoom(heater, aircon, walls, cols, open)
        omega2 = Omega2.Kitchen

        i = 0

        return om1, om2, om3, om4


class Problem(Solver):

    def __init__(self, heater, aircon, walls):
        """
        Params:
        -------
        heater : float
            the value prescribed to walls with a Heater
        aircon : float
            the value prescribed to walls with a Window
        wall : float
            the value prescribed to walls without a window or heater

        """
        self.heater = heater
        self.aircon = aircon
        self.wall = walls


    def __call__(self, cols, iters, open = False):
        """
        Params:
        -------

        Returns:
        -------
        None

        """
        self.OM1, self.OM2, self.OM3, self.OM4 = DirNeuIter(self, heater, aircon, walls, cols, iterations, open)

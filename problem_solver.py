#!/usr/bin/env python3
# @Author: Thomas Turner <thomas>
# @Date:   2020-08-03T12:07:54+02:00
# @Email:  thomas.benjamin.turner@gmail.com
# @Last modified by:   thomas
# @Last modified time: 2020-08-28T17:01:11+02:00



from mpi4py import MPI
import numpy as np
import room_kitchen, room_bathroom, room_livingroom, room_entryway, plot_domain

class Solver:
    """
    Performs parallel computation of the 2D-Heat Equation on five processes.
    Implements blocking communication from MPI4py to generate solutions by
    iteratively prescribing Dirichlet and Neumann boundary conditions.

    Attributes:
    ----------
    None
``
    Methods:
    -------
    dirichelt_neumann_iteration(self, heater, aircon, cols, iters, open, on_off)
        Implements Dirichlet/Neumann iteration to solve the 2D-Heat Equation

    """
    def dirichelt_neumann_iteration(self, heater, aircon, walls, cols, iters, open, on_off):
        """
        Solves the 2D-Heat Equation iteratively perscribing Dirichlet and
        Neumann boundary conditions to the different domains.

        Uses five processes (0-4) to produce a convergent solution.



        Params:
        -------
        heater : int
            The temperature presecribed to sources of heat in Dirichlet BCs

        aircon : int
            The temperature presecribed as cold source/heatsink in Dirichlet BCs

        walls : int
            The temperature presecribed as neutral wall temperature in Dir BCs.

        cols : int
            The number of columns of interior gridpoints to be solved for.

        open : bool
            Boolean value passed into the problem determin whether or not the patio
            door is open(True) or closed(False) as default.

        on_off : bool
            Boolean value indicating whether the oven/stove is on or off thereby
            acting as a source of heat or a source of cold due to constant vaccum
            being applied at that location.

        Returns:
        -------
        om1, om2, om3, om4 : ndarray, ndarray, ndarray, ndarray
            Matrices containing the computed heat values.
        """

        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        nprocessors = comm.Get_size()
        iterations = iters


        livingroom = room_livingroom.LivingRoom(heater, aircon, walls, cols, open)
        kitchen = room_kitchen.Kitchen(heater, aircon, walls, cols, open, on_off)
        entryway = room_entryway.Entry(heater, aircon, walls, int(cols/2))
        bathroom = room_bathroom.BathRoom(aircon, walls, int(cols/2))

        i = 0

        while i != iterations:
            if i == 0:
                if rank == 1:
                    g_13, g_12 = livingroom.temp_gradient_calc(livingroom.twestt,livingroom.twestb)
                    comm.send(g_12, dest=2, tag=12)
                    comm.send(g_13, dest=3, tag=13)

                if rank == 4:
                    g_43 = bathroom.temp_gradient_calc(bathroom.teastt)
                    comm.send(g_43, dest=3, tag=43)

                if rank == 2:
                    temp2_old = kitchen.temperature_matrix #needed to relax temperature calculations
                    data2_in = comm.recv(source=1, tag=12)
                    kitchen.get_temperature_matrix(data2_in)
                    #relaxation step necessary for producing convergent solution.
                    kitchen.temperature_matrix = (0.8)*kitchen.temperature_matrix + (0.2)*temp2_old
                    data2_out = kitchen.get_neumann_temps()
                    comm.send(data2_out, dest=1, tag=21)

                if rank == 3:
                    temp3_old = entryway.temperature_matrix #needed to relax temperature calculations
                    data31_in = comm.recv(source=1, tag=13)
                    data34_in = comm.recv(source=4, tag=43)
                    entryway.get_temperature_matrix(data34_in, data31_in)
                    #relaxation step necessary for producing convergent solution.
                    entryway.temperature_matrix = (0.8)*entryway.temperature_matrix + (0.2)*temp3_old
                    data31_out, data34_out = entryway.get_neumann_temps()
                    comm.send(data31_out, dest=1, tag=31)
                    comm.send(data34_out, dest=4, tag=34)

                i +=1
            else:
                if rank == 1:
                    btemp2 = comm.recv(source=2, tag=21)
                    btemp3 = comm.recv(source=3, tag=31)
                    g_13, g_12 = livingroom.temp_gradient_calc(btemp3, btemp2)
                    comm.send(g_12, dest=2, tag=12)
                    comm.send(g_13, dest=3, tag=13)


                if rank == 4:
                    btemp4 = comm.recv(source=3, tag=34)
                    g_43 = bathroom.temp_gradient_calc(btemp4)
                    comm.send(g_43, dest=3, tag=43)

                if rank == 2:
                    data2_in = comm.recv(source=1, tag=12)
                    temp2_old = kitchen.temperature_matrix #needed to relax temperature calculations
                    kitchen.get_temperature_matrix(data2_in)
                    #relaxation step necessary for producing convergent solution.
                    kitchen.temperature_matrix = (0.8)*kitchen.temperature_matrix + (0.2)*temp2_old
                    data2_out = kitchen.get_neumann_temps()
                    comm.send(data2_out, dest=1, tag=21)

                if rank == 3:
                    data31_in = comm.recv(source=1, tag=13)
                    data34_in = comm.recv(source=4, tag=43)
                    temp3_old = entryway.temperature_matrix #needed to relax temperature calculations
                    entryway.get_temperature_matrix(data34_in, data31_in)
                    #relaxation step necessary for producing convergent solution.
                    entryway.temperature_matrix = (0.8)*entryway.temperature_matrix + (0.2)*temp3_old
                    data31_out, data34_out = entryway.get_neumann_temps()
                    comm.send(data31_out, dest=1, tag=31)
                    comm.send(data34_out, dest=4, tag=34)
                i += 1
        #due to blocking communication processes 1 and 4 must first receive
        if rank == 1:
            btemp2 = comm.recv(source=2, tag=21)
            btemp3 = comm.recv(source=3, tag=31)
            send_temp1 = livingroom.temperature_matrix[1:-1,1:-1] #remove exterior points
            comm.send(send_temp1, dest=0, tag=10)

        if rank == 2:
            send_temp2 = kitchen.temperature_matrix[1:-1,1:-1] #remove exterior points
            comm.send(send_temp2, dest=0, tag=20)

        if rank == 3:
            send_temp3 = entryway.temperature_matrix[1:-1,1:-1] #remove exterior points
            comm.send(send_temp3, dest=0, tag=30)

        if rank == 4:
            btemp4 = comm.recv(source=3, tag=34)
            send_temp4 = bathroom.temperature_matrix[1:-1,1:-1] #remove exterior points
            comm.send(send_temp4, dest=0, tag=40)

        if rank == 0:
            om1 = comm.recv(source=1, tag=10)
            om2 = comm.recv(source=2, tag=20)
            om3 = comm.recv(source=3, tag=30)
            om4 = comm.recv(source=4, tag=40)

        return om1, om2, om3, om4


class Problem(Solver):

    def __init__(self, heater, aircon, walls):
        """
        Sets up the problem parameters.

        Params:
        -------
        heater : float
            the value prescribed to walls with a Heater
        aircon : float
            the value prescribed to walls with a Window
        wall : float
            the value prescribed to walls without a window or heater

        OM1 : ndarray
            computed temperature from the living room

        OM2 : ndarray
            computed temperature from the kitchen

        OM3 : ndarray
            computed temperature from the entryway/mudroom

        OM4 : ndarray
            computed temperature from the bathroom

        """
        self.heater = heater
        self.aircon = aircon
        self.wall = walls



    def __call__(self, cols, iters, open = False, on_off = False):
        """
        Performs the algorithm and produces the solutions.

        Params:
        -------
        cols : int
            The number of columns of interior gridpoints to be solved for.

        iters : int
            The number of iterations for the algorithm to perform.

        open : bool
            Boolean value passed into the problem determin whether or not the patio
            door is open(True) or closed(False) as default.

        on_off : bool
            Boolean value indicating whether the oven/stove is on or off thereby
            acting as a source of heat or a source of cold due to constant vaccum
            being applied at that location.

        Returns:
        -------
        None

        """
        self.open = open
        self.on_off = on_off
        self.OM1, self.OM2, self.OM3, self.OM4 = self.dirichelt_neumann_iteration(self.heater, self.aircon, self.wall, cols, iters, open, on_off)

    def img_creator(self):
        """
        Takes the final solutions from __call__ that performs the iterative algorithm
        for the 2D heat equation using parallel processing.

        Params:
        -------
        None

        Returns:
        --------
        None

        """
        apartment = plot_domain.Plotter(self.OM1, self.OM2, self.OM3, self.OM4, self.wall)
        apartment.room_image(self.heater,self.open,self.on_off)





"""
mpirun -n 5 python3 problem_solver.py
"""
tester = Problem(35, 8, 22)
tester(20,10,True,True)
tester.img_creator()

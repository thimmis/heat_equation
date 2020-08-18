#!/usr/bin/env python3
# @Author: Thomas Turner <thomas>
# @Date:   2020-08-03T12:07:54+02:00
# @Email:  thomas.benjamin.turner@gmail.com
# @Last modified by:   thomas
# @Last modified time: 2020-08-18T11:14:21+02:00



from mpi4py import MPI
from numpy.linalg import solve
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import *
import room_kitchen, room_bathroom, room_livingroom, room_entryway, plot_domain

class Solver:


    def dirichelt_neumann_iteration(self, heater, aircon, walls, cols, iters, open):
        """Solves the 2D-Heat Equation iteratively perscribing Dirichlet and
        Neumann boundary conditions to the different domains.

        Uses 5 processes (0-4) to produce a convergent solution.



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
        iterations = iters


        livingroom = room_livingroom.LivingRoom(heater, aircon, walls, cols, open)
        kitchen = room_kitchen.Kitchen(heater, aircon, walls, cols, open)
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
                    temp2_old = kitchen.temperature_matrix
                    data2_in = comm.recv(source=1, tag=12)
                    kitchen.get_temperature_matrix(data2_in)
                    kitchen.temperature_matrix = 0.8*kitchen.temperature_matrix + 0.2*temp2_old
                    data2_out = kitchen.get_neumann_temps()
                    comm.send(data2_out, dest=1, tag=21)

                if rank == 3:
                    temp3_old = entryway.temperature_matrix
                    data31_in = comm.recv(source=1, tag=13)
                    data34_in = comm.recv(source=4, tag=43)
                    entryway.get_temperature_matrix(data34_in, data31_in)
                    entryway.temperature_matrix = 0.8*entryway.temperature_matrix + 0.2*temp3_old
                    data31_out, data34_out = entryway.get_neumann_temps()
                    comm.send(data31_out, dest=1, tag=31)
                    comm.send(data34_out, dest=4, tag=34)

                i +=1
            else:
                if rank == 1:
                    btemp2 = comm.recv(source=2, tag=21)
                    btemp3 = comm.recv(source=3, tag=31)
                    g_13, g_12 = livingroom.temp_gradient_calc(btemp3, -btemp2)
                    comm.send(g_12, dest=2, tag=12)
                    comm.send(g_13, dest=3, tag=13)


                if rank == 4:
                    btemp4 = comm.recv(source=3, tag=34)
                    g_43 = bathroom.temp_gradient_calc(btemp4)
                    comm.send(g_43, dest=3, tag=43)

                if rank == 2:
                    data2_in = comm.recv(source=1, tag=12)
                    temp2_old = kitchen.temperature_matrix
                    kitchen.get_temperature_matrix(data2_in)
                    kitchen.temperature_matrix = (0.8)*kitchen.temperature_matrix + (0.2)*temp2_old
                    data2_out = kitchen.get_neumann_temps()
                    comm.send(data2_out, dest=1, tag=21)

                if rank == 3:
                    data31_in = comm.recv(source=1, tag=13)
                    data34_in = comm.recv(source=4, tag=43)
                    temp3_old = entryway.temperature_matrix
                    entryway.get_temperature_matrix(data34_in, data31_in)
                    entryway.temperature_matrix = (0.8)*entryway.temperature_matrix + (0.2)*temp3_old
                    data31_out, data34_out = entryway.get_neumann_temps()
                    comm.send(data31_out, dest=1, tag=31)
                    comm.send(data34_out, dest=4, tag=34)
                i += 1
        if rank == 1:
            btemp2 = comm.recv(source=2, tag=21)
            btemp3 = comm.recv(source=3, tag=31)
            livingroom.temp_gradient_calc(btemp3,btemp2)
            send_temp1 = livingroom.temperature_matrix[1:-1,1:-1]
            comm.send(send_temp1, dest=0, tag=10)

        if rank == 2:
            send_temp2 = kitchen.temperature_matrix[1:-1,1:-1]
            comm.send(send_temp2, dest=0, tag=20)

        if rank == 3:
            send_temp3 = entryway.temperature_matrix[1:-1,1:-1]
            comm.send(send_temp3, dest=0, tag=30)

        if rank == 4:
            btemp4 = comm.recv(source=3, tag=34)
            bathroom.temp_gradient_calc(btemp4)
            send_temp4 = bathroom.temperature_matrix[1:-1,1:-1]
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
        self.OM1, self.OM2, self.OM3, self.OM4 = self.dirichelt_neumann_iteration(self.heater, self.aircon, self.wall, cols, iters, open)

    def img_creator(self,resolution):
        apartment = plot_domain.Plotter(self.OM1, self.OM2, self.OM3, self.OM4, self.wall)
        apartment.room_image(resolution)


tester = Problem(40, 10, 25)
tester(20,10)
tester.img_creator(20)

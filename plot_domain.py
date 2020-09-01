#!/usr/bin/env python3
# @Author: Thomas Turner <thomas>
# @Date:   2020-08-03T12:07:54+02:00
# @Email:  thomas.benjamin.turner@gmail.com
# @Last modified by:   thomas
# @Last modified time: 2020-09-01T13:48:46+02:00

import matplotlib.pyplot as plt
from matplotlib import *
import numpy as np

class Plotter(object):
    """
    Takes the computed temperature solution matrices. Flips the individual
    matrices based on orientation and location, applies padding where needed
    and stacks them all together.

    Once the information has been pieced together it takes the heater temperature
    and applies that as the resolution to the filled contour plot, as well as
    the boolean arguments used and saves the figure as an image with the conditions
    as the suffix.

    Attributes:
    -----------
    room : ndarray
        The final matrix containing all temperature values in all rooms to be
        plotted.

    Methods:
    --------
    room_setup(self, om1, om2, om3, om4,temp)
        Stacks, pads and flips all of the different matrices into one.


    room_image(self,resolution,input1, input2)
        Plots the room and saves the image.



    """

    def __init__(self, om1, om2, om3, om4,temp):
        """
        Takes the individual room solutions and pieces them together using
        room_setup.

        Params:
        -------
        om1 : ndarray
            The solution matrix computed on process 1 (living room).

        om2 : ndarray
            The solution matrix computed on process 2 (kitchen).

        om3 : ndarray
            The solution matrix computed on process 3 (entryway).

        om4 : ndarray
            The solution matrix computed on process 4 (bathroom).

        temp : int
            The wall temperature to determine the padding where the closets
            are located.

        """
        self.room = self.room_setup(om1, om2, om3, om4,temp)

    def room_setup(self, om1, om2, om3, om4,temp):
        """
        Takes the computed temperature solution matrices. Flips the individual
        matrices based on orientation and location, applies padding where needed
        and stacks them all together.

        Params:
        -------
        om1 : ndarray
            The solution matrix computed on process 1 (living room).

        om2 : ndarray
            The solution matrix computed on process 2 (kitchen).

        om3 : ndarray
            The solution matrix computed on process 3 (entryway).

        om4 : ndarray
            The solution matrix computed on process 4 (bathroom).

        temp : int
            The wall temperature to determine the padding where the closets
            are located.

        Returns:
        --------
        apartment : ndarray
            The final solution matrix for all of the rooms.


        """
        tleft = np.hstack((om4,np.vstack((om3,np.full(om3.shape,temp)))))
        left = np.vstack((tleft,om2))
        apartment = np.hstack((left,om1))
        apartment = np.flip(apartment, 0)
        return apartment

    def room_image(self,resolution,input1, input2):
        """
        Creates and saves a filled contour plot of the domain. The resolution is
        determined by the maximum temperature input. The suffix of the image is
        appended with the boolean conditions passed in as arguments.

        Params:
        -------
        resolution : int
            The number used to determine how many different temperatures to
            show in the contour plot.

        input1 : str
            The boolean that tells if the door is open or closed.

        input2 : str
            The boolean that tells if the oven is on or off.

        Returns:
        --------
        None
        """
        suffix = str(input1)+'_'+str(input2) #force everything to be a string for concatenatation
        fig, ax2 = plt.subplots()
        CS = plt.contourf(self.room,resolution,cmap = plt.cm.inferno)
        Cbar = fig.colorbar(CS)
        plt.savefig('Temperature_' + suffix +'.png')

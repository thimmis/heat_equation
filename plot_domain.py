#!/usr/bin/env python3
# @Author: Thomas Turner <thomas>
# @Date:   2020-08-03T12:07:54+02:00
# @Email:  thomas.benjamin.turner@gmail.com
# @Last modified by:   thomas
# @Last modified time: 2020-08-17T14:58:01+02:00

import matplotlib.pyplot as plt
from matplotlib import *
import numpy as np
from datetime import datetime

class Plotter(object):
    """

    Attributes:
    -----------

    Methods:
    --------


    """

    def __init__(self, om1, om2, om3, om4,temp):
        """


        Params:
        -------

        """
        self.room = self.room_setup(om1, om2, om3, om4,temp)

    def room_setup(self, om1, om2, om3, om4,temp):
        """

        Params:
        -------

        Returns:
        --------


        """
        x,toss = om2.shape
        tleft = np.hstack((om4,np.vstack((om3,np.full(om3.shape,temp)))))
        left = np.vstack((tleft,om2))
        apartment = np.hstack((left,om1))
        apartment = np.flip(apartment, 0)
        return apartment

    def room_image(self,resolution):
        now = datetime.now()
        suffix = now.strftime("%Y%m%d%M%S")
        fig, ax2 = plt.subplots()
        CS = plt.contourf(self.room,resolution,cmap = plt.cm.inferno)
        Cbar = fig.colorbar(CS)
        plt.savefig('Temperature_' + suffix +'.png')

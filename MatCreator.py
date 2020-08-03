# @Author: Thomas Turner <thomas>
# @Date:   2020-06-17T13:27:38+02:00
# @Email:  thomas.benjamin.turner@gmail.com
# @Last modified by:   thomas
<<<<<<< current
# @Last modified time: 2020-08-03T12:02:04+02:00
=======
# @Last modified time: 2020-08-03T11:58:49+02:00
>>>>>>> before discard


from scipy import *
import numpy as np


class MCreate:

    def __init__(self, dimension):
        self.dimension = int(dimension/3)
        self.A = self.kronstamp(self.dimension)

    def kronstamp(self,dimension):
        S = np.array([[-4,1,0],[1,-4,1],[0,1,-4]])
        return np.kron(np.eye(dimension),S) + np.kron(np.diag(np.ones(dimension-1),-1),np.eye(3)) +np.kron(np.diag(np.ones(dimension-1),1),np.eye(3))

# HEAT EQUATION SOLVER

## Motivation
While living in Sweden during my Master's degree I noticed that our meager student accomodations would get very cold, so I wrote a steady state heat equation solver to model the temperature of our apartment with different boundary conditions.

## How does it work?
The model uses MPI4py to discretize the domain into different sections. With blocking we iteratively apply Neumann/Dirichlet boundary conditions until we reach a convergent solution, i.e. smooth solution curves across boundaries.

### Selecting boundary conditions:
To run the solver we can pass in the boundary condition temperatures as well as additional flags that allow us to control certain aspects of the boundaries. In particular we can set the external wall to a lower temperature,  we can open the windows to provide a draft or we can turn on the oven and bake something to understand how the diffuse temeperature varies.

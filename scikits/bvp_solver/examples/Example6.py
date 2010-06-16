# generated using print bvp_solver.get_template(num_ODE = 5, num_parameters =0, num_left_boundary_conditions = 4, boundary_conditions_derivative = True)
import scikits.bvp_solver
import numpy
from numpy import array

#all None values in the code below are dummy values which must be replaced with real values

def function( X, y):


    Y3MY1 =  y[2] - y[0]
    Y3MY5 = y[2] - y[4]
    return array([.5 * y[0] * Y3MY1 / y[1]      , #evaluate ODE number 0
                  -.5 * Y3MY1      , #evaluate ODE number 1
                  (.9 - 1000. * Y3MY5 - .5 * y[2] * Y3MY1) / y[3]      , #evaluate ODE number 2
                  .5 * Y3MY1      , #evaluate ODE number 3
                  100. * Y3MY5      ]) #evaluate ODE number 4

def boundary_conditions(Ya, Yb):

    BCa = array([Ya[0] - 1.0      , #evaluate left BC number 0
                 Ya[1] - 1.0      , #evaluate left BC number 1
                 Ya[2] - 1.0      , #evaluate left BC number 2
                 Ya[3] + 10.0      ]) #evaluate left BC number 3

    BCb = array([Yb[2] - Yb[4]      ]) #evaluate right BC number 0

    return BCa, BCb


def boundary_conditions_derivative( Ya, Yb):

    #evaluate left boundary conditions derivative with respect to variables
    #increasing differentiation index
    # (Ya,i) ----->
    dBCadYa = array([[1.0      ,0      ,0      ,0      ,0      ], #left BC number 0
                     [0      ,1.0      ,0      ,0      ,0      ], #left BC number 1
                     [0      ,0      ,1.0      ,0      ,0      ], #left BC number 2
                     [0      ,0      ,0      ,1.0      ,0      ]]) #left BC number 3

    #evaluate right boundary conditions derivative with respect to variables
    #increasing differentiation index
    # (Yb,i) ----->
    dBCbdYb = array([[0      ,0      ,1.0      ,0      ,-1      ]]) #right BC number 0

    return dBCadYa, dBCbdYb

def guess_y (X):
    return numpy.array([1.0,
                        1.0,
                        1.0 + 8.9 * X - 4.5 * X**2,
                        -10.0,
                        .91 + 9 * X - 4.5 * X**2])

problem_definition = scikits.bvp_solver.ProblemDefinition(num_ODE = 5,
                                                  num_parameters = 0,
                                                  num_left_boundary_conditions = 4,
                                                  boundary_points = (0.0  , 1.0),
                                                  function = function,
                                                  boundary_conditions = boundary_conditions,
                                                  boundary_conditions_derivative = boundary_conditions_derivative)
solution = scikits.bvp_solver.solve(bvp_problem = problem_definition,
                            solution_guess = guess_y)

import pylab

x = numpy.linspace(0,1.0)
pylab.plot(x, solution(x)[0,:],'-')
pylab.plot(x, solution(x)[1,:],'-')
pylab.plot(x, solution(x)[2,:],'-')
pylab.plot(x, solution(x)[3,:],'-')
pylab.plot(x, solution(x)[4,:],'-')

pylab.show()
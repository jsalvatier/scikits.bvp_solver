# based on the Fortran code found at http://cs.stmarys.ca/~muir/BVP_SOLVER_Files/swave.f90

'''
Created on Apr 29, 2009

@author: johnsalvatier
'''
from __future__ import division
import scikits.bvp_solver
import numpy
import pylab
from numpy import array

#all None values in the code below are dummy values which must be replaced with real values

eps = array(.01)
gamma = array(1.4)

def function( t, y):


    term1 = 1.0/ (eps * (1.0 + t**2))
    term2 = .5 + .5 * gamma - 2.0 * eps * t
    
    F = array([y[1]      , #evaluate ODE number 0
                  (term1/y[0]) * (term2 * y[0] * y[1] - y[1]/y[0] -
                                  (2.0 * t/(1.0 + t**2)) * (1.0 - .5 *(gamma - 1.0) * y[0] **2))       ]) #evaluate ODE number 1
        
    return F 

def boundary_conditions(Ya, Yb):

    BCa = array([Ya[0] - .9129      ]) #evaluate left BC number 0

    BCb = array([Yb[0] - .375      ]) #evaluate right BC number 0
    return BCa, BCb

def function_derivative( t, y):

    term1 = 1.0/ (eps * (1.0 + t**2))
    term2 = .5 + .5 * gamma - 2.0 * eps * t

    #evaluate function derivative with respect to variables    #increasing differentiation index
    #(Yi) ---->

    dFdY = array([[0.0      ,1.0      ], #ODE number 0
                  [term1*( 2.0*y[1]/(y[0]**3) + 2.0*t/((1.0+t**2)*y[0]**2) + (t/(1.0+t**2))*(gamma-1.0))      ,
                   (term1/y[0])*( term2*y[0] - 1.0/y[0])      ]]) #ODE number 1
    
    return dFdY

def boundary_conditions_derivative( Ya, Yb):

    #evaluate left boundary conditions derivative with respect to variables
    #increasing differentiation index
    # dBC1/dYa1----->dBC1/dYai
    dBCadYa = array([[1.0      ,0.0      ]]) #left BC number 0

    #evaluate right boundary conditions derivative with respect to variables
    #increasing differentiation index
    # (Yb,i) ----->
    dBCbdYb = array([[1.0      ,0.0      ]]) #right BC number 0

    return dBCadYa, dBCbdYb


problem_definition = scikits.bvp_solver.ProblemDefinition(num_ODE = 2,
                                                  num_parameters = 0,
                                                  num_left_boundary_conditions = 1,
                                                  boundary_points = (0.0, 1.0),
                                                  function = function,
                                                  boundary_conditions = boundary_conditions,
                                                  function_derivative = function_derivative,
                                                  boundary_conditions_derivative = boundary_conditions_derivative)


slope = array(.375 - .9129)
x = array([0.0 ,  0.11111,  0.22222,  0.33333,  0.44444,  0.55555,  0.66666,  0.77778,  0.88888,  1.0   ]) #
#x = numpy.linspace(0,1,20)


solution = scikits.bvp_solver.solve(bvp_problem = problem_definition,
                                        initial_mesh = x,
                                        solution_guess = [ .9129 + slope *x, x*0 + slope     ],
                                        trace = 0)


points = numpy.linspace(0,1, 200)
pylab.plot(points, solution(points)[0,:],'-')
pylab.plot(points, solution(points)[1,:],'-')
pylab.show()


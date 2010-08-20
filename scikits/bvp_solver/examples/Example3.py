#based off of this Fortran example http://cs.stmarys.ca/~muir/BVP_SOLVER_Files/lubrication.f90
# unknown parameters
# function derivative callback
# boundary conditions derivative callback
# saving of a solution

import scikits.bvp_solver
import numpy
import pylab

"""
 This is Example 3.5.2 of Solving ODEs with Matlab, a nonlinear
 eigenvalue problem of lubrication theory studied in section
 6.1 of H.B. Keller, Numerical Methods for Two-Point Boundary-
 Value Problems, Dover, New York, 1992.
 
 Problem specification:
 y' = (sin(x)**2 - p * sin(x) **4/y)/eps
 0 = y(0) - 1
 0 = y(L) - 1
 
"""

print ("test #2 (unknown parameters)")
eps = .1
def function2(X , Y, P):
    return numpy.array([ ( numpy.sin(X)**2 - P[0]*numpy.sin(X)**4/Y[0] ) / eps])


def dfunction2(X,Y, P):

    dFdY = numpy.zeros((1,1))
    dFdY[0,0] = ( P[0]*numpy.sin(X)**4/Y[0]**2 ) / eps

    dFdP = numpy.zeros((1,1))
    dFdP[0,0] = (- numpy.sin(X)**4/Y[0] ) / eps

    return dFdY, dFdP

def boundary_conditions2(YA,YB, P):

    BCA= numpy.zeros(1)
    BCB= numpy.zeros(1)

    BCA[0] = YA[0] - 1.0
    BCB[0] = YB[0] - 1.0

    return BCA, BCB

def dbconditions2(YA,YB,P):

    DYA = numpy.zeros((1,1))
    DYA[0,0] = 1

    DYB = numpy.zeros((1,1))
    DYB[0,0] = 1

    DAP = numpy.zeros((1,1))
    DBP = numpy.zeros((1,1))
    return DYA, DYB, DAP, DBP


problem2 =scikits.bvp_solver.ProblemDefinition(num_ODE = 1,
                                      num_parameters = 1,
                                      num_left_boundary_conditions = 1,
                                      boundary_points = (-numpy.pi/2.0, numpy.pi/2.0),
                                      function = function2,
                                     boundary_conditions = boundary_conditions2,
                                     function_derivative = dfunction2,
                                     boundary_conditions_derivative = dbconditions2)

solution =scikits.bvp_solver.solve(problem2,
                            solution_guess = .5,
                            parameter_guess = 1,
                            initial_mesh = numpy.linspace(problem2.boundary_points[0],problem2.boundary_points[1], 21),
                            trace =2)
print (solution.solution)
print solution.parameters
pylab.plot(solution.mesh, solution.solution[0,:],'-')

pylab.show()
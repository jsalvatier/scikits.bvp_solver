#based off of this fortran example http://cs.stmarys.ca/~muir/BVP_SOLVER_Files/mathieu.f90
# use of function derivative callback
# function callback for a guess
# a singular term

import scikits.bvp_solver
import numpy
import pylab
print ("example #4")

"""
 This is the MAT4BVP example of Matlab.  It finds the fourth
 eigenvalue of Mathieu's equation

      y'' + (lambda - 2*q*cos(2*x))*y = 0
   
 on the interval [0, pi] with boundary conditions y'(0) = 0, 
 y'(pi) = 0 when the parameter q = 5.

 Special codes for Sturm-Liouville problems can compute specific
 eigenvalues. The general-purpose code BVP_SOLVER can only compute 
 an eigenvalue near to a guessed value. We can make it much more 
 likely that we compute the eigenfunction corresponding to the fourth 
 eigenvalue by supplying a guess that has the correct qualitative 
 behavior. The eigenfunction y(x) is determined only to a constant 
 multiple, so the normalizing condition y(0) = 1 is used to specify 
 a particular solution.
"""


Q = 5.0

def function(X,Y,P):

    LAMBDA = P[0]
    return numpy.array([Y[1],
                        -(LAMBDA - 2.0*Q*numpy.cos(2.0*X))*Y[0]])

def function_derivative(X,Y,P):

    LAMBDA = P[0]
    dFdY = numpy.array([[0.0                              ,1.0 ],
                       [-(LAMBDA - 2.0*Q*numpy.cos(2.0*X)) ,0.0 ]])
    dFdP = numpy.array([[0],
                        [-Y[0]]])

    return dFdY, dFdP


def boundary_conditions (Ya,Yb,P):
    BCa = numpy.array([Ya[1],
                       Ya[0] - 1.0])

    BCb = numpy.array([Yb[1]])
    return BCa, BCb

def boundary_conditions_derivative(YA,YB,P):

    dBCdYa = numpy.array([[ 0.0 , 1.0],
                         [ 1.0 , 0.0]])

    dBCdYb = numpy.array([[ 0.0 , 1.0]])

    dBCAdP = numpy.zeros((2,1))
    dBCBdP = numpy.zeros((1,1))

    return dBCdYa, dBCdYb, dBCAdP, dBCBdP


problem = scikits.bvp_solver.ProblemDefinition(num_ODE = 2,
                                  num_parameters = 1,
                                  num_left_boundary_conditions = 2,
                                  boundary_points = (0, numpy.pi),
                                  function = function,
                                  boundary_conditions = boundary_conditions,
                                  function_derivative = function_derivative,
                                  boundary_conditions_derivative = boundary_conditions_derivative)

def guess(X):

    return numpy.array([numpy.cos(4.0*X)   ,
                        -4.0*numpy.sin(4.0*X)])

solution = scikits.bvp_solver.solve(problem,
                            solution_guess = guess,
                            parameter_guess = numpy.array([15.0]),
                            trace = 2)

x = numpy.linspace(problem.boundary_points[0],problem.boundary_points[1], 100)
y = solution(x)
print "lambda = " + str(solution.parameters[0])
print "mat4bvp gives lambda ~= 17.097"

pylab.plot(x, y[0,:],'-')
pylab.plot(x, y[1,:],'-')

pylab.show()
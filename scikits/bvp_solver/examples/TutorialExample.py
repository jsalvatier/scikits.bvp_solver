"""
First we import scikits.bvp_solver and numpy because the callback functions must return numpy arrays.
"""
import scikits.bvp_solver
import numpy

# next we define the important constants
T10 = 130
T2Ahx = 70
Ahx = 5
U = 1.0

"""
Then we define the callback function which evaluates the ODEs. The first argument of this function gets the independent variable,
in our case A, and the second argument gets an array of the values of the dependent variables, in our case, T_1 and T_2.
If the boundary value problem includes unknown parameters, then a third argument gets the values of the unknown parameters,
but since this problem does not deal with unknown parameters, the function can only have 2 arguments. In this problem
we will have the first variable of the dependent variable array represent stream 1 (the hot liquid).

This function must return a numpy array that contains the value of the first derivative of each dependent variable. These
must be in the same order as the dependent variable array.
"""

def function(a , T):
    q = (T[0] - T[1]) * U           # calculate the heat transfer from stream 1 to stream 2
    return numpy.array([-q ,        # evaluate dT1/dA
                        q/-2.0])    # evaluate dT2/dA

"""
To finish the problem definition, we define the define the callback function which evaluates the difference between the actual
boundary conditions and the required boundary conditions. The first argument this function receives is an array of the values
of the dependent variables (here T_1 and T_2) at the left boundary condition. The second argument this function receives is
an array of the values of the dependent variables (here T_1 and T_2) at the right boundary condition.

This function must return two numpy arrays that contains the difference between the actual boundary conditions and the required
boundary conditions. The first return array must contain these differences for all boundary conditions on the left, and second
return array must contain these differences for all boundary conditions on the right. The sizes of these arrays must add up to
the total number of ODEs plus the number of unknown parameters. Both these arrays must be in the same order as the dependent
variable arrays (which is the same as in the first function).
"""
def boundary_conditions(Ta,Tb):

    return (numpy.array([Ta[0] - T10]),  #evaluate the difference between the temperature of the hot stream on the
                                         #left and the required boundary condition
            numpy.array([Tb[1] - T2Ahx]))#evaluate the difference between the temperature of the cold stream on the
                                         #right and the required boundary condition

"""
Next we create the ProblemDefinition object, by passing the relevant information and callbacks to its constructor.
"""
problem = scikits.bvp_solver.ProblemDefinition(num_ODE = 2,
                                      num_parameters = 0,
                                      num_left_boundary_conditions = 1,
                                      boundary_points = (0, Ahx),
                                      function = function,
                                      boundary_conditions = boundary_conditions)

solution = scikits.bvp_solver.solve(problem,
                            solution_guess = ((T10 + T2Ahx)/2.0,
                                              (T10 + T2Ahx)/2.0))
"""
The "solve" function returns a Solution object which can be passed an array of points at which to evaluate the solution.
"""
A = numpy.linspace(0,Ahx, 45)
T = solution(A)
print T

"""
We can plot the solution using pylab with the following code
"""
import pylab
pylab.plot(A, T[0,:],'-')
pylab.plot(A, T[1,:],'-')
pylab.show()

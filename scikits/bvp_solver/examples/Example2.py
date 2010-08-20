#based off of this fortran example http://cs.stmarys.ca/~muir/BVP_SOLVER_Files/esi.f90
# use of function derivative callback
# function callback for a guess
# a singular term

import scikits.bvp_solver
import numpy
import pylab
print ("test #1")

"""
 U.M. Ascher and R.D. Russell, Reformulation of boundary
 value problems into `standard' form, SIAM Review 23 (1981)
 238-254 use a BVP from electromagnetic self-interaction
 theory to discuss problems with singular coefficients set
 on an infinite interval.  After some preparation they solve

    u'' + 4*u'/t + (t*u - 1)*u = 0

    u'(0) = 0 
    u(L) + u'(L) = 0
    
 The problem is set on
 an infinite interval, so some experimentation is necessary
 to verify that a sufficiently large L has been specified.
 They present results for u(0) when L = 5,8,10,20.  They
 use the initial guess u(t) = 2 for 0 <= t <= 1.5 and 
 u(t) = 2*exp(1.5 - t) for t > 1.5.  We use this guess for
 the first L and thereafter use the solution for one L as
 the guess for the next, extending to the right with the
 value from their guess (which has the right asymptotic
 behavior).
 """



def function1(X , Y):
    return numpy.array([Y[1], -(X*Y[0] - 1.0)*Y[0]])

def dfunction1(X,Y):
    PD = numpy.zeros((2,2))
    PD[0,1] = 1.0
    PD[1,0] = 1.0 - 2*X*Y[0]
    return PD

def boundary_conditions1(YA,YB):
    BCA= numpy.zeros(1)
    BCB= numpy.zeros(1)

    BCA[0] = YA[1]
    BCB[0] = YB[0] + YB[1]

    return BCA, BCB

def guess_y1(X):
    Y = numpy.zeros(2)

    if X <= 1.5:
      Y[0] = 2
      Y[1] = 0
    else:
      Y[0] = 2.0*numpy.exp(1.5 - X)
      Y[1] = - Y[0]

    return Y

singular_term = numpy.zeros((2,2))
singular_term[1,1] = -4.0


L = [5, 8 , 10, 20]

solutionList = []
problemList = []
for i in range(4):
    problemList.append(scikits.bvp_solver.ProblemDefinition(num_ODE = 2,
                                    num_parameters = 0,
                                    num_left_boundary_conditions = 1,
                                    boundary_points = (0, L[i]),
                                    function = function1,
                                    boundary_conditions = boundary_conditions1,
                                    function_derivative = dfunction1))
    if i == 0:
        solutionList.append(scikits.bvp_solver.solve(problemList[i],
                                         solution_guess = guess_y1,
                                         singular_term = singular_term,
                                         max_subintervals = 300))
    else:
        solutionList.append(solutionList[i-1].extend(0, L[i]))
        solutionList[i] = scikits.bvp_solver.solve(problemList[i],
                                         solution_guess = solutionList[i],
                                         singular_term = singular_term,
                                         max_subintervals = 300)

for i in range(4):
    x = numpy.linspace(problemList[i].boundary_points[0],problemList[i].boundary_points[1], 45)
    y = solutionList[i](x)

    pylab.subplot(1, 4, i + 1)
    pylab.plot(x, y[0,:],'-')
    pylab.plot(x, y[1,:],'-')

pylab.show()

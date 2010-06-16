#based off of this fortran example http://cs.stmarys.ca/~muir/BVP_SOLVER_Files/esi.f90

# Tests by confirming that the results match those reported in BVP_SOLVER's example
# USES:
# use of function derivative callback
# function callback for a guess
# a singular term

import scikits.bvp_solver
import numpy
from numpy.testing import assert_almost_equal
import nose
import example2data as data
import testing

#Results obtained by Ascher/Russell.

AR = numpy.array([2.03464, 2.11986, 2.11997, 2.11997 ])

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


def test_example2():

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
                                             singular_term = singular_term))
        else:
            solutionList.append(solutionList[i-1].extend(0, L[i]))

            #make sure you can't evaluate an extended solution
            nose.tools.assert_raises(ValueError, solutionList[i], numpy.array([0]) )

            solutionList[i] = scikits.bvp_solver.solve(problemList[i],
                                             solution_guess = solutionList[i],
                                             singular_term = singular_term)
        if i == 2:

            testing.assert_solution_matches_data_eval(solutionList[i], data, 4)
            ## possibly look into this, the mesh does not line up correctly
            testing.assert_solution_matches_data_calculated(solutionList[i], data, 2)

    for i in range(4):
        # evaluate at 0 and compare with  Ascher/Russell
        # assert that all the values are the same
        assert_almost_equal(solutionList[i](0.0)[0,0], AR[i], 4)
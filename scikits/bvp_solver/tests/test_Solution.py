#based off of this Fortran example http://cs.stmarys.ca/~muir/BVP_SOLVER_Files/lubrication.f90
# unknown parameters
# function derivative callback
# boundary conditions derivative callback
# saving of a solution

import scikits.bvp_solver
import numpy
from numpy.testing import assert_almost_equal
import nose

def functionGood(X , Y, P):
    return numpy.array([ 0])


def dfunctionGood(X,Y, P):

    dFdY = numpy.zeros((1,1))
    dFdP = numpy.zeros((1,1))

    return dFdY, dFdP

def boundary_conditionsGood(YA,YB, P):

    BCA= YA[0] - numpy.ones(1)
    BCB= YB[0] - numpy.zeros(1)

    return BCA, BCB


def dbconditionsGood(YA,YB,P):

    DYA = numpy.ones((1,1))

    DYB = numpy.ones((1,1))

    DAP = numpy.zeros((1,1))
    DBP = numpy.zeros((1,1))

    return DYA, DYB, DAP, DBP



def test_failure():
    problem1 = scikits.bvp_solver.ProblemDefinition(num_ODE = 1,
                                        num_parameters = 1,
                                        num_left_boundary_conditions = 1,
                                        boundary_points = (-numpy.pi/2.0, numpy.pi/2.0),
                                        function = functionGood,
                                        boundary_conditions = boundary_conditionsGood,
                                        function_derivative = dfunctionGood,
                                        boundary_conditions_derivative = dbconditionsGood)

    #should fail because it is impossible to solve this problem
    nose.tools.assert_raises(ValueError, scikits.bvp_solver.solve, problem1,
                                solution_guess = 0,
                                initial_mesh = numpy.linspace(problem1.boundary_points[0],problem1.boundary_points[1], 21),
                                parameter_guess = 1)

    solution = scikits.bvp_solver.solve(problem1,
                                solution_guess = 0,
                                initial_mesh = numpy.linspace(problem1.boundary_points[0],problem1.boundary_points[1], 21),
                                parameter_guess = 1,
                                error_on_fail = False)

    nose.tools.assert_raises(ValueError, solution, 0)
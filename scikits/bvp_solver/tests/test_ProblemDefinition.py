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

def dfunctionBad1(X,Y, P):

    dFdY = numpy.ones((1,1))
    dFdP = numpy.zeros((1,1))

    return dFdY, dFdP

def dfunctionBad2(X,Y, P):

    dFdY = numpy.zeros((1,1))
    dFdP = numpy.ones((1,1))

    return dFdY, dFdP

def boundary_conditionsGood(YA,YB, P):

    BCA= numpy.zeros(1)
    BCB= numpy.zeros(1)

    return BCA, BCB


def dbconditionsGood(YA,YB,P):

    DYA = numpy.zeros((1,1))

    DYB = numpy.zeros((1,1))

    DAP = numpy.zeros((1,1))
    DBP = numpy.zeros((1,1))

    return DYA, DYB, DAP, DBP

def dbconditionsBad1(YA,YB,P):

    DYA = numpy.ones((1,1))

    DYB = numpy.zeros((1,1))

    DAP = numpy.zeros((1,1))
    DBP = numpy.zeros((1,1))

    return DYA, DYB, DAP, DBP

def dbconditionsBad2(YA,YB,P):

    DYA = numpy.zeros((1,1))

    DYB = numpy.ones((1,1))

    DAP = numpy.zeros((1,1))
    DBP = numpy.zeros((1,1))

    return DYA, DYB, DAP, DBP

def dbconditionsBad3(YA,YB,P):

    DYA = numpy.zeros((1,1))

    DYB = numpy.zeros((1,1))

    DAP = numpy.ones ((1,1))
    DBP = numpy.zeros((1,1))

    return DYA, DYB, DAP, DBP

def dbconditionsBad4(YA,YB,P):

    DYA = numpy.zeros((1,1))

    DYB = numpy.zeros((1,1))

    DAP = numpy.zeros((1,1))
    DBP = numpy.ones((1,1))

    return DYA, DYB, DAP, DBP

def functionNoParamGood(X , Y):
    return numpy.array([ 0])


def dfunctionNoParamGood(X,Y):

    dFdY = numpy.zeros((1,1))

    return dFdY
def dfunctionNoParamBad1(X,Y):

    dFdY = numpy.ones((1,1))


    return dFdY

def boundary_conditionsNoParamGood(YA,YB):

    BCA= numpy.zeros(1)
    BCB= numpy.zeros(0)

    return BCA, BCB
def boundary_conditionsNoParamGood2(YA,YB):

    BCA= numpy.zeros(0)
    BCB= numpy.zeros(1)

    return BCA, BCB

def dbconditionsNoParamGood(YA,YB):

    DYA = numpy.zeros((1,1))

    DYB = numpy.zeros((0,1))


    return DYA, DYB

def dbconditionsNoParamBad1(YA,YB):

    DYA = numpy.ones((1,1))

    DYB = numpy.zeros((0,1))


    return DYA, DYB

def dbconditionsNoParamBad2(YA,YB):

    DYA = numpy.zeros((0,1))

    DYB = numpy.ones((1,1))


    return DYA, DYB


def test_creation():
    scikits.bvp_solver.ProblemDefinition(num_ODE = 1,
                                    num_parameters = 1,
                                    num_left_boundary_conditions = 1,
                                    boundary_points = (-numpy.pi/2.0, numpy.pi/2.0),
                                    function = functionGood,
                                    boundary_conditions = boundary_conditionsGood,
                                    function_derivative = dfunctionGood,
                                    boundary_conditions_derivative = dbconditionsGood)


def test_paramPlusNoDerivs():
    scikits.bvp_solver.ProblemDefinition(num_ODE = 1,
                                    num_parameters = 1,
                                    num_left_boundary_conditions = 1,
                                    boundary_points = (-numpy.pi/2.0, numpy.pi/2.0),
                                    function = functionGood,
                                    boundary_conditions = boundary_conditionsGood)

def test_BadBCDerivatives():
    problem1 = scikits.bvp_solver.ProblemDefinition(num_ODE = 1,
                                        num_parameters = 1,
                                        num_left_boundary_conditions = 1,
                                        boundary_points = (-numpy.pi/2.0, numpy.pi/2.0),
                                        function = functionGood,
                                        boundary_conditions = boundary_conditionsGood,
                                        function_derivative = dfunctionGood,
                                        boundary_conditions_derivative = dbconditionsBad1)

    nose.tools.assert_raises(ValueError, scikits.bvp_solver.solve, problem1,
                                solution_guess = 0,
                                initial_mesh = numpy.linspace(problem1.boundary_points[0],problem1.boundary_points[1], 21),
                                parameter_guess = 1)

    problem2 = scikits.bvp_solver.ProblemDefinition(num_ODE = 1,
                                        num_parameters = 1,
                                        num_left_boundary_conditions = 1,
                                        boundary_points = (-numpy.pi/2.0, numpy.pi/2.0),
                                        function = functionGood,
                                        boundary_conditions = boundary_conditionsGood,
                                        function_derivative = dfunctionGood,
                                        boundary_conditions_derivative = dbconditionsBad2)

    nose.tools.assert_raises(ValueError, scikits.bvp_solver.solve, problem2,
                                solution_guess = 0,
                                initial_mesh = numpy.linspace(problem2.boundary_points[0],problem2.boundary_points[1], 21),
                                parameter_guess = 1)

    problem3 = scikits.bvp_solver.ProblemDefinition(num_ODE = 1,
                                        num_parameters = 1,
                                        num_left_boundary_conditions = 1,
                                        boundary_points = (-numpy.pi/2.0, numpy.pi/2.0),
                                        function = functionGood,
                                        boundary_conditions = boundary_conditionsGood,
                                        function_derivative = dfunctionGood,
                                        boundary_conditions_derivative = dbconditionsBad3)

    nose.tools.assert_raises(ValueError, scikits.bvp_solver.solve, problem3,
                                solution_guess = 0,
                                initial_mesh = numpy.linspace(problem3.boundary_points[0],problem3.boundary_points[1], 21),
                                parameter_guess = 1)

    problem4 = scikits.bvp_solver.ProblemDefinition(num_ODE = 1,
                                        num_parameters = 1,
                                        num_left_boundary_conditions = 1,
                                        boundary_points = (-numpy.pi/2.0, numpy.pi/2.0),
                                        function = functionGood,
                                        boundary_conditions = boundary_conditionsGood,
                                        function_derivative = dfunctionGood,
                                        boundary_conditions_derivative = dbconditionsBad4)

    nose.tools.assert_raises(ValueError, scikits.bvp_solver.solve, problem4,
                                solution_guess = 0,
                                initial_mesh = numpy.linspace(problem4.boundary_points[0],problem4.boundary_points[1], 21),
                                parameter_guess = 1)
def test_BadFuncDerivatives():
    problem1 = scikits.bvp_solver.ProblemDefinition(num_ODE = 1,
                                        num_parameters = 1,
                                        num_left_boundary_conditions = 1,
                                        boundary_points = (-numpy.pi/2.0, numpy.pi/2.0),
                                        function = functionGood,
                                        boundary_conditions = boundary_conditionsGood,
                                        function_derivative = dfunctionBad1,
                                        boundary_conditions_derivative = dbconditionsGood)

    nose.tools.assert_raises(ValueError, scikits.bvp_solver.solve, problem1,
                                solution_guess = 0,
                                initial_mesh = numpy.linspace(problem1.boundary_points[0],problem1.boundary_points[1], 21),
                                parameter_guess = 1)

    problem2 = scikits.bvp_solver.ProblemDefinition(num_ODE = 1,
                                        num_parameters = 1,
                                        num_left_boundary_conditions = 1,
                                        boundary_points = (-numpy.pi/2.0, numpy.pi/2.0),
                                        function = functionGood,
                                        boundary_conditions = boundary_conditionsGood,
                                        function_derivative = dfunctionBad2,
                                        boundary_conditions_derivative = dbconditionsGood)

    nose.tools.assert_raises(ValueError, scikits.bvp_solver.solve, problem2,
                                solution_guess = 0,
                                initial_mesh = numpy.linspace(problem2.boundary_points[0],problem2.boundary_points[1], 21),
                                parameter_guess = 1)

def test_BadBCDerivativesNoParam():
    problem = scikits.bvp_solver.ProblemDefinition(num_ODE = 1,
                                                   num_parameters =0,
                                        num_left_boundary_conditions = 1,
                                        boundary_points = (-numpy.pi/2.0, numpy.pi/2.0),
                                        function = functionNoParamGood,
                                        boundary_conditions = boundary_conditionsNoParamGood,
                                        function_derivative = dfunctionNoParamGood,
                                        boundary_conditions_derivative = dbconditionsNoParamBad1)

    nose.tools.assert_raises(ValueError, scikits.bvp_solver.solve, problem,
                                solution_guess = 0,
                                initial_mesh = numpy.linspace(problem.boundary_points[0],problem.boundary_points[1], 21))

    problem2 = scikits.bvp_solver.ProblemDefinition(num_ODE = 1,
                                                   num_parameters =0,
                                        num_left_boundary_conditions = 0,
                                        boundary_points = (-numpy.pi/2.0, numpy.pi/2.0),
                                        function = functionNoParamGood,
                                        boundary_conditions = boundary_conditionsNoParamGood2,
                                        function_derivative = dfunctionNoParamGood,
                                        boundary_conditions_derivative = dbconditionsNoParamBad2)

    nose.tools.assert_raises(ValueError, scikits.bvp_solver.solve, problem2,
                                solution_guess = 0,
                                initial_mesh = numpy.linspace(problem2.boundary_points[0],problem2.boundary_points[1], 21))

def test_BadFuncDerivativesNoParam():
    problem1 = scikits.bvp_solver.ProblemDefinition(num_ODE = 1,
                                                   num_parameters =0,
                                        num_left_boundary_conditions = 1,
                                        boundary_points = (-numpy.pi/2.0, numpy.pi/2.0),
                                        function = functionNoParamGood,
                                        boundary_conditions = boundary_conditionsNoParamGood,
                                        function_derivative = dfunctionNoParamBad1,
                                        boundary_conditions_derivative = dbconditionsNoParamGood)

    nose.tools.assert_raises(ValueError, scikits.bvp_solver.solve, problem1,
                                solution_guess = 0,
                                initial_mesh = numpy.linspace(problem1.boundary_points[0],problem1.boundary_points[1], 21))

def test_parameters():
    problem = scikits.bvp_solver.ProblemDefinition(num_ODE = 1,
                                        num_parameters = 1,
                                        num_left_boundary_conditions = 1,
                                        boundary_points = (-numpy.pi/2.0, numpy.pi/2.0),
                                        function = functionGood,
                                        boundary_conditions = boundary_conditionsGood,
                                        function_derivative = dfunctionGood,
                                        boundary_conditions_derivative = dbconditionsGood)

    nose.tools.assert_raises(ValueError, scikits.bvp_solver.solve, problem,
                                solution_guess = 0,
                                initial_mesh = numpy.linspace(problem.boundary_points[0],problem.boundary_points[1], 21),
                                parameter_guess = [1,2,3] )
    nose.tools.assert_raises(ValueError, scikits.bvp_solver.solve, problem,
                            solution_guess = 0,
                            initial_mesh = numpy.linspace(problem.boundary_points[0],problem.boundary_points[1], 21) )

def test_properties():
    problem = scikits.bvp_solver.ProblemDefinition(num_ODE = 1,
                                    num_parameters = 1,
                                    num_left_boundary_conditions = 1,
                                    boundary_points = (-numpy.pi/2.0, numpy.pi/2.0),
                                    function = functionGood,
                                    boundary_conditions = boundary_conditionsGood,
                                    function_derivative = dfunctionGood,
                                    boundary_conditions_derivative = dbconditionsGood)

    assert (problem.num_ODE == 1)
    assert (problem.num_parameters == 1)
    assert (problem.num_left_boundary_conditions == 1)
    assert_almost_equal(problem.boundary_points, numpy.array([-numpy.pi/2.0, numpy.pi/2.0]))
    assert(problem.function == functionGood)
    assert(problem.function_derivative == dfunctionGood)
    assert(problem.boundary_conditions == boundary_conditionsGood)
    assert(problem.boundary_conditions_derivative == dbconditionsGood)




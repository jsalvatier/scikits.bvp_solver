'''
Created on Aug 30, 2009

@author: johnsalvatier
'''

from numpy.testing import assert_almost_equal

def assert_solution_matches_data_eval( solution, data, precision = 4):

    yEval, yDerivativeEval = solution(data.xEval, eval_derivative = True)

    # if there is y data in then test it
    if 'yEval' in data.__dict__:
        assert_almost_equal(yEval, data.yEval, precision)

    # if there is derivative data, then test it
    if 'yDerivativeEval' in data.__dict__:
        assert_almost_equal(yDerivativeEval, data.yDerivativeEval, precision)

def assert_solution_matches_data_calculated( solution, data, precision = 4):

    # test the equivalence of the raw stored data
    assert_almost_equal(solution.mesh, data.xSol, precision)
    assert_almost_equal(solution.solution, data.ySol, precision)
    assert_almost_equal(solution.work, data.workSol, precision)
    assert_almost_equal(solution.iwork, data.iworkSol, precision)
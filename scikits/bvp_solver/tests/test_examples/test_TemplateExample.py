'''
Created on Mar 27, 2009

@author: johnsalvatier
'''
from scikits.bvp_solver import get_template
import numpy
from numpy import array

# test to see if the template example runs without errors
def test_templateExample():
    get_template(num_ODE = 3,
                 num_parameters = 1,
                 num_left_boundary_conditions = 1,
                 function_derivative = True,
                 boundary_conditions_derivative = True,
                 singular_term = True)


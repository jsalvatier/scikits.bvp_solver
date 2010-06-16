'''
Created on Mar 27, 2009

@author: johnsalvatier
'''
from scikits.bvp_solver import get_template
import numpy
from numpy import array
import codeop

# test to see if the template example runs without errors
def test_templateExample():
    source1 = get_template(num_ODE = 3,
                 num_parameters = 1,
                 num_left_boundary_conditions = 1,
                 function_derivative = True,
                 boundary_conditions_derivative = True,
                 singular_term = True)
    # test to see if source1 has no syntax errors
    codeop.compile_command(source1)

    source2 = get_template(num_ODE = 3,
                 num_parameters = 0,
                 num_left_boundary_conditions = 1,
                 function_derivative = False,
                 boundary_conditions_derivative = False,
                 singular_term = False)

    # test to see if source2 has no syntax errors
    codeop.compile_command(source2)

    source3 = get_template(num_ODE = 1,
                 num_parameters = 3,
                 num_left_boundary_conditions = 0,
                 function_derivative = False,
                 boundary_conditions_derivative = False,
                 singular_term = False)

    # test to see if source3 has no syntax errors
    codeop.compile_command(source3)

    source4 = get_template(num_ODE = 1,
                 num_parameters = 3,
                 num_left_boundary_conditions = 4,
                 function_derivative = True,
                 boundary_conditions_derivative = True,
                 singular_term = True)

    # test to see if source4 has no syntax errors
    codeop.compile_command(source4)
'''
Created on Mar 27, 2009

@author: johnsalvatier
'''
tab = 6
def spaces (size):
    l = []
    for i in range(size):
        l.append(' ')
    return l
def spacing (size):
    return ''.join(spaces(size))


def array1d( rows,indent, rowComment):
    l = ['array([']
    for i in range(rows):
        if i != 0:

            l.extend(spaces(indent+len('array([') - 2)  )

        l.append('None')
        l.append(listspace)
        if i == rows - 1:
                # if it's the last one, we need to close the array
                l.append('])')
                 #', rowComment ,' ', str(i), '\n'])
        else:
                l.append(',')

        l.extend([' #', rowComment ,' ', str(i), '\n'])

    return l

def list1d(len):
    l = ['[']
    for i in range (len -1):
        l.append('None')
        l.append(listspace)
        l.append(',')

    l.append('None')
    l.append(listspace)
    l.append(']')
    return l



def array2d( rows, col, indent, rowComment):
    l = ['array([']
    for i in range(rows):
        if i != 0:
            l.extend(spaces(indent+ len ('array([') -2))

        l.extend(list1d(col))
        if i == rows - 1:
                # if it's the last one, we need to close the array
                l.append('])')
                 #', rowComment ,' ', str(i), '\n'])
        else:
                l.append(',')

        l.extend([' #', rowComment ,' ', str(i), '\n'])

    return l

def get_template(
                 num_ODE,
                 num_parameters,
                 num_left_boundary_conditions,
                 function_derivative = False,
                 boundary_conditions_derivative = False,
                 singular_term = False,
                 default_object = 'None',
                 array_space = 6
                ):



    """Generates code for solving a boundary value problem using bvp_solver.

    Dummy values (default is None) are used for the places where the user must fill in their own code.

    Suggested use (in ipython or other console)::

        >> print (get_template(parameters))

    and then copy-paste the result into your script.

    Parameters
    ----------
    num_ODE : int
        Number of first order ordinary differential equations in the problem.
    num_parameters : int
        Number of unknown parameters in the problem.
    num_left_boundary_conditions : int
        Number of boundary conditions enforced on the left boundary.
    function_derivative : logical
        Indicates whether a function for the ODE derivatives will be supplied.
    boundary_conditions_derivative : logical
        Indicates whether a function for the boundary conditions derivatives will be supplied.
    singular_term : logical
        Indicates whether a singular term will be supplied.
    default_object : string
        Default string object to appear in arrays. Default is 'None'; '0.0' may also be convienient.
    array_space : int
        Number of spaces between elements in array definitions. More are recommended if terms in the arrays will be long.

    Returns
    -------
    code : string
        Python code skeleton for solving the boundary value problem.
    """

    global listspace
    listspace = spacing(array_space)
    global listobject
    listobject = default_object

    num_right_boundary_conditions = num_ODE + num_parameters - num_left_boundary_conditions

    if num_parameters > 0:
        paramArg = ', P'
    else:
        paramArg = ''

    t = ['import scikits.bvp_solver\n',
         'import numpy\n',
         'from numpy import array\n\n',

    # FUNCTION
         '#all None values in the code below are dummy values which must be replaced with real values\n\n',
         'def function( X, Y', paramArg, '):\n\n',

            '\treturn ']

    t.extend(array1d(num_ODE, tab + len('return '), 'evaluate ODE number'))

    #BOUNDARY CONDITIONS
    t.extend(['\ndef boundary_conditions(Ya, Yb', paramArg, '):\n\n',

                '\tBCa = '])

    t.extend(array1d(num_left_boundary_conditions,tab + len('BCa = '), 'evaluate left BC number'))

    t.append('\n\tBCb = ')
    t.extend(array1d(num_right_boundary_conditions,tab + len('BCb = '), 'evaluate right BC number'))
    t.append('\treturn BCa, BCb')

    t.append('\n\n')

    #FUNCTION_DERIVATIVE
    if function_derivative == True:
         t.extend(['def function_derivative( X, Y', paramArg, '):\n\n',
            '\t#evaluate function derivative with respect to variables',
            '\t#increasing differentiation index\n',
            '\t#(Yi) ---->\n',
            '\tdFdY = '])
         t.extend(array2d(num_ODE,num_ODE,tab +  len('dFdY = '), 'ODE number'))

         if num_parameters > 0:
             t.extend(['\n\t#increasing differentiation index\n',
                       '\t# (Pi) ---->\n',
                       '\tdFdP = '])
             t.extend(array2d(num_ODE,num_parameters,tab +  len('dFdP = '), 'ODE number'))

             t.append('\n\treturn dFdY, dFdP\n')
         else:
             t.append('\n\treturn dFdY\n')

    #BOUNDARY_CONDITIONS_DERIVATIVE
    if boundary_conditions_derivative == True:
         t.extend(['\ndef boundary_conditions_derivative( Ya, Yb', paramArg, '):\n\n',
                   '\t#evaluate left boundary conditions derivative with respect to variables\n',
                   '\t#increasing differentiation index\n',
                   '\t# (Ya,i) ----->\n',
                   '\tdBCadYa = '])
         t.extend(array2d(num_left_boundary_conditions,num_ODE,tab +  len('dBCadYa = '), 'left BC number'))

         if num_parameters > 0:
             t.extend(['\n\t#evaluate left boundary conditions derivative with respect to unknown parameters\n',
                       '\t#increasing differentiation index\n',
                       '\t# (Pi) ----->\n'])
             t.append('\tdBCadP = ')
             t.extend(array2d(num_left_boundary_conditions,num_parameters,tab +  len('dBCadP = '), 'left BC number'))

         t.extend(['\n\t#evaluate right boundary conditions derivative with respect to variables\n',
                   '\t#increasing differentiation index\n',
                   '\t# (Yb,i) ----->\n'])
         t.append('\tdBCbdYb = ')
         t.extend(array2d(num_right_boundary_conditions,num_ODE,tab +  len('dBCbdYb = '), 'right BC number'))

         if num_parameters > 0:
             t.extend(['\n\t#evaluate right boundary conditions derivative with respect to unknown parameters\n',
                       '\t#increasing differentiation index\n',
                       '\t# (Pi) ----->\n'])
             t.append('\tdBCbdP = ')
             t.extend(array2d(num_right_boundary_conditions,num_parameters,tab +  len('dBCbdP = '), 'right BC number'))

         if num_parameters > 0:
             t.append('\n\treturn dBCadYa, dBCbdYb, dBCadP, dBCbdP\n')
         else:
             t.append('\n\treturn dBCadYa, dBCbdYb\n')

    s = len('problem_definition = bvp_solver.ProblemDefinition(')
    if function_derivative == True:
        fderive_param = ',\n' + spacing(s) + 'function_derivative = function_derivative'
    else:
        fderive_param = ''

    if boundary_conditions_derivative == True:
        bcderive_param = ',\n' + spacing(s) + 'boundary_conditions_derivative = boundary_conditions_derivative'
    else:
        bcderive_param = ''


    t.extend(['\n\nproblem_definition = scikits.bvp_solver.ProblemDefinition(num_ODE = ', str(num_ODE), ',\n',
                                                         spacing(s),'num_parameters = ', str(num_parameters), ',\n',
                                                         spacing(s),'num_left_boundary_conditions = ', str(num_left_boundary_conditions), ',\n',
                                                         spacing(s),'boundary_points = (None, None),\n',
                                                         spacing(s),'function = function,\n',
                                                         spacing(s),'boundary_conditions = boundary_conditions',
                                                         fderive_param,
                                                         bcderive_param, ')'])

    s = len('solution = scikits.bvp_solver.solve(')

    if singular_term == True:
        t.append('\n\n\nsingular_term = ')
        t.extend(array2d(num_ODE, num_ODE, len('singular_term = '), ''))

        singular_term = ',\n' + spacing(s) + 'singular_term = singular_term'
    else:
        singular_term = ''


    if num_parameters > 0:
        paramguessString = list1d(num_parameters)
        parameter_guess = ',\n' + spacing(s) + 'parameter_guess = ' + ''.join(paramguessString)
    else:
        parameter_guess = ''


    yguessString = list1d(num_ODE)
    yguessString = ''.join(yguessString)


    t.extend(['\nsolution = scikits.bvp_solver.solve(bvp_problem = problem_definition,\n',
                                   spacing(s),'solution_guess = ', yguessString,
                                               parameter_guess,
                                               singular_term,')'])
    return ''.join(t)
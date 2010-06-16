import scikits.bvp_solver
import numpy
from numpy import array

#all None values in the code below are dummy values which must be replaced with real values

def function( X, Y, P):

    return array([None      , #evaluate ODE number 0
                  None      , #evaluate ODE number 1
                  None      ]) #evaluate ODE number 2

def boundary_conditions(Ya, Yb, P):

    BCa = array([None      ]) #evaluate left BC number 0

    BCb = array([None      , #evaluate right BC number 0
                 None      , #evaluate right BC number 1
                 None      ]) #evaluate right BC number 2
    return BCa, BCb

def function_derivative( X, Y, P):

    #evaluate function derivative with respect to variables    #increasing differentiation index
    #(Yi) ---->
    dFdY = array([[None      ,None      ,None      ], #ODE number 0
                  [None      ,None      ,None      ], #ODE number 1
                  [None      ,None      ,None      ]]) #ODE number 2

    #increasing differentiation index
    # (Pi) ---->
    dFdP = array([[None      ], #ODE number 0
                  [None      ], #ODE number 1
                  [None      ]]) #ODE number 2

    return dFdY, dFdP

def boundary_conditions_derivative( Ya, Yb, P):

    #evaluate left boundary conditions derivative with respect to variables
    #increasing differentiation index
    # (Ya,i) ----->
    dBCadYa = array([[None      ,None      ,None      ]]) #left BC number 0

    #evaluate left boundary conditions derivative with respect to unknown parameters
    #increasing differentiation index
    # (Pi) ----->
    dBCadP = array([[None      ]]) #left BC number 0

    #evaluate right boundary conditions derivative with respect to variables
    #increasing differentiation index
    # (Yb,i) ----->
    dBCbdYb = array([[None      ,None      ,None      ], #right BC number 0
                     [None      ,None      ,None      ], #right BC number 1
                     [None      ,None      ,None      ]]) #right BC number 2

    #evaluate right boundary conditions derivative with respect to unknown parameters
    #increasing differentiation index
    # (Pi) ----->
    dBCbdP = array([[None      ], #right BC number 0
                    [None      ], #right BC number 1
                    [None      ]]) #right BC number 2

    return dBCadYa, dBCbdYb, dBCadP, dBCbdP


problem_definition = scikits.bvp_solver.ProblemDefinition(num_ODE = 3,
                                              num_parameters = 1,
                                              num_left_boundary_conditions = 1,
                                              boundary_points = (None, None),
                                              function = function,
                                              boundary_conditions = boundary_conditions,
                                              function_derivative = function_derivative,
                                              boundary_conditions_derivative = boundary_conditions_derivative)


singular_term = array([[None      ,None      ,None      ], # 0
                     [None      ,None      ,None      ], # 1
                     [None      ,None      ,None      ]]) # 2

solution = scikits.bvp_solver.solve(bvp_problem = problem_definition,
                                    solution_guess = [None      ,None      ,None      ],
                                    parameter_guess = [None      ],
                                    singular_term = singular_term)
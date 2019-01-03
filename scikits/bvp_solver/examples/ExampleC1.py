import numpy as np 
import matplotlib.pyplot as plt
import scikits.bvp_solver as bvp
#import ComplexProblemDefinition as cbvp
import sys

def function(X,Y,P):

    k = P[0]
    return np.array([ 1j*k*Y[1],
                      1j*k*Y[0]])

def function_derivative(X,Y,P):
    dFdY = np.array([[     0.0 , 1j*P[0] ],
                        [ 1j*P[0] ,     0.0 ]])
    dFdP = np.array([[ 1j*Y[1]],
                     [ 1j*Y[0]]])

    return dFdY, dFdP

def boundary_conditions(Ya,Yb,P):
    BCa = np.array([Ya[0],
                    Ya[1] - 1.0])

    BCb = np.array([Yb[0]])
    return BCa, BCb

def boundary_conditions_derivative(Ya, Yb, P):
    dBCa = np.array([[1.0, 0.0],
                     [0.0, 1.0]])

    dBCb = np.array([[1.0, 0.0]])

    dBPa = np.zeros((2,1), dtype=complex)
    dBPb = np.zeros((1,1), dtype=complex)
    return dBCa, dBCb, dBPa, dBPb

if __name__=='__main__':

    P0 = 5.0
    if len(sys.argv)==2:
        P0 = float(sys.argv[1])

    problem = bvp.ComplexProblemDefinition(num_ODE_c = 2,
                                    num_parameters_c = 1,
                                    num_left_boundary_conditions_c = 2,
                                    boundary_points = (0, 2.0*np.pi),
                                    function_c = function,
                                    boundary_conditions_c = boundary_conditions,
                                    function_derivative_c = function_derivative,
                                    boundary_conditions_derivative_c = boundary_conditions_derivative
                                    )

    def guess(X, P0):
        return np.array( [np.sin(P0*X), 0.0, 0.0, -np.cos(P0*X)] )

    initmesh = np.linspace(problem.boundary_points[0],problem.boundary_points[1],512)
    solution = bvp.solve(problem,
                        initial_mesh = initmesh,
                        solution_guess = lambda x: guess(x, P0),
                        parameter_guess = np.array([P0, 0.0]),
                        trace = 2)

    x = np.linspace(problem.boundary_points[0],problem.boundary_points[1], 100)
    y = solution(x)
    print("lambda = " + str(solution.parameters[0]))

    solution.save("test_ExampleC1.sol")


    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x, y[0,:], c='C0')
    ax.plot(x, y[1,:], c='C0', ls='--')
    ax.plot(x, y[2,:], c='C1')
    ax.plot(x, y[3,:], c='C1', ls='--')
    plt.show()
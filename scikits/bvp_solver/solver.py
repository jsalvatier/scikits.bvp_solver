import bvp_solverf
import numpy
from Solution import Solution
import tools

def solve(bvp_problem,
          solution_guess,
          initial_mesh = None,
          parameter_guess = None,
          max_subintervals = 300,
          singular_term = None,
          tolerance = 1.0e-6,
          method = 4,
          trace = 0,
          error_on_fail = True):
    """Attempts to solve the supplied boundary value problem starting from the user supplied guess for the solution using BVP_SOLVER.

    Parameters
    ----------
    bvp_problem : :class:`ProblemDefinition`
        Defines the boundary value problem to be solved.
    solution_guess : :class:`Solution`, constant, array of values or function
        A guess for the solution.
    initial_mesh : castable to floating point ndarray
        Points on the x-axis to use for the supplied solution guess, default is 10 evenly spaced points. Must not be supplied if solution_guess is a :class:`Solution` object.
    parameter_guess : castable to floating point ndarray, shape (num_parameters)
        Guesses for the unknown parameters. Must not be supplied if solution_guess is a :class:`Solution` object.
    max_subintervals : int
        Maximum number of points on the mesh before an error is returned.
    singular_term : castable to floating point ndarray, shape(num_ODE, num_ODE)
        Matrix that defines the singular term for the problem if one exist.
    tolerance : positive float
        Tolerance for size of defect of approximate solution.
    method : {2, 4, 6}
        Order of Runge-Kutta to use.
    trace : {0, 1, 2}
        Indicates verbosity of output. 0 for no output, 1 for some output, 2 for full output.
    error_on_fail : logical
        Indicates whether an exception should be raised if solving fails.

    Returns
    -------
    sol : :class:`Solution`
        Approximate problem solution.

    Raises
    ------
    ValueError
        If bvp_problem failed validation in some way.
    ValueError
        If solving fails.
    """

    init_solution = 0

    if isinstance(solution_guess, Solution):
        if not (initial_mesh == None and
                parameter_guess == None):
            raise ValueError("values for initial mesh and parameter_guess must not be given if solution_guess is a Solution object")
        init_solution = solution_guess
    else:

        if initial_mesh == None:
            initial_mesh = numpy.linspace(bvp_problem.boundary_points[0],bvp_problem.boundary_points[1] , 10)

        # here we call one of the BVP_GUESS_i routines that make up BVP_INIT to set up the solution
        if ( not callable(solution_guess)):   # in this case the initial solution passed was not a function

            #try to cast the solution to an array
            solution_guess = numpy.array(solution_guess)

            # if the solution_guess is just an array the size of ODE
            if solution_guess.shape == (bvp_problem.num_ODE,) or (solution_guess.shape == () and bvp_problem.num_ODE == 1):

                bvp_solverf.bvp.guess_1_wrap(nparam_in = bvp_problem.num_parameters,
                                           leftbc_in = bvp_problem.num_left_boundary_conditions,
                                           x_in = tools.farg(initial_mesh),
                                           y_in = tools.farg(solution_guess),
                                           parameters_in = tools.farg(parameter_guess),
                                           mxnsub_in = max_subintervals,
                                           node_in = bvp_problem.num_ODE)

                init_solution = Solution.from_arg_list(bvp_solverf.bvp)

            else:

                tools.argShapeTest(solution_guess, (bvp_problem.num_ODE,initial_mesh.shape[0]), "solution guess")

                bvp_solverf.bvp.guess_2_wrap(nparam_in = bvp_problem.num_parameters,
                                             leftbc_in = bvp_problem.num_left_boundary_conditions,
                                             x_in = tools.farg(initial_mesh),
                                             y_in = tools.farg(solution_guess),
                                             parameters_in = tools.farg(parameter_guess),
                                             mxnsub_in = max_subintervals,
                                             node_in = bvp_problem.num_ODE)
                init_solution = Solution.from_arg_list(bvp_solverf.bvp)

        else:
            y_in = numpy.zeros((bvp_problem.num_ODE,1))

            bvp_solverf.bvp.guess_1_wrap(nparam_in = bvp_problem.num_parameters,
                                           leftbc_in = bvp_problem.num_left_boundary_conditions,
                                           x_in = tools.farg(initial_mesh),
                                           y_in = y_in,
                                           parameters_in = tools.farg(parameter_guess),
                                           mxnsub_in = max_subintervals,
                                           node_in = bvp_problem.num_ODE)

            init_solution = Solution.from_arg_list(bvp_solverf.bvp)
            for i,v in enumerate(init_solution._mesh):
                init_solution._solution[:, i] = solution_guess(v)


    if not (method == 2 or method == 4 or method == 6 ):
        raise ValueError ("method must be either 2, 4 or 6 but got " + str(method) )

    if (tolerance < 0):
        raise ValueError("tolerance must be nonnegative")

    singular = not (singular_term is None)

    # check to see if the singular term is of the right size
    singular_term = tools.preparg(singular_term)
    if singular and not (singular_term.shape == (bvp_problem.num_ODE, bvp_problem.num_ODE)):
        raise ValueError("singular_term has the wrong shape/size. Expected: " +
                         (bvp_problem.num_ODE, bvp_problem.num_ODE)+
                         " but got :" +
                          singular_term.shape)

    # test the problem specifications with the initial solution
    bvp_problem.test(init_solution)
                                    
                                
                    
    bvp_solverf.bvp.bvp_solver_wrap(node_in = bvp_problem.num_ODE,
                                    npar_in = bvp_problem.num_parameters,
                                    leftbc_in = bvp_problem.num_left_boundary_conditions,
                                    npts_in = len(init_solution.mesh),
                                    info_in = init_solution.successIndicator,
                                    mxnsub_in = max_subintervals,
                                    x_in = tools.farg(init_solution.mesh),
                                    y_in = tools.farg(init_solution.solution),
                                    parameters_in = tools.farg(init_solution.parameters),
                                    work_in = tools.farg(init_solution.work),
                                    iwork_in = tools.farg(init_solution.iwork),
                                    fsub = bvp_problem._function,
                                    fsubp = bvp_problem._functionp,
                                    bcsub = bvp_problem._boundary_conditions,
                                    bcsubp = bvp_problem._boundary_conditionsp,

                                    singular = singular,

                                    hasdfdy = bvp_problem.has_function_derivative,
                                    dfdy = bvp_problem._function_derivative,
                                    dfdyp = bvp_problem._function_derivativep,

                                    hasdbcdy = bvp_problem.has_boundary_conditions_derivative,
                                    dbcdy = bvp_problem._boundary_conditions_derivative,
                                    dbcdyp = bvp_problem._boundary_conditions_derivativep,
                                    #optional arguments
                                    method = method,
                                    tol = tolerance,

                                    trace = trace,
                                    # we never want the actual program to shut down from the fortran side, we should throw an error
                                    stop_on_fail = False,

                                    singularterm = tools.farg(singular_term))

    calculatedSolution = Solution.from_arg_list(bvp_solverf.bvp)

    # check to see if there was a problem with the solution
    if error_on_fail and calculatedSolution._successIndicator == -1:
        raise ValueError("Boundary value problem solving failed. Run with trace = 1 or 2 for more information.")

    return calculatedSolution

def _guess_3_wrap(node_in,
                 nparam_in ,
                 leftbc_in,
                 x_in,
                 fcn,
                 parameters_in,
                 mxnsub_in):
    y_in = 0.0
    sol = bvp_solverf.bvp.guess_1_wrap(node_in,nparam_in, leftbc_in, x_in, y_in,
            parameters_in, mxnsub_in)
    

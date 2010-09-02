import numpy
import tools
import pickle

class ProblemDefinition:
    """Defines a boundary value problem.
    """
    def __init__(self,
                 num_ODE,
                 num_parameters,
                 num_left_boundary_conditions,
                 boundary_points,
                 function,
                 boundary_conditions,
                 function_derivative = None,
                 boundary_conditions_derivative = None):
        """
        Parameters
        ----------
        num_ODE : int
            Number of first order ordinary differential equations in the problem.
        num_parameters : int
            Number of unknown parameters in the problem.
        num_left_boundary_conditions : int
            Number of boundary conditions enforced on the left boundary.
        boundary_points : arraylike, shape(2)
            Array that defines the two boundary points on the x axis.
        function : function (see definition below)
            A function which calculates the value of the ODE equations.

                function (X, Y[, P]):

                    Parameters:
                        X : float
                            scalar value of x at which to evaluate the ODEs

                        Y : ndarray, shape(num_ODE)
                            current value of all variables

                        P : ndarray, shape(num_parameters)
                            value of all unknown parameters (only included if num_parameters > 0)
                    Returns:
                        ODE : ndarray, shape(num_ODE)
                            array of all ODEs

        boundary_conditions : function (see definition below)
            A function which calculates the difference between the boundary conditions and the actual variables currently calculated.

                boundary_conditions(YA, YB[, P]):

                    Parameters:
                        YA : ndarray, shape(num_ODE)
                            value of all variables at the left boundary

                        YB : ndarray, shape(num_ODE)
                            value of all variables at the right boundary

                        P : ndarray, shape(num_parameters)
                            value of all unknown parameters (only used if num_parameters > 0)

                    Returns:
                        BCA : ndarray, shape(num_left_boundary_conditions)
                            difference between the boundary condition and variables at the left boundary

                        BCB : ndarray, shape(num_ODE + num_parameters - num_left_boundary_conditions)
                            array of the difference between the boundary condition and variables at the right boundary

        function_derivative : optional function (see definition below)
            A function which returns the partial derivatives of the function argument.

                function (X, Y[, P]):

                    Parameters:
                        X : float
                            scalar value of x at which to evaluate the ODEs

                        Y : ndarray, shape(num_ODE)
                            current value of all variables

                        P : ndarray, shape(num_parameters)
                            value of all unknown parameters (only included if num_parameters > 0)
                    Returns:
                        dODE : ndarray, shape(num_ODE, num_ODE)
                            array of partial derivative of all ODEs with respect to all variables; index of ODEs is first, index of variables is second

                        dOdP : ndarray, shape(num_ODE, num_parameters)
                            array of partial derivative of all ODEs with respect to all unknown parameters; index of ODEs is first, index of parameters is second must not be returned if the problem does not include unknown parameters


        boundary_conditions_drivative : optional function (see definition below)
            A function which returns the partial derivatives of the boundary_conditions argument.
                boundary_conditions(YA, YB[, P]):

                    Parameters:
                        YA : ndarray, shape(num_ODE)
                            value of all variables at the left boundary

                        YB : ndarray, shape(num_ODE)
                            value of all variables at the right boundary

                        P : ndarray, shape(num_parameters)
                            value of all unknown parameters (only used if num_parameters > 0)

                    Returns:
                        dBCA : ndarray, shape(num_left_boundary_conditions, num_ODE)
                            partial derivatives of the difference between the left boundary condition and the actual variables at the left boundary; boundary condition index is first and variable index is second

                        dBCB : ndarray, shape(num_ODE + num_parameters - num_left_boundary_conditions, num_ODE)
                            partial derivatives of the difference between the right boundary condition and the actual variables at the right boundary; boundary condition index is first and variable index is second

                        dBPA : ndarray, shape(num_left_boundary_conditions, num_parameters)
                            partial derivatives of the difference between the left boundary condition and the unknown parameters; boundary condition index is first and parameter index is second

                        dBPB : ndarray, shape(num_ODE + num_parameters - num_left_boundary_conditions, num_parameters)
                            partial derivatives of the difference between the right boundary condition and the unknown parameters; boundary condition index is first and parameter index is second

        """

        self._num_ODE = num_ODE
        self._num_parameters = num_parameters
        self._num_left_boundary_conditions = num_left_boundary_conditions
        self._boundary_points = tools.preparg(boundary_points)
        self._function = function
        self._function_store = function

        self._boundary_conditions =boundary_conditions
        self._boundary_conditions_store =boundary_conditions

        self._function_derivative_store = function_derivative
        self._function_derivative = function_derivative

        self._boundary_conditions_derivative_store = boundary_conditions_derivative
        self._boundary_conditions_derivative = boundary_conditions_derivative


        #figure out whether the user has supplied derivatives for the function and boundary conditions
        self.has_function_derivative = not self._function_derivative is None
        self.has_boundary_conditions_derivative = not self._boundary_conditions_derivative is None

        if self._num_parameters == 0:

            # if don't have unknown parameters then give all the arguments for callbacks with parameter arguments dummy functions
            self._functionp = fp_dummy
            self._boundary_conditionsp = bcp_dummy
            self._function_derivativep = fderivep_dummy
            self._boundary_conditions_derivativep = bcderivep_dummy

            # also assign dummy arguments to optional derivatives if they were not supplied
            if self._function_derivative is None:
                self._function_derivative = f_dummy

            if self._boundary_conditions_derivative is None:
                self._boundary_conditions_derivative = bcderive_dummy
        else:

            #also assign all the supplied arguments to the callback arguments with parameter arguments in them
            self._functionp = function
            self._boundary_conditionsp = self._boundary_conditions
            self._function_derivativep = self._function_derivative
            self._boundary_conditions_derivativep = self._boundary_conditions_derivative

            # if have unknown parameters then give all the arguments for callbacks without parameter arguments dummy functions
            self._function = f_dummy
            self._boundary_conditions = bc_dummy
            self._function_derivative = fderive_dummy
            self._boundary_conditions_derivative = bcderive_dummy

            # give dummy arguments in those cases when no derivative arugment is supplied
            if self._function_derivativep is None:
                self._function_derivativep = fp_dummy


            if self._boundary_conditions_derivativep is None:
                self._boundary_conditions_derivativep = bcderivep_dummy

    @property
    def num_ODE(self):
        """
        read only, int
            Number of first order ordinary differential equations in the problem.
        """
        return self._num_ODE

    @property
    def num_parameters(self):
        """
        read only, int
            Number of unknown parameters in the problem.
        """
        return self._num_parameters

    @property
    def num_left_boundary_conditions(self):
        """
        read only, int
            Number of boundary conditions enforced on the left boundary.
        """
        return self._num_left_boundary_conditions

    @property
    def boundary_points(self):
        """
        read only, ndarray, shape(2)
            Array that containing the location of the two boundary points on the x axis.
        """
        return self._boundary_points

    @property
    def function(self):
        """
        read only, function
            Function which calculates the value of the ODE equations.
        """
        return self._function_store

    @property
    def boundary_conditions(self):
        """
        read only, function
            Function which calculates the difference between the actual boundary conditions and the desired boundary conditions.
        """
        return self._boundary_conditions_store

    @property
    def function_derivative(self):
        """
        read only, function
            A function which returns the partial derivatives of the ODEs.
        """
        return self._function_derivative_store

    @property
    def boundary_conditions_derivative(self):
        """
        read only, function
            A function which returns the partial derivatives of the boundary conditions.
        """
        return self._boundary_conditions_derivative_store


    def test(self, test_solution):
        """Test that the boundary value problem definition is self consistent, and
        tests whether test_solution is consistent with the bvp definition.
        This requires some legal values for the parameters and thus requires a test solution.
        
        Parameters
        ----------
        test_solution : :class:`Solution`
            solution to be tested with
        """

        if self._num_parameters > 0:
            tools.argShapeTest(test_solution.parameters, (self._num_parameters,),
                               "parameter array",
                               "Should be (num_parameters,)")

        tools.argShapeTest(test_solution.solution, (self._num_ODE,len(test_solution.mesh)),
                               "solution array",
                               "Should be (num_ODE,length(mesh))")

        if not (self._boundary_points.shape == (2,)):
            raise ValueError("This boundary value problem definition must be given exactly two boundary points, but got: "
                              + self._boundary_points +" as the boundary values")

        #at this point we want to check the call backs to make sure they take the right arguments and return the right things
        if self._num_parameters ==0:
            f = self._function(test_solution.mesh[0],
                         test_solution.solution[:,0])

            tools.argShapeTest(f, (self._num_ODE,),
                               "function callback return",
                               "Should be (num_ODE,)")

            bca,bcb = self._boundary_conditions(test_solution.solution[:,0],
                                     test_solution.solution[:,-1])

            tools.argShapeTest(bca, (self._num_left_boundary_conditions,),
                               "Boundary conditions callback first return",
                               "Should be (num_left_boundary_conditions,)")

            tools.argShapeTest(bcb, (self._num_ODE + self._num_parameters - self._num_left_boundary_conditions,),
                               "Boundary conditions callback second return",
                               "Should be (num_ODE + num_parameters - num_left_boundary_conditions,)")

            if self.has_function_derivative:
                df = self._function_derivative(test_solution.mesh[0],
                                         test_solution.solution[:,0])

                tools.argShapeTest(df, (self._num_ODE, self._num_ODE),
                                   "function derivative callback first return",
                                   "Should be (num_ODE, num_ODE)")

            if self.has_boundary_conditions_derivative:
                dbca, dbcb = self._boundary_conditions_derivative(test_solution.solution[:,0],
                                                            test_solution.solution[:,0])

                tools.argShapeTest(dbca, (self._num_left_boundary_conditions, self._num_ODE),
                                   "Boundary conditions derivative callback first return",
                                   "Should be (num_left_boundary_conditions, num_ODE)")

                tools.argShapeTest(dbcb, (self._num_ODE + self._num_parameters - self._num_left_boundary_conditions, self._num_ODE),
                                   "Boundary conditions derivative callback second return",
                                   "Should be (num_ODE + num_parameters - num_left_boundary_conditions, num_ODE)")

        else: ## if unknown parameters are used, things should be a little different
            f = self._functionp(test_solution.mesh[0],
                         test_solution.solution[:,0],
                         test_solution.parameters)

            tools.argShapeTest(f, (self._num_ODE,),
                               "function callback return",
                               "Should be (num_ODE,)")

            bca, bcb = self._boundary_conditionsp(test_solution.solution[:,0],
                                     test_solution.solution[:,0],
                                     test_solution.parameters)

            tools.argShapeTest(bca, (self._num_left_boundary_conditions,),
                               "Boundary conditions callback first return",
                               "Should be (num_left_boundary_conditions,)")

            tools.argShapeTest(bcb, (self._num_ODE + self._num_parameters - self._num_left_boundary_conditions,),
                               "Boundary conditions callback second return",
                               "Should be (num_ODE + num_parameters - num_left_boundary_conditions,)")

            if self.has_function_derivative:
                df, dfp = self._function_derivativep(test_solution.mesh[0],
                                         test_solution.solution[:,0],
                                         test_solution.parameters)

                tools.argShapeTest(df, (self._num_ODE, self._num_ODE),
                                   "function derivative callback first return",
                                   "Should be (num_ODE, num_ODE)")

                tools.argShapeTest(dfp, (self._num_ODE, self._num_parameters),
                                   "function derivative callback second return",
                                   "Should be (num_ODE ,num_parameters)")

            if self.has_boundary_conditions_derivative:
                dbca, dbcb, dbcap, dbcbp = self._boundary_conditions_derivativep(test_solution.solution[:,0],
                                                            test_solution.solution[:,0],
                                                            test_solution.parameters)

                tools.argShapeTest(dbca, (self._num_left_boundary_conditions, self._num_ODE),
                             "boundary conditions derivative callback first return",
                             "Should be (num_left_boundary_conditions, num_ODE)")

                tools.argShapeTest(dbcb, (self._num_ODE + self._num_parameters - self._num_left_boundary_conditions, self._num_ODE),
                             "boundary conditions derivative callback second return",
                             "Should be (num_ODE + num_parameters - num_left_boundary_conditions, num_ODE)")

                tools.argShapeTest(dbcap, (self._num_left_boundary_conditions, self._num_parameters),
                             "boundary conditions derivative callback third return",
                             "Should be (num_left_boundary_conditions, num_parameters)")

                tools.argShapeTest(dbcbp, (self._num_ODE + self._num_parameters - self._num_left_boundary_conditions, self._num_parameters),
                             "boundary conditions derivative callback fourth return",
                             "Should be (num_ODE + num_parameters - num_left_boundary_conditions, num_ODE, num_parameters)")

        # test derivatives
        step = 1e-8
        places = 4

        # chose the point near the middle of the test_solution to check the derivatives
        middlePoint = numpy.round(test_solution.mesh.size * .61, 0)


        if self.has_function_derivative:
            if self._num_parameters > 0:

                func_derivative_calc, func_param_derivative_calc = self._function_derivativep(test_solution.mesh[middlePoint], test_solution.solution[:, middlePoint], test_solution.parameters)
                func_derivative_num = numpy.zeros((self._num_ODE,self._num_ODE))
                func_param_derivative_num = numpy.zeros( (self._num_ODE, self._num_parameters))



                for i in range(self._num_ODE):
                    delta = numpy.zeros(self._num_ODE)
                    delta[i] += step

                    point1 = self._functionp(test_solution.mesh[middlePoint], test_solution.solution[:, middlePoint], test_solution.parameters)
                    point2 = self._functionp(test_solution.mesh[middlePoint], test_solution.solution[:, middlePoint] + delta, test_solution.parameters)

                    func_derivative_num[:, i] = (point2 - point1)/step

                for i in range(self._num_parameters):
                    delta = numpy.zeros(self._num_parameters)
                    delta[i] += step

                    point1 = self._functionp(test_solution.mesh[middlePoint], test_solution.solution[:, middlePoint], test_solution.parameters)
                    point2 = self._functionp(test_solution.mesh[middlePoint], test_solution.solution[:, middlePoint], test_solution.parameters + delta)

                    func_param_derivative_num[:, i] = (point2 - point1)/step

                # now compare the actual derivatives with the calculated ones

                difference = func_derivative_calc - func_derivative_num
                if not (numpy.round(difference, places) == 0).all():
                    raise ValueError("analytical derivative matrix does not match numerical derivative matrix.\n Analytical is:\n" + str(func_derivative_calc)
                                     + "\n Numerical is:\n" + str(func_derivative_num))

                parameterDifference = func_param_derivative_calc - func_param_derivative_num
                if not (numpy.round(parameterDifference, places) == 0).all():
                    raise ValueError("analytical derivative (with respect to parameters) matrix does not match numerical derivative matrix.\n Analytical is:\n" + str(func_param_derivative_calc)
                                     + "\n Numerical is:\n" + str(func_param_derivative_num))
            else:

                func_derivative_calc = self._function_derivative(test_solution.mesh[middlePoint], test_solution.solution[:, middlePoint])
                func_derivative_num = numpy.zeros((self._num_ODE,self._num_ODE))

                for i in range(self._num_ODE):
                    delta = numpy.zeros(self._num_ODE)
                    delta[i] += step

                    point1 = self._function(test_solution.mesh[middlePoint], test_solution.solution[:, middlePoint])
                    point2 = self._function(test_solution.mesh[middlePoint], test_solution.solution[:, middlePoint] + delta)

                    func_derivative_num[:, i] = (point2 - point1)/step

                # now compare the actual derivatives with the calculated ones
                difference = func_derivative_calc - func_derivative_num
                if not (numpy.round(difference, places) == 0).all():
                    raise ValueError("analytical derivative matrix does not match numerical derivative matrix.\n Analytical is:\n" + str(func_derivative_calc)
                                     + "\n Numerical is:\n" + str(func_derivative_num))

            if self.has_boundary_conditions_derivative:

                endPoint = test_solution.mesh.size - 1

                if self._num_parameters > 0:

                    bcA_derivative_calc,bcB_derivative_calc, bcA_param_derivative_calc, bcB_param_derivative_calc = self._boundary_conditions_derivativep(test_solution.solution[:, 0],
                                                                                                                                                        test_solution.solution[:,endPoint],
                                                                                                                                                        test_solution.parameters)

                    bcA_derivative_num = numpy.zeros((self._num_left_boundary_conditions, self._num_ODE))
                    bcB_derivative_num = numpy.zeros((self._num_ODE + self._num_parameters - self._num_left_boundary_conditions, self._num_ODE))
                    bcA_param_derivative_num = numpy.zeros((self._num_left_boundary_conditions, self._num_parameters))
                    bcB_param_derivative_num = numpy.zeros( (self._num_ODE + self._num_parameters - self._num_left_boundary_conditions, self._num_parameters))

                    # calculate analytic derivatives for
                    for i in range(self._num_ODE):
                        delta = numpy.zeros(self._num_ODE)
                        delta[i] += step

                        point1A, point1B = self._boundary_conditionsp(test_solution.solution[:, 0],test_solution.solution[:, endPoint], test_solution.parameters)
                        point2A, point2B = self._boundary_conditionsp(test_solution.solution[:, 0] + delta ,test_solution.solution[:, endPoint] + delta, test_solution.parameters)

                        bcA_derivative_num[:, i] = (point2A - point1A)/step
                        bcB_derivative_num[:, i] = (point2B - point1B)/step

                    for i in range(self._num_parameters):
                        delta = numpy.zeros(self._num_parameters)
                        delta[i] += step

                        point1A, point1B = self._boundary_conditionsp(test_solution.solution[:, 0],test_solution.solution[:, endPoint], test_solution.parameters)
                        point2A, point2B = self._boundary_conditionsp(test_solution.solution[:, 0] ,test_solution.solution[:, endPoint], test_solution.parameters + delta)

                        bcA_param_derivative_num[:, i] = (point2A - point1A)/step
                        bcB_param_derivative_num[:, i] = (point2B - point1B)/step

                    # now compare the actual derivatives with the calculated ones

                    if not  (numpy.round(bcA_derivative_calc - bcA_derivative_num, places) == 0).all():
                        raise ValueError("analytical derivative matrix for the left boundary condition does not match numerical derivative matrix.\n Analytical is:\n" + str(bcA_derivative_calc)
                                     + "\n Numerical is:\n" + str(bcA_derivative_num))

                    if not  (numpy.round(bcB_derivative_calc - bcB_derivative_num, places) == 0).all():
                        raise ValueError("analytical derivative matrix for the right boundary condition does not match numerical derivative matrix.\n Analytical is:\n" + str(bcB_derivative_calc)
                                     + "\n Numerical is:\n" + str(bcB_derivative_num))

                    if not  (numpy.round(bcA_param_derivative_calc - bcA_param_derivative_num, places) == 0).all():
                        raise ValueError("analytical derivative matrix (with respect to parameters) for the left boundary condition does not match numerical derivative matrix.\n Analytical is:\n" + str(bcA_param_derivative_calc)
                                     + "\n Numerical is:\n" + str(bcA_param_derivative_num))

                    if not  (numpy.round(bcB_param_derivative_calc - bcB_param_derivative_num, places) == 0).all():
                        raise ValueError("analytical derivative matrix (with respect to parameters) for the left boundary condition does not match numerical derivative matrix.\n Analytical is:\n" + str(bcB_param_derivative_calc)
                                     + "\n Numerical is:\n" + str(bcB_param_derivative_num))

                else:
                    bcA_derivative_calc,bcB_derivative_calc = self._boundary_conditions_derivative(test_solution.solution[:, 0],
                                                                                                   test_solution.solution[:,endPoint])

                    bcA_derivative_num = numpy.zeros((self._num_left_boundary_conditions, self._num_ODE))
                    bcB_derivative_num = numpy.zeros((self._num_ODE - self._num_left_boundary_conditions, self._num_ODE))

                    # calculate analytic derivatives for
                    for i in range(self._num_ODE):
                        delta = numpy.zeros(self._num_ODE)
                        delta[i] += step

                        point1A, point1B = self._boundary_conditions(test_solution.solution[:, 0],test_solution.solution[:, endPoint])
                        point2A, point2B = self._boundary_conditions(test_solution.solution[:, 0] + delta ,test_solution.solution[:, endPoint] + delta)

                        bcA_derivative_num[:, i] = (point2A - point1A)/step
                        bcB_derivative_num[:, i] = (point2B - point1B)/step

                    # now compare the actual derivatives with the calculated ones


                    if not  (numpy.round(bcA_derivative_calc - bcA_derivative_num, places) == 0).all():
                        raise ValueError("analytical derivative matrix for the left boundary condition does not match numerical derivative matrix.\n Analytical is:\n" + str(bcA_derivative_calc)
                                     + "\nNumerical is:\n" + str(bcA_derivative_num))

                    if not  (numpy.round(bcB_derivative_calc - bcB_derivative_num, places) == 0).all():
                        raise ValueError("analytical derivative matrix for the right boundary condition does not match numerical derivative matrix.\n Analytical is:\n" + str(bcB_derivative_calc)
                                     + "\nNumerical is:\n" + str(bcB_derivative_num))


# these are dummy functions for when either one of the derivative functions is not defined
# and to replace the fuctions that either do or do not take parameter arguments
#something always needs to be passed to BVP with the right arguments
def f_dummy(T, Y):
    pass

def fp_dummy(T,P, Y):
    pass

def bc_dummy(Ya, Yb):
    pass
def bcp_dummy(Ya, Yb, P):
    pass

def fderive_dummy(T, Y):
    pass
def fderivep_dummy(T,P, Y):
    pass

def bcderive_dummy(Ya, Yb):
    pass
def bcderivep_dummy(Ya, Yb, P):
    pass
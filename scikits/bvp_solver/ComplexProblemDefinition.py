import numpy as np
from . import ProblemDefinition

class ComplexProblemDefinition(ProblemDefinition):
    """Defines a complex boundary value problem.
        This is wrapper for complex differential equations:
            dy/dx = f(x,y)
            where x is real and y(x) and f can be complex. The original solver takes only real values, 
            so we decompose the equations into real and imaginary parts. We use the following for a complex
            vector y:

            y=[y0, y1, ..., yN1] -> [y0.real,y1.real,...,yN1.real, y0.imag, y1.imag, ..., yN1.imag],

            where yN1 is the N-th element of y.
            
            Linear ODEs:
            For N linear ODEs, y' = My, where M
            is a NxN matrix, we have 2 x N real equations. That is 
            dy_big/dx = M_big y_big,
            where y_big = [y_r, y_i] and 
            M = [ M_r  -M_i ]
                [ M_i   M_r ].

        Note:
        1. The class ``ComplexProblemDefinition" basically replaces ProblemDefinition and allows users to
            define functions with complex arguments. 
        2. The original ``solve" function is used. Therefore, the ``parameter_guess" argument for ``solve"
            should be real. 
    """

    def __init__(self,
                num_ODE_c,
                num_parameters_c,
                num_left_boundary_conditions_c,
                boundary_points,
                function_c,
                boundary_conditions_c,
                function_derivative_c = None,
                boundary_conditions_derivative_c = None):
        
        self._num_ODE_c = num_parameters_c
        self._num_parameters_c = num_parameters_c
        self._num_left_boundary_conditions_c = num_left_boundary_conditions_c
        self._function_c = function_c
        self._boundary_conditions_c = boundary_conditions_c
        self._function_derivative_c = function_derivative_c
        self._boundary_conditions_derivative_c = boundary_conditions_derivative_c

        if self._num_parameters_c == 0:
            self._function_r = self._function_dummy
            self._boundary_conditions_r = self._boundary_conditions_dummy
            if self._function_derivative_c is None:
                self._function_derivative_r = None
            else:
                self._function_derivative_r = self._function_derivative_dummy
            if self._boundary_conditions_derivative_c is None:
                self._boundary_conditions_derivative_r = None
            else:
                self._boundary_conditions_derivative_r = self._boundary_conditions_derivative_dummy
            
        else:
            self._function_r = self._functionp_dummy
            self._boundary_conditions_r = self._boundary_conditionsp_dummy

            if self._function_derivative_c is None:
                self._function_derivative_r = None
            else:
                self._function_derivative_r = self._function_derivativep_dummy

            if self._boundary_conditions_derivative_c is None:
                self._boundary_conditions_derivative_r = None
            else:
                self._boundary_conditions_derivative_r = self._boundary_conditionsp_derivative_dummy

        ProblemDefinition.__init__(self,
                num_ODE = num_ODE_c * 2,
                num_parameters = num_parameters_c * 2,
                num_left_boundary_conditions = num_left_boundary_conditions_c * 2,
                boundary_points = boundary_points,
                function = self._function_r,
                boundary_conditions = self._boundary_conditions_r,
                function_derivative = self._function_derivative_r,
                boundary_conditions_derivative = self._boundary_conditions_derivative_r)

    def _function_dummy(self, T, Y):
        """dummy real function without parameter
        """
        dYc = self._function_c(T, real_to_complex(Y))
        return complex_to_real(dYc)

    def _functionp_dummy(self, T, P, Y):
        """dummy real function with parameter
        """
        dYc = self._function_c(T, real_to_complex(P), real_to_complex(Y))
        return complex_to_real(dYc)

    def _boundary_conditions_dummy(self, Ya, Yb):
        """dummy real BC without parameter
        """
        BCac, BCbc = self._boundary_conditions_c(
            real_to_complex(Ya), 
            real_to_complex(Yb))
        return complex_to_real(BCac), complex_to_real(BCbc)

    def _boundary_conditionsp_dummy(self, Ya, Yb, P):
        """dummy real BC with parameter
        """
        BCac, BCbc = self._boundary_conditions_c(
            real_to_complex(Ya), 
            real_to_complex(Yb),
            real_to_complex(P) )
        return complex_to_real(BCac), complex_to_real(BCbc)

    def _function_derivative_dummy(self, T, Y):
        """function derivative w/o parameter
        """
        dFdYc = self._function_derivative_c(T, real_to_complex(Y))
        return complex_to_real_matrix(dFdYc)

    def _function_derivativep_dummy(self, T, Y, P):
        """function derivative w/ parameter
        """
        dFdYc, dFdPc = self._function_derivative_c(T, real_to_complex(Y), real_to_complex(P))
        return complex_to_real_matrix(dFdYc), complex_to_real_matrix(dFdPc)

    def _boundary_conditions_derivative_dummy(self, Ya, Yb):
        """boundary conditions derivative w/o parameter
        """
        dBCa, dBCb = self._boundary_conditions_derivative_c(
                                real_to_complex(Ya), 
                                real_to_complex(Yb))
        return complex_to_real_matrix(dBCa), complex_to_real_matrix(dBCb)
    
    def _boundary_conditionsp_derivative_dummy(self, Ya, Yb, P):
        """boundary conditions derivative w/ parameter
        """

        dBCa, dBCb, dBPa, dBPb = self._boundary_conditions_derivative_c(
                                    real_to_complex(Ya), 
                                    real_to_complex(Yb),
                                    real_to_complex(P)
                                    )
        return complex_to_real_matrix(dBCa), complex_to_real_matrix(dBCb), \
            complex_to_real_matrix(dBPa), complex_to_real_matrix(dBPb)

# some helper functions
def complex_to_real(Yc):
    """convert a complex vector into real vector 
    """
    return np.hstack([Yc.real, Yc.imag])

def real_to_complex(Y):
    """convert a real vector into a complex vector
    """
    N = Y.shape[0]
    return Y[:int(N/2)] + 1j*Y[int(N/2):]

def complex_to_real_matrix(Mc):
    """compute an ``enlarged" real matrix
    """
    dim1, dim2 = Mc.shape
    M = np.zeros((dim1*2, dim2*2), dtype=float)
    M[:dim1,:dim2] = Mc.real
    M[:dim1,dim2:] =-Mc.imag
    M[dim1:,:dim2] = Mc.imag
    M[dim1:,dim2:] = Mc.real
    return M

def real_to_complex_matrix(M):
    dim1, dim2 = M.shape
    Mc = M[:int(dim1/2),:int(dim2/2)] - 1j * M[:int(dim1/2),int(dim2/2):]
    return Mc






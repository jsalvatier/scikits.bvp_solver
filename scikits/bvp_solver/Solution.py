import numpy
import tools
import bvp_solverf
import pickle


class Solution:
    """Stores the results of a :func:`solve` operation and provides a mechanism to evaluate valid results.
    """

    def __init__(self,
                 mesh,   #the current mesh
                 solution,   #the current solution
                 parameters = None,
                 work = None,
                 iwork = None,
                 yerror = None,
                 successIndicator = None,
                 extended = False):
        """
        """
        self._mesh = tools.preparg(mesh)
        self._solution = tools.preparg(solution)

        self._parameters = tools.preparg(parameters)
        self._work = tools.preparg(work)
        self._iwork = tools.preparg(iwork)
        self._yerror = yerror
        self._successIndicator = successIndicator
        self._extended = extended

    @property
    def mesh(self):
        """
        read only, ndarray, shape(N)
            Final mesh used by the solver for finding the approximate solution.
        """
        return self._mesh

    @property
    def solution(self):
        """
        read only, ndarray, shape(num_ODE, N)
            Approximate solution evaluated at the mesh points.
        """
        return self._solution

    @property
    def parameters(self):
        """
        read only, ndarray, shape(num_parameters)
            Approximate solution parameters.
        """
        return self._parameters

    @property
    def work(self):
        """
        read only, ndarray
            Working array generated during solving.
        """
        return self._work

    @property
    def iwork(self):
        """
        read only, ndarray
            Working integer array generated during solving.
        """
        return self._iwork

    @property
    def yerror(self):
        """
        read only, ndarray, shape(N)
            Approximate defect size evaluated at the mesh points.
        """
        return self._yerror

    @property
    def successIndicator(self):
        """
        read only, int
            Indicates success or failure of the solving operation.
        """
        return self._successIndicator

    @property
    def extended(self):
        """
        read only, logical
            Indicates whether the solution has been extended.
        """
        return self._extended

    def __call__(self, points, eval_derivative = False):
        """Evaluates the approximate solution and optionally the first derivative at an array of points.

        Parameters
        ----------
        points : castable to floating point ndarray, shape (N)
            Array of points where the approximate solution and derivative should be evaluated.
        eval_derivative : logical
            Determines whether the first derivative should be returned.

        Returns
        -------
        S : floating point ndarray, shape(num_ODE,N)
            Vector containing the approximate solution evaluated at points. Variable index is first, point index is second.
        D : floating point ndarray, shape(num_ODE,N)
            Vector of the first derivative to the approximate solution evaluated at points. Variable index is first, point index is second (only returned if eval_derivative = True).

        Raises
        ------
        ValueError
            If the approximate solution cannot be evaluated.
        ValueError
            If the Solution is the result of an :meth:`extend` operation.
        """
        if self._successIndicator == -1:
            raise ValueError("Solution is the result of a failed run, cannot evaluate")

        if self._extended == True:
            raise ValueError("""this solution is the result of extending a previous solution and cannot be evaluated.
             If you really want to know what the solution looks like, look at .mesh and .solution""")

        # if any of the points are more than a certain tolerance outside of the bounds, something has gone wrong
        tol = 1e-12
        dist = tol * (self._mesh[-1] - self._mesh[0])
        if (points < self._mesh[0] - dist).any() or (points > self._mesh[-1] + dist).any():
            raise ValueError("some points are outside the bounds of the solution")

        npar = 0
        if not (self._parameters is None):
            npar = len(self._parameters)

        bvp_solverf.bvp.bvp_eval_wrap(eval_derivative = eval_derivative,
                        points = tools.farg(points),
                        node_in = self._solution.shape[0],
                        npar_in = npar,
                        leftbc_in = 1, # this value doesn't matter
                        npts_in = len(self._mesh),
                        info_in = self._successIndicator,
                        mxnsub_in = 300, # nor does this
                        x_in = tools.farg(self._mesh),
                        y_in = tools.farg(self._solution),
                        parameters_in = tools.farg(self._parameters),
                        work_in = tools.farg(self._work),
                        iwork_in = tools.farg(self._iwork))

        #would prefer not to copy arrays here
        # but the results end up getting screwed up when this function is called again for some other purpose
        # this is probably because the old arrays "deallocated" by Fortran

        if eval_derivative:
            return bvp_solverf.bvp.evaluated.copy(), bvp_solverf.bvp.evaluated_d.copy()
        else:
            return bvp_solverf.bvp.evaluated.copy()


    def extend(self, new_left, new_right, order = 0, new_parameters = None):
        """Extends the solution to a new domain using polynomial extrapolation.

        Parameters
        ----------
        new_left : float
            Location of the new left boundary point.
        new_right : float
            Location of the new right boundary point.
        order : int
            Order of the polynomial to use for extrapolation
        new_parameters : castable to floating point ndarray
            Value of new parameters to use.

        Returns
        -------
        extended : Solution
            Extended :class:`Solution`.
        """
        new_parameters = tools.preparg(new_parameters)
        if new_parameters is None:
            new_parameters = self._parameters

        npar = 0
        if not (self._parameters is None):
            npar = len(self._parameters)

        bvp_solverf.bvp.bvp_extend_wrap(node_in = self._solution.shape[0],
                                        npar_in = npar,
                                        leftbc_in = 1, # this value doesn't matter
                                        npts_in = len(self._mesh),
                                        info_in = self._successIndicator,
                                        mxnsub_in = 300, # nor does this
                                        x_in = tools.farg(self._mesh),
                                        y_in = tools.farg(self._solution),
                                        parameters_in = tools.farg(self._parameters),
                                        work_in = tools.farg(self._work),
                                        iwork_in = tools.farg(self._iwork),
                                        anew = new_left,
                                        bnew = new_right,
                                        order = order,
                                        p = tools.farg(new_parameters),
                                        max_num_subintervals = 300  # nor does this
                                        )

        result = self.from_arg_list(bvp_solverf.bvp)
        result._extended = True
        return result

    @staticmethod
    def from_arg_list( bvp_object):
        """
        Gets the results of the fortran code from the fortran object, and returns a Soltuion object created from them.
        This is necessary because 2d arrays cannot be passed back to Python so the fortran code must
        store the results in a temporary area where we can get them later.
        """

        #would prefer not to copy arrays here
        # but the results end up getting screwed up when this function is called again for some other purpose
        # this is probably because the old arrays "deallocated" by Fortran
        new = Solution(mesh = tools.fromf(bvp_object.x),
                      solution = tools.fromf(bvp_object.y),
                      parameters = tools.fromf(bvp_object.parameters),
                      work = tools.fromf(bvp_object.work),
                      iwork = tools.fromf(bvp_object.iwork),
                      yerror = tools.fromf(bvp_object.yerror),
                      successIndicator = tools.fromf(bvp_object.info))
        return new

    @staticmethod
    def load(path):
        """Loads a Solution object from the file at the given path.

        Parameters
        ----------
        path : string
            Path of the file containing the solution to load.

        Returns
        -------
        sol : Solution
            Solution loaded from the file.
        """
        loadfile = open(path, "r")
        solution = pickle.load(loadfile)
        loadfile.close()
        return solution

    def save(self, path):
        """Saves the Solution object to a file given by the path.

        Parameters
        ----------
        path : string
            Path of the file to store the Solution in.
        """
        savefile = open(path, "w")
        pickle.dump(self, savefile)
        savefile.close()
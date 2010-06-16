.. currentmodule:: scikits.bvp_solver

============================
:func:`scikits.bvp_solver.get_template`
============================

Because the problem specification format for solving boundary value problems somewhat complex and arbitrary scikits.bvp_solver comes with a template generating function, :func:`get_template`, that will generate a code skeleton for a problem that can be filled in by the user.

The easiest way to use get_template is to call it in a shell. For example in IPython::

	In [1]: from scikits.bvp_solver import get_template

	In [2]: print get_template(num_ODE = 3, num_parameters = 1, num_left_boundary_conditions = 1,
	function_derivative = True, boundary_conditions_derivative = True, singular_term = True)

This prints a string containing code that can be copy-pasted into a script and then edited with appropriate formulas and variables. The code skeleton that these parameters generates can be seen :doc:`here <examples/examples.templateExampleOutput>`.

.. autofunction:: get_template
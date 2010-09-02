.. currentmodule:: scikits.bvp_solver

Welcome
=====================================

:mod:`scikits.bvp_solver` is a python package for solving two point boundary value problems which is based on a modified version of the `BVP_SOLVER <http://cs.stmarys.ca/~muir/BVP_SOLVER_Webpage.shtml>`_ Fortran package. If you have any questions, comments or suggestions about this tutorial, the examples or bvp_solver itself, please e-mail them to the mailing list or to me at jsalvati@u.washington.edu.

To join the mailing list send an e-mail to scikits-bvp_solver+subscribe@googlegroups.com

Installing and learning to use :mod:`scikits.bvp_solver`
--------------------------------------------------------
:mod:`scikits.bvp_solver` is available through `PyPi <http://pypi.python.org/pypi/scikits.bvp_solver>`_. The easiest way to learn how to install and use :mod:`scikits.bvp_solver` is to read the :doc:`tutorial <tutorial>`. It is also helpful to look at the :doc:`examples <examples/examples>`, and to read about the :doc:`template generator <scikits.bvp_solver.get_template>`, which will generate a code skeleton for a boundary value problem which can then be filled in. Using the template generator reduces the busywork of solving a boundary value problem. 

Documentation
-------------

.. toctree::
   :maxdepth: 1

   tutorial
   examples/examples
   core

The `BVP_SOLVER webpage <http://cs.stmarys.ca/~muir/BVP_SOLVER_Webpage.shtml>`_ has several more Fortran examples which should translate easily as well as a `paper <http://cs.stmarys.ca/~muir/JNAIAM_Shampine_Muir_Xu2006.pdf>`_ on the solver which describes its usage and capabilities in greater detail.  

Compilation Help
-----------------
:mod:`scikits.bvp_solver` requires the gfortran compiler; it may work with other f90 compilers, but this has not been tested.

Compiling on Windows
---------------------
To install on Windows (tested on Windows 7):

#. Get MinGW with gfortran `here <http://www.equation.com/servlet/equation.cmd?fa=fortran>`_

#. Compile from source using ``python setup.py config --compiler=mingw32 build --compiler=mingw32 install``
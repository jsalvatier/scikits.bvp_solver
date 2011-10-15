#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# This script solves the concentration of a reactant by solving the axial 
# dispersion model for a first order reaction. The model is represented 
# by a second order differential equation:                               
#                                                                         
# uS*(dCi/dx) = DL*(d^2Ci/dx^2) + Ri                                      
#                                                                         
# uS = superficial velocity, DL = axial dispersion coefficient,      
# Ri = reaction rate and Ci = concentration of i                         
#                                                                         
# The ODE obeys the following boundary conditions:                        
#                                                                         
# Ci0 = Ci - DL/uS*(dCi/dx)         at        x=(+0)                      
# dCi/dx = 0                               at        x=L                         
#                                                                         
# where L is the total length of the reactor.                             
# For convinience the ODE is scaled by making vi = Ci/Ci0 and z = x/L      
# and by making dvi = -dXi the dependent variable represents the          
# conversion of the reactant from Ci = Ci0*(1-Xi)                         
#                                                                         
# The ODE is now:                                                         
#                                                                         
# d^2vi/dz^2 = Pe(dvi/dz + Ri*L/Ci0/uS)      where       Pe =uS*L/DL            
#                                                                         
# with Boundary Conditions (BC) given by:                                 
#                                                                         
# dvi/dz = Pe*vi              at         z = (+0)                         
# dvi/dz = 0.0                at         z = (1)                          
#                                                                         
# The problem is to find the concentration of the reactant at the outlet  
# for a firts order reaction Ri = k*Ci0*(1-vi) with the following
# parameters Pe = 16, Ci0 = 1 and k*L/uS = kt = 2.0.
# The problem was taken from E. B. Nauman, Chemical Reactor
# Design, Optimization, and Scaleup, McGraw-Hill, 2001.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

import scikits.bvp_solver
from numpy import zeros, array, linspace
import pylab as pl

# Definition of the problem parameters

PeL = 16.0
Ci0 = 1.0
kt = 2.0

# By converting the original problem into a set of 2 firts order
# ODE's. Reaction rate is then Ri = -k*Ci0*(1-vi) the problem is defined as:

dvdz = zeros(2)
def systemode(z, v):
    dvdz[0] = v[1]
    dvdz[1] = PeL*(v[1]-kt*(1.0-v[0]))
    return dvdz

# Defining the BC on the left (z=0) and the right (z=1) side:

BCleft =  zeros(1)
BCright = zeros(1)
def boundary_conditions(va, vb):
    BCleft[0] = va[1] - PeL * va[0]
    BCright[0] = vb[1]
    return (BCleft), (BCright)

# Problem definition: a set of 2 first order ODE with a BC on the left side.
# Because the total length was scaled, the solution is reduced within an
# interval of (0.0 , 1.0)

problem = scikits.bvp_solver.ProblemDefinition(num_ODE = 2,
                                      num_parameters = 0,
                                      num_left_boundary_conditions = 1,
                                      boundary_points = (0.0, 1.0),
                                      function = systemode,
                                      boundary_conditions = boundary_conditions)

# By using conversion instead of concentration the solution is narrowed to
# an interval where the minimum value is 0.0 and the maximum is 1.0.
# Thus (0.0 , 1.0) is a sensible guess to initialize the problem

guess = array ([0.0 , 1.0])

# Once the problem and the BC's are being defined the solution of the problem
# is done by calling

solution = scikits.bvp_solver.solve(problem, solution_guess = guess, max_subintervals = 300)

# Getting the solution from z =0 to z = 1.0

z = linspace (0.0, 1.0, 100, endpoint = True )
v = solution(z)

# Comparison of the results

print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
print 'The analytical solution of the problem at the outlet is: '
print 'Final Conversion=  0.83605373297390484' 
print 'Final Concentration = 0.16394626702609519'
print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
print 'The numerical solution of the problem at the outlet is: '
print 'Final Conversion= ', v[0,99]
print 'Final Concentration = ', 1.0- v[0,99]
print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

# Plotting of the solution along z
pl.figure()
pl.subplot(1,2,1)
pl.plot(z, v[1,:], 'b-', z, v[0,:], 'g-')
pl.legend(('dv/dz[0]', 'dv/dz[1]'))
pl.xlabel('Normalized Reactor Length')
pl.ylabel('dv/dz')
pl.title('System of ODEs for the Axial Dispersion Model')
pl.grid('True')

pl.subplot(1,2,2)
pl.plot(z, v[0,:], 'g-')
pl.xlabel('Normalized Reactor Length')
pl.ylabel('Conversion of Reactant')
pl.title('Conversion produced by the Axial Dispersion Model')
pl.grid('True')
pl.show()

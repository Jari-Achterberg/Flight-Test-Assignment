import numpy as np
import control as ctrl
import scipy.io as sio
import matplotlib.pyplot as plt
from math import pi, sin, cos
import cmath

from Cit_par import statespacematrix
hp0    = 1500       	     # pressure altitude in the stationary flight condition [m]
V0 = 0.000000000000000000000000000001  # true airspeed in the stationary flight condition [m/sec]
alpha0 =  0.04           # angle of attack in the stationary flight condition [rad]
th0    =  0.08           # pitch angle in the stationary flight condition [rad]
m      = 10000           # Aircraft mass

print(statespacematrix(hp0,V0,alpha0,th0,m))



from Verification import get_eigenvalues, get_coefficients_3x3_matrix,abcd_formula, abc_formula, get_short_period_eigenvalues,get_phugoid_eigenvalues,get_aperiodic_roll_eigenvalue,get_dutch_roll_eigenvalues, get_spiral_motion_eigenvalue, get_dutch_roll_with_aperiodic_roll_eigenvalues

A= np.array([[2,-1,-1,0],
             [-1,3,-1,-1],
             [-1,-1,3,-1],
             [0,-1,-1,2]])
print(get_eigenvalues(A))    #should be 0,2,4   check

#get_coefficients_3x3_matrix(B)

abcd_formula(1, 0, -9, 0)  #x^{3}-9x=0 -> -3,03  check

print(abc_formula(1, 2, -8))  #2,4  check

print(get_short_period_eigenvalues(-1, 0, 0.5, 2, 3,0, -1,1)) #2,4  check

print(get_phugoid_eigenvalues(-2, -4, 0.5, 2))  #2,4  check

print(get_aperiodic_roll_eigenvalue(8, 2, 0.25))  #4  check

print(get_dutch_roll_eigenvalues(1,1,-80,-80,-80))   #behaves as it should check

print(get_spiral_motion_eigenvalue(2, 2, 1, 2, 2, 2, 0, 0, 0.25)) #4  check

print(get_dutch_roll_with_aperiodic_roll_eigenvalues(2, 2, 1, 2, -80, -80, -80, -80, -80, -80))  #behaves as it should check

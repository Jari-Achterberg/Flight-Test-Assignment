import numpy as np
import control as ctrl
from Cit_par import statespacematrix
import scipy.io as sio
import matplotlib.pyplot as plt

# Stationary flight condition
hp0    = 1500      	     # pressure altitude in the stationary flight condition [m]
V0     = 82.3            # true airspeed in the stationary flight condition [m/sec]
alpha0 =  0.04           # angle of attack in the stationary flight condition [rad]
th0    =  0.08           # pitch angle in the stationary flight condition [rad]

# Aircraft mass
m      = 10000            # mass [kg]

A_sym,B_sym,A_asym,B_asym = statespacematrix(hp0,V0,alpha0,th0)     #Calling state space matrix

#-----------Importing matlab flight data-------------------

mat_contents = sio.loadmat('FTISxprt-20200305_flight3.mat')
#Format is the following
elevator_dte = mat_contents['flightdata'][0][0][1][0][0][0]
vane_AOA = mat_contents['flightdata'][0][0][0][0][0][0]
time = mat_contents['flightdata'][0][0][-1][0][0][0]
pitch_rate = mat_contents['flightdata'][0][0][27][0][0][0]

#plt.plot(time.T,pitch_rate)
#---------------------------------------------------------

#Phugoid
C = A_sym[2,:]
D = B_sym[2,:]

def State_Space_output(A_sym,B_sym,C,D):
    sys = ctrl.ss(A_sym,B_sym,C,D)
    return sys

T, yout = ctrl.step_response(State_Space_output(A_sym,B_sym,C,D), T=np.arange(0, 250, 0.1))
plt.plot(T,yout)
plt.show()
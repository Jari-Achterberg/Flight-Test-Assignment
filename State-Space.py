import numpy as np
import control as ctrl
from Cit_par import statespacematrix
import scipy.io as sio
import matplotlib.pyplot as plt

# Stationary flight condition
hp0    = 1500      	     # pressure altitude in the stationary flight condition [m]
V0     = 82.3            # true airspeed in the stationary flight condition [m/sec]
alpha0 =  0.04          # angle of attack in the stationary flight condition [rad]
th0    =  0.08           # pitch angle in the stationary flight condition [rad]

# Aircraft mass
m      = 17000            # mass [kg]

A_sym,B_sym,A_asym,B_asym = statespacematrix(hp0,V0,alpha0,th0)     #Calling state space matrix

#-----------Importing matlab flight data-------------------

mat_contents = sio.loadmat('FTISxprt-20200305_flight3.mat')

#Importing values from the matlab file
elevator_dte = mat_contents['flightdata'][0][0][1][0][0][0]         #Deflection of the trim tab elevator in degree
vane_AOA = mat_contents['flightdata'][0][0][0][0][0][0]
time = mat_contents['flightdata'][0][0][-1][0][0][0]
pitch_rate = mat_contents['flightdata'][0][0][27][0][0][0]
pitch = mat_contents['flightdata'][0][0][22][0][0][0]
deltae = mat_contents['flightdata'][0][0][17][0][0][0]               #Deflection of the elevator in degree

plt.figure()
plt.plot(time.T,pitch)
plt.show()
#---------------------------------------------------------

#Phugoid

#Generating an output vector with pitch angle and pitch rate
C = np.array([[1,0,0,0],
              [0,1,0,0],
              [0,0,1,0],
              [0,0,0,1]])
D = np.array([[0],
              [0],
              [0],
              [0]])

sys = ctrl.ss(A_sym,B_sym,C,D)
ctrl.damp(sys)

#Getting the output parameters for a step input

T, u_sim = ctrl.step_response(sys, T=np.arange(0,600,0.1), X0=0.0, output=0)
T, aoa_sim = ctrl.step_response(sys, T=np.arange(0,600,0.1), X0=0.0, output=1)
T, pitch_sim = ctrl.step_response(sys, T=np.arange(0,600,0.1), X0=0.0, output=2)
T, pitch_rate_sim = ctrl.step_response(sys, T=np.arange(0,600,0.1), X0=0.0, output=3)

plt.figure(1)
plt.subplot(2,2,1)
plt.plot(T, pitch_sim)
plt.xlabel('Kutje')
plt.subplot(2,2,2)
plt.plot(T, pitch_rate_sim)
plt.subplot(2,2,3)
plt.plot(T, u_sim)
plt.subplot(2,2,4)
plt.plot(T, aoa_sim)
plt.show(1)
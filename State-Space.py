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
elevator_dte = mat_contents['flightdata'][0][0][1][0][0][0]             #Deflection of the trim tab elevator in degree
vane_AOA = mat_contents['flightdata'][0][0][0][0][0][0]                 #Angle of attack in degree
time = mat_contents['flightdata'][0][0][-1][0][0][0]                    #Time in seconds
pitch_rate = mat_contents['flightdata'][0][0][27][0][0][0]              #Pitch rate in deg per seconds
pitch = mat_contents['flightdata'][0][0][22][0][0][0]                   #Pitch in degrees
deltae = mat_contents['flightdata'][0][0][17][0][0][0]                  #Deflection of the elevator in degree
Vtrue = mat_contents['flightdata'][0][0][42][0][0][0]                   #True airspeed in knots

Vtrue = Vtrue * 0.51444444444444444                                     #True airspeed in m/s

#plt.figure(2)
#plt.plot(time.T,Vtrue)
#plt.show(2)


plt.figure(2)
plt.plot(time.T,vane_AOA, 'b', time.T, pitch, 'r')
plt.show(2)

#plt.figure()
#plt.plot(time.T, pitch)
#plt.show()

plt.figure(3)
plt.plot(time.T,deltae)
plt.show(3)


#---------------------------------------------------------

#---------------------Phugoid------------------------------

#Getting the output parameters for a step input

xinit = np.array([Vtrue[25360]*np.cos(vane_AOA[25360]*(np.pi/180)), vane_AOA[25360]*(np.pi/180), pitch[25360]*(np.pi/180), pitch_rate[25360]*(np.pi/180)])
#xinit = np.array([0,0,0,0])
Udeltae = deltae[25360:26910]*(np.pi/180)

#Generating an output vector with pitch angle and pitch rate
C = np.array([[1,0,0,0]])
D = np.array([[0]])

sys = ctrl.ss(A_sym,B_sym,C,D)
ctrl.damp(sys)

T, u_sim, xout = ctrl.forced_response(sys, T=np.arange(0,155,0.1), U=Udeltae.T,  X0=xinit)

C = np.array([[0,1,0,0]])
D = np.array([[0]])

sys = ctrl.ss(A_sym,B_sym,C,D)

T, aoa_sim, xout = ctrl.forced_response(sys, T=np.arange(0,155,0.1), U=Udeltae.T, X0=xinit)

C = np.array([[0,0,1,0]])
D = np.array([[0]])

sys = ctrl.ss(A_sym,B_sym,C,D)

T, pitch_sim, xout = ctrl.forced_response(sys, T=np.arange(0,155,0.1), U=Udeltae.T, X0=xinit)

C = np.array([[0,0,0,1]])
D = np.array([[0]])

sys = ctrl.ss(A_sym,B_sym,C,D)

T, pitch_rate_sim, xout = ctrl.forced_response(sys, T=np.arange(0,155,0.1), U=Udeltae.T, X0=xinit)
# OR USE FORCED RESPONSE -> PHUGOID -> short, sharp deflection followed by a return to the centered position

#Converting horizontal velocity to flight velocity
V_sim = u_sim / np.cos(aoa_sim)

#Converting radians to degree
pitch_sim = pitch_sim[:] * (180/np.pi)
pitch_rate_sim = pitch_rate_sim[:] * (180/np.pi)
aoa_sim = aoa_sim[:] * (180/np.pi)

plt.figure(1)
plt.subplot(2,2,1)
plt.plot(T, pitch_sim, T, pitch[25360:26910])
plt.grid()
plt.legend(['Simulation','Flight Data'],loc=4)
plt.ylabel('Pitch angle [degree]')
plt.xlabel('Time [s]')
plt.subplot(2,2,2)
plt.plot(T, pitch_rate_sim, T, pitch_rate[25360:26910])
plt.grid()
plt.legend(['Simulation','Flight Data'],loc=4)
plt.ylabel('Pitch rate [degree/sec]')
plt.xlabel('Time [s]')
plt.subplot(2,2,3)
plt.plot(T, V_sim, T, Vtrue[25360:26910])
plt.grid()
plt.legend(['Simulation','Flight Data'],loc=4)
plt.ylabel('Flight velocity [m/s]')
plt.xlabel('Time [s]')
plt.subplot(2,2,4)
plt.plot(T, aoa_sim, T, vane_AOA[25360:26910])
plt.grid()
plt.legend(['Simulation','Flight Data'], loc=4)
plt.ylabel('Angle of attack [degree]')
plt.xlabel('Time [s]')
plt.show(1)

#------------------------Short Period----------------------------------------

import numpy as np
import control as ctrl
from Cit_par import statespacematrix
import scipy.io as sio
import matplotlib.pyplot as plt

#-----------Importing matlab flight data-------------------

mat_contents = sio.loadmat('FTISxprt-20200305_flight3.mat')
vane_AOA = mat_contents['flightdata'][0][0][0][0][0][0]                 #Angle of attack in degree
time = mat_contents['flightdata'][0][0][-1][0][0][0]                    #Time in seconds
pitch_rate = mat_contents['flightdata'][0][0][27][0][0][0]              #Pitch rate in deg per seconds
pitch = mat_contents['flightdata'][0][0][22][0][0][0]                   #Pitch in degrees
deltae = mat_contents['flightdata'][0][0][17][0][0][0]                  #Deflection of the elevator in degree
Vtrue = mat_contents['flightdata'][0][0][42][0][0][0]                   #True airspeed in knots
Pressure_Altitude = mat_contents['flightdata'][0][0][37][0][0][0]       #Pressure altitude in ft

#Conversion to SI units
Pressure_Altitude = Pressure_Altitude * 0.3048                          #Pressure altitude in m
Vtrue = Vtrue * 0.51444444444444444                                     #True airspeed in m/s

#---------------------Phugoid------------------------------

#Short, sharp deflection followed by a return to the centered position
startvalue =25380                   #Starts at index 25380
endvalue = 26930                    #Ends at index 26930

# Stationary flight condition
hp0    = Pressure_Altitude[startvalue]      # pressure altitude in the stationary flight condition [m]
V0     = Vtrue[startvalue]                  # true airspeed in the stationary flight condition [m/sec]
alpha0 = vane_AOA[startvalue]*(np.pi/180)   # angle of attack in the stationary flight condition [rad]
th0    = pitch[startvalue]*(np.pi/180)      # pitch angle in the stationary flight condition [rad]

# Aircraft mass
m      = 6406.644            # mass [kg]

A_sym,B_sym,A_asym,B_asym = statespacematrix(hp0[0],V0[0],alpha0[0],th0[0],m)     #Calling state space matrix

C = np.identity(4)
D = np.array([[0],
              [0],
              [0],
              [0]])

sys = ctrl.ss(A_sym,B_sym,C,D)
xinit = np.array([0,2*(np.pi/180),0,0])
T ,yout  = ctrl.initial_response(sys, T=np.arange(0,100,0.1), X0=xinit)

#Converting horizontal velocity to flight velocity
V_sim = yout[0] / np.cos(yout[1])

#Converting radians to degree

pitch_sim = yout[2] * (180/np.pi)
pitch_rate_sim = yout[3] * (180/np.pi)
aoa_sim = yout[1] * (180/np.pi)

plt.figure(1)
plt.title('Response curves for change in angle of attack')
plt.subplot(2,2,1)
plt.plot(T, pitch_sim)
plt.grid()
plt.ylabel('Pitch angle change [degree]')
plt.xlabel('Time [s]')
plt.subplot(2,2,2)
plt.plot(T, pitch_rate_sim)
plt.grid()
plt.ylabel('Pitch rate change [degree/sec]')
plt.xlabel('Time [s]')
plt.subplot(2,2,3)
plt.plot(T, V_sim)
plt.grid()
plt.ylabel('Flight velocity change [m/s]')
plt.xlabel('Time [s]')
plt.subplot(2,2,4)
plt.plot(T, aoa_sim)
plt.grid()
plt.ylabel('Angle of attack change [degree]')
plt.xlabel('Time [s]')
plt.show(1)
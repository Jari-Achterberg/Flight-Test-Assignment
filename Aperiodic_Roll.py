import numpy as np
import control as ctrl
from Cit_par import statespacematrix
import scipy.io as sio
import matplotlib.pyplot as plt

#-----------Importing matlab flight data-------------------

mat_contents = sio.loadmat('FTISxprt-20200305_flight3.mat')
vane_AOA = mat_contents['flightdata'][0][0][0][0][0][0]                 #Angle of attack in degree
time = mat_contents['flightdata'][0][0][-1][0][0][0]                    #Time in seconds
pitch = mat_contents['flightdata'][0][0][22][0][0][0]                   #Pitch in degrees
Vtrue = mat_contents['flightdata'][0][0][42][0][0][0]                   #True airspeed in knots
Pressure_Altitude = mat_contents['flightdata'][0][0][37][0][0][0]       #Pressure altitude in ft
deltar = mat_contents['flightdata'][0][0][18][0][0][0]                  #Deflection of the rudder in degree
deltaa = mat_contents['flightdata'][0][0][16][0][0][0]                  #Deflection of the aileron in degree
rollangle = mat_contents['flightdata'][0][0][21][0][0][0]               #Roll angle in degrees
rollrate =  mat_contents['flightdata'][0][0][26][0][0][0]               #Body roll rate in deg/s
yawrate = mat_contents['flightdata'][0][0][28][0][0][0]                 #Body yaw rate in deg/s

#Conversion to SI units
Pressure_Altitude = Pressure_Altitude * 0.3048                          #Pressure altitude in m
Vtrue = Vtrue * 0.51444444444444444                                     #True airspeed in m/s

#---------------------Aperiod roll------------------------------

#We have a pulse shaped aileron deflection on t = 2915 [s] for aperiodic roll

startvalue =29060                   #Starts at index 29060
endvalue = 29460                    #We want to see 40 seconds after start

# Stationary flight condition
hp0    = Pressure_Altitude[startvalue]      # pressure altitude in the stationary flight condition [m]
V0     = Vtrue[startvalue]                  # true airspeed in the stationary flight condition [m/sec]
alpha0 = vane_AOA[startvalue]*(np.pi/180)   # angle of attack in the stationary flight condition [rad]
th0    = pitch[startvalue]*(np.pi/180)      # pitch angle in the stationary flight condition [rad]

# Aircraft mass
m      = 17000            # mass [kg]

A_sym,B_sym,A_asym,B_asym = statespacematrix(hp0[0],V0[0],alpha0[0],th0[0],m)     #Calling state space matrix

xinit = np.array([0,0,0,0])
#xinit = np.array([0,rollangle[startvalue]*(np.pi/180),rollrate[startvalue]*(np.pi/180),yawrate[startvalue]*(np.pi/180)])

#Input from the rudder and aileron
Udeltar = -deltar[startvalue:endvalue]*(np.pi/180)      #Minus because of other sign convention?
Udeltaa = -deltaa[startvalue:endvalue]*(np.pi/180)      #Minus because of other sign convention?

#Generating an output vector with sideslip

C = np.array([[1,0,0,0]])
D = np.array([[0,0]])
sys = ctrl.ss(A_asym,B_asym,C,D)

T, sideslip_sim, xout = ctrl.forced_response(sys, T=np.arange(0,(endvalue-startvalue)/10,0.1), U=np.concatenate((Udeltaa.T, Udeltar.T), axis=0),  X0=xinit)
sideslip_sim = sideslip_sim*(180/np.pi)         #from radians to degree

#Generating an output vector with rolling angle

C = np.array([[0,1,0,0]])
D = np.array([[0,0]])
sys = ctrl.ss(A_asym,B_asym,C,D)

T, rollangle_sim, xout = ctrl.forced_response(sys, T=np.arange(0,(endvalue-startvalue)/10,0.1), U=np.concatenate((Udeltaa.T, Udeltar.T), axis=0),  X0=xinit)
rollangle_sim = rollangle_sim*(180/np.pi)         #from radians to degree
rollangle_sim = rollangle_sim + rollangle[startvalue]       #implementing initial condition

#Generating an output vector with roll rate

C = np.array([[0,0,1,0]])
D = np.array([[0,0]])
sys = ctrl.ss(A_asym,B_asym,C,D)

T, rollrate_sim, xout = ctrl.forced_response(sys, T=np.arange(0,(endvalue-startvalue)/10,0.1), U=np.concatenate((Udeltaa.T, Udeltar.T), axis=0),  X0=xinit)
rollrate_sim = rollrate_sim*(180/np.pi)         #from radians/s to degree/s
rollrate_sim = rollrate_sim + rollrate[startvalue]       #implementing initial condition

#Generating an output vector with roll rate

C = np.array([[0,0,0,1]])
D = np.array([[0,0]])
sys = ctrl.ss(A_asym,B_asym,C,D)

T, yawrate_sim, xout = ctrl.forced_response(sys, T=np.arange(0,(endvalue-startvalue)/10,0.1), U=np.concatenate((Udeltaa.T, Udeltar.T), axis=0),  X0=xinit)
yawrate_sim = yawrate_sim*(180/np.pi)                   #from radians/s to degree/s
yawrate_sim = yawrate_sim + yawrate[startvalue]       #implementing initial condition

plt.figure()
plt.title('Response curves for a pulse-shaped aileron deflection, aperiodic roll')
plt.subplot(2,2,1)
plt.plot(T, sideslip_sim)
plt.grid()
plt.legend(['Simulation','Flight Data'],loc=4)
plt.ylabel('Side slip [degree]')
plt.xlabel('Time [s]')
plt.subplot(2,2,2)
plt.plot(T, rollangle_sim, T, rollangle[startvalue:endvalue])
plt.grid()
plt.legend(['Simulation','Flight Data'],loc=4)
plt.ylabel('Roll angle [degree]')
plt.xlabel('Time [s]')
plt.subplot(2,2,3)
plt.plot(T, rollrate_sim, T, rollrate[startvalue:endvalue])
plt.grid()
plt.legend(['Simulation','Flight Data'],loc=4)
plt.ylabel('Roll rate [degree/s]')
plt.xlabel('Time [s]')
plt.subplot(2,2,4)
plt.plot(T, yawrate_sim, T, yawrate[startvalue:endvalue])
plt.grid()
plt.legend(['Simulation','Flight Data'], loc=4)
plt.ylabel('Yaw rate [degree/s]')
plt.xlabel('Time [s]')
plt.show()

plt.figure()
plt.plot(T,deltar[startvalue:endvalue], 'b', T, deltaa[startvalue:endvalue], 'r')
plt.show()
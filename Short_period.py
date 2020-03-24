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

#---------------------Short period------------------------------

#Short, sharp deflection followed by a return to the centered position
startvalue = 27305                  #Starts at index 25385 (Phugoid)/27305 (Short Period)
endvalue = 27405                    #Ends at index 25485 (Phugoid)/27405 (Short Period)

# Stationary flight condition
hp0    = Pressure_Altitude[startvalue]      # pressure altitude in the stationary flight condition [m]
V0     = Vtrue[startvalue]                  # true airspeed in the stationary flight condition [m/sec]
alpha0 = vane_AOA[startvalue]*(np.pi/180)   # angle of attack in the stationary flight condition [rad]
th0    = pitch[startvalue]*(np.pi/180)      # pitch angle in the stationary flight condition [rad]

# Aircraft mass
m      = 6385.191            # mass [kg]

A_sym,B_sym,A_asym,B_asym = statespacematrix(hp0[0],V0[0],alpha0[0],th0[0],m)     #Calling state space matrix

c = 2.0569
print("symmetric short period: ", np.linalg.eigvals(A_sym)[0:2])

#------------------------Short Period----------------------------------------

#xinit = np.array([Vtrue[25360]*np.cos(vane_AOA[25360]*(np.pi/180)), vane_AOA[25360]*(np.pi/180), pitch[25360]*(np.pi/180), pitch_rate[25360]*(np.pi/180)])
xinit = np.array([0,0,0,0])
Udeltae = deltae[startvalue:endvalue]*(np.pi/180)

#Generating an output vector with velocity

C = np.array([[1,0,0,0]])
D = np.array([[0]])
sys = ctrl.ss(A_sym,B_sym,C,D)

#Using input data of flight test
T, u_sim, xout = ctrl.forced_response(sys, T=np.arange(0,(endvalue-startvalue)/10,0.1), U=Udeltae.T,  X0=xinit)
u_sim = u_sim + Vtrue[startvalue]*np.cos(vane_AOA[startvalue]*(np.pi/180))            #For the initial condition

#Using step input on elevator
input = np.zeros(len(Udeltae))
input[0:100]=-1.35*(np.pi/180)

#Using a unit pulse input
T_, u_sim_, xout_ = ctrl.forced_response(sys, T=np.arange(0,(endvalue-startvalue)/10,0.1), U=input,  X0=xinit)
u_sim_ = u_sim_ + Vtrue[startvalue]*np.cos(vane_AOA[startvalue]*(np.pi/180))            #For the initial condition

#Generating an output vector with angle of attack

C = np.array([[0,1,0,0]])
D = np.array([[0]])

sys = ctrl.ss(A_sym,B_sym,C,D)

T, aoa_sim, xout = ctrl.forced_response(sys, T=np.arange(0,(endvalue-startvalue)/10,0.1), U=Udeltae.T, X0=xinit)
aoa_sim = aoa_sim + vane_AOA[startvalue]*(np.pi/180)                             #For the initial condition

#Using a unit pulse input
T_, aoa_sim_, xout_ = ctrl.forced_response(sys, T=np.arange(0,(endvalue-startvalue)/10,0.1), U=input, X0=xinit)
aoa_sim_ = aoa_sim_ + vane_AOA[startvalue]*(np.pi/180)                             #For the initial condition

#Generating an output vector with pitch angle

C = np.array([[0,0,1,0]])
D = np.array([[0]])

sys = ctrl.ss(A_sym,B_sym,C,D)

T, pitch_sim, xout = ctrl.forced_response(sys, T=np.arange(0,(endvalue-startvalue)/10,0.1), U=Udeltae.T, X0=xinit)
pitch_sim = pitch_sim + pitch[startvalue]*(np.pi/180)                            #For the initial condition

#Using a unit pulse input
T_, pitch_sim_, xout_ = ctrl.forced_response(sys, T=np.arange(0,(endvalue-startvalue)/10,0.1), U=input, X0=xinit)
pitch_sim_ = pitch_sim_ + pitch[startvalue]*(np.pi/180)                            #For the initial condition

#Generating an output vector with pitch rate

C = np.array([[0,0,0,1]])
D = np.array([[0]])

sys = ctrl.ss(A_sym,B_sym,C,D)

T, pitch_rate_sim, xout = ctrl.forced_response(sys, T=np.arange(0,(endvalue-startvalue)/10,0.1), U=Udeltae.T, X0=xinit)
pitch_rate_sim = pitch_rate_sim + pitch_rate[startvalue]*(np.pi/180)             #For the initial condition

#Using a unit pulse input
T_, pitch_rate_sim_, xout_ = ctrl.forced_response(sys, T=np.arange(0,(endvalue-startvalue)/10,0.1), U=input, X0=xinit)
pitch_rate_sim_ = pitch_rate_sim_ + pitch_rate[startvalue]*(np.pi/180)             #For the initial condition

#Converting horizontal velocity to flight velocity
V_sim = u_sim / np.cos(aoa_sim)
V_sim_ = u_sim_ / np.cos(aoa_sim_)

# Converting radians to degree
pitch_sim = pitch_sim[:] * (180/np.pi)
pitch_rate_sim = pitch_rate_sim[:] * (180/np.pi)
aoa_sim = aoa_sim[:] * (180/np.pi)
#
a, b = 6, -1.2
ePower1 = a*np.exp(b*T) - 0.32*T + 3.7
ePower2 = -a*np.exp(b*T) - 0.32*T + 3.7
Thalf = 0.533
P = 3.84 - 0.7  # P = 2*np.pi/1.3

plt.plot(T, ePower2, T, ePower1, T, pitch_rate[startvalue:endvalue])
plt.grid()
plt.legend(['epower','epower2','Flight Data'],loc=4)
plt.ylabel('Pitch rate [degree/sec]')
plt.xlabel('Time [s]')
plt.show()
####
eigenvalue_real = np.log(0.5) / Thalf
eigenvalue_imag = 2*np.pi/P
print("lambda = ", eigenvalue_real, "+ i", eigenvalue_imag)
print(V0)

pitch_sim_ = pitch_sim_[:] * (180/np.pi)
pitch_rate_sim_ = pitch_rate_sim_[:] * (180/np.pi)
aoa_sim_ = aoa_sim_[:] * (180/np.pi)

plt.figure(1)
plt.title('Response curves for a step elevator deflection, short period response')
plt.subplot(2,2,1)
plt.plot(T, pitch_sim, T, pitch[startvalue:endvalue])
plt.grid()
plt.legend(['Simulation','Flight Data'],loc=4)
plt.ylabel('Pitch angle [degree]')
plt.xlabel('Time [s]')
plt.subplot(2,2,2)
plt.plot(T, pitch_rate_sim, T, pitch_rate[startvalue:endvalue])
plt.grid()
plt.legend(['Simulation','Flight Data'],loc=4)
plt.ylabel('Pitch rate [degree/sec]')
plt.xlabel('Time [s]')
plt.subplot(2,2,3)
plt.plot(T, V_sim, T, Vtrue[startvalue:endvalue])
plt.grid()
plt.legend(['Simulation','Flight Data'],loc=4)
plt.ylabel('Flight velocity [m/s]')
plt.xlabel('Time [s]')
plt.subplot(2,2,4)
plt.plot(T, aoa_sim, T, vane_AOA[startvalue:endvalue])
plt.grid()
plt.legend(['Simulation','Flight Data'], loc=4)
plt.ylabel('Angle of attack [degree]')
plt.xlabel('Time [s]')
plt.show(1)

plt.figure(2)
plt.plot(T, deltae[startvalue:endvalue])
plt.show(2)

#For unit impulse!
plt.figure(3)
plt.title('Response curves for a step elevator deflection, short period response')
plt.subplot(2,2,1)
plt.plot(T, pitch_sim_)
plt.grid()
plt.ylabel('Pitch angle [degree]')
plt.xlabel('Time [s]')
plt.subplot(2,2,2)
plt.plot(T, pitch_rate_sim_)
plt.grid()
plt.ylabel('Pitch rate [degree/sec]')
plt.xlabel('Time [s]')
plt.subplot(2,2,3)
plt.plot(T, V_sim_)
plt.grid()
plt.ylabel('Flight velocity [m/s]')
plt.xlabel('Time [s]')
plt.subplot(2,2,4)
plt.plot(T, aoa_sim_)
plt.grid()
plt.ylabel('Angle of attack [degree]')
plt.xlabel('Time [s]')
plt.show(3)
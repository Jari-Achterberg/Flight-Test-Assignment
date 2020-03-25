import numpy as np
import control as ctrl
from Cit_par import statespacematrix
import scipy.io as sio
import matplotlib.pyplot as plt

# -----------Importing matlab flight data-------------------
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

# --------------------------------Spiral----------------------------------------------
# Spiral starts at t=3159 [s] and ends at t=3306 [s]
startvalue = 31590
endvalue = 33060

# Stationary flight condition
hp0    = Pressure_Altitude[startvalue]      # pressure altitude in the stationary flight condition [m]
V0     = Vtrue[startvalue]                  # true airspeed in the stationary flight condition [m/sec]
alpha0 = vane_AOA[startvalue]*(np.pi/180)   # angle of attack in the stationary flight condition [rad]
th0    = pitch[startvalue]*(np.pi/180)      # pitch angle in the stationary flight condition [rad]

# Aircraft mass
m      = 6332.428            # mass [kg]

A_sym,B_sym,A_asym,B_asym = statespacematrix(hp0[0],V0[0],alpha0[0],th0[0],m)     #Calling state space matrix
print("asymmetric spiral: ", np.linalg.eigvals(A_asym))

xinit = np.array([0,0,0,0])
#xinit = np.array([0,rollangle[startvalue]*(np.pi/180),rollrate[startvalue]*(np.pi/180),yawrate[startvalue]*(np.pi/180)])

#Input from the rudder and aileron
Udeltar = deltar[startvalue:endvalue]*(np.pi/180)
Udeltaa = deltaa[startvalue:endvalue]*(np.pi/180)

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

#Generating an output vector with yaw rate
C = np.array([[0,0,0,1]])
D = np.array([[0,0]])
sys = ctrl.ss(A_asym,B_asym,C,D)

T, yawrate_sim, xout = ctrl.forced_response(sys, T=np.arange(0,(endvalue-startvalue)/10,0.1), U=np.concatenate((Udeltaa.T, Udeltar.T), axis=0),  X0=xinit)
yawrate_sim = yawrate_sim*(180/np.pi)                   #from radians/s to degree/s
yawrate_sim = yawrate_sim + yawrate[startvalue]       #implementing initial condition

a, b = 11.8, +0.0092
TE = np.arange(-150, 147, 0.1)
ePower1 = a*np.exp(b*TE)
ePower2 = -a*np.exp(b*TE)
Thalf = -75.34 # -71.68 # = - 76.7
P = 10**99  # P = 2*np.pi/1.3

plt.figure()
plt.grid()
plt.plot(TE, ePower2, T, rollangle[startvalue:endvalue])
plt.legend(['e','Spiral Motion'],loc=1)
plt.ylabel('Roll angle [degree]')
plt.xlabel('Time [s]')
plt.show()

eigenvalue_real = np.log(0.5) / Thalf
eigenvalue_imag = 2*np.pi/P
print("lambda = ", eigenvalue_real, "+ i", eigenvalue_imag)

'''
plt.figure()
plt.grid()
plt.plot(T,rollangle_sim,T,rollangle[startvalue:endvalue])
plt.legend(['Simulation','Spiral Motion'],loc=1)
plt.ylabel('Roll angle [degree]')
plt.xlabel('Time [s]')
plt.show()

plt.figure()
plt.plot(T,deltar[startvalue:endvalue], 'b', T, deltaa[startvalue:endvalue], 'r')
plt.show()
#----
'''

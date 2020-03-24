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

#---------------------Dutch roll------------------------------

#We have a pulse shaped rudder deflection on t = 2792 [s] for dutch roll
#We have a pulse shaped rudder deflection on t = 2837 [s] for dutch roll with yaw damper

startvalue =27830                   #Starts at index 25360
endvalue = 28000                    #We want to see ... seconds after start

# Stationary flight condition
hp0    = Pressure_Altitude[startvalue]      # pressure altitude in the stationary flight condition [m]
V0     = Vtrue[startvalue]                  # true airspeed in the stationary flight condition [m/sec]
alpha0 = vane_AOA[startvalue]*(np.pi/180)   # angle of attack in the stationary flight condition [rad]
th0    = pitch[startvalue]*(np.pi/180)      # pitch angle in the stationary flight condition [rad]

# Aircraft mass
m      = 6378.821            # mass [kg]

A_sym,B_sym,A_asym,B_asym = statespacematrix(hp0[0],V0[0],alpha0[0],th0[0],m)     #Calling state space matrix

print("asymmetric dutchroll: ", np.linalg.eigvals(A_asym)[1:3])

xinit = np.array([0,0,0,0])
#xinit = np.array([0,rollangle[startvalue]*(np.pi/180),rollrate[startvalue]*(np.pi/180),yawrate[startvalue]*(np.pi/180)])

#Input from the rudder and aileron
Udeltar = -deltar[startvalue:endvalue]*(np.pi/180)
Udeltaa = -deltaa[startvalue:endvalue]*(np.pi/180)

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

#Generating an output vector with yaw rate

C = np.array([[0,0,0,1]])
D = np.array([[0,0]])
sys = ctrl.ss(A_asym,B_asym,C,D)

T, yawrate_sim, xout = ctrl.forced_response(sys, T=np.arange(0,(endvalue-startvalue)/10,0.1), U=np.concatenate((Udeltaa.T, Udeltar.T), axis=0),  X0=xinit)
yawrate_sim = yawrate_sim*(180/np.pi)                   #from radians/s to degree/s
yawrate_sim = yawrate_sim + yawrate[startvalue]       #implementing initial condition

#T, yawrate_sim1= ctrl.impulse_response(sys, T=np.arange(0,40,0.1), X0=xinit)

#----------------------------Dutch roll with yaw damping--------------------------------

startvalue1 = 28271
endvalue1 = startvalue1 + (endvalue-startvalue)


a, b = 55, -0.23
ePower1 = a*np.exp(b*T)     # - 0.32*T + 3.7
ePower2 = -a*np.exp(b*T)    # - 0.32*T + 3.7
Thalf = 3.01
# peak 1 and 5
amp1 = 3.62
amp2 = 16.52
P = (amp2 - amp1)/4
print("T: ", Thalf)
print("P: ", P)
plt.plot(T, ePower1, T, ePower2, T, yawrate[startvalue:endvalue])
plt.grid()
plt.legend(['e1','e2','Dutch Roll Flight'],loc=4)
plt.ylabel('Roll rate [degree/s]')
plt.xlabel('Time [s]')
plt.show()

print(V0)
eigenvalue_real = np.log(0.5) / Thalf
eigenvalue_imag = 2*np.pi/P
print("lambda = ", eigenvalue_real, "+ i", eigenvalue_imag)

'''
plt.figure(1)
plt.title('Response curves for a pulse-shaped rudder deflection, dutch roll')
plt.subplot(2,2,1)
plt.plot(T, sideslip_sim)
plt.grid()
plt.legend(['Simulation','Dutch Roll Flight', 'Dutch Roll Yaw Damper'],loc=4)
plt.ylabel('Side slip [degree]')
plt.xlabel('Time [s]')
plt.subplot(2,2,2)
plt.plot(T, rollangle_sim, T, rollangle[startvalue:endvalue], T, rollangle[startvalue1:endvalue1])
plt.grid()
plt.legend(['Simulation','Dutch Roll Flight', 'Dutch Roll Yaw Damper'],loc=4)
plt.ylabel('Roll angle [degree]')
plt.xlabel('Time [s]')
plt.subplot(2,2,3)
plt.plot(T, rollrate_sim, T, rollrate[startvalue:endvalue], T, rollrate[startvalue1:endvalue1])
plt.grid()
plt.legend(['Simulation','Dutch Roll Flight', 'Dutch Roll Yaw Damper'],loc=4)
plt.ylabel('Roll rate [degree/s]')
plt.xlabel('Time [s]')
plt.subplot(2,2,4)
plt.plot(T, yawrate_sim, T, yawrate[startvalue:endvalue], T, yawrate[startvalue1:endvalue1])
plt.grid()
plt.legend(['Simulation','Dutch Roll Flight', 'Dutch Roll Yaw Damper'], loc=4)
plt.ylabel('Yaw rate [degree/s]')
plt.xlabel('Time [s]')
plt.show()

plt.figure(2)
plt.grid()
plt.plot(T[0:250],rollrate[startvalue:endvalue],'k', linestyle='--')
plt.plot(T[0:250],yawrate[startvalue:endvalue], 'b',  linestyle= '-')
plt.legend(['Roll rate','Yaw rate'],loc=4)
plt.ylabel('p,r [degree/s]')
plt.xlabel('Time [s]')
plt.show(2)

plt.figure()
plt.plot(T,deltar[startvalue1:endvalue1], 'b', T, deltaa[startvalue1:endvalue1], 'r')
plt.show()
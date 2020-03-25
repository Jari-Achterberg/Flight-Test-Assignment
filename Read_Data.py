import numpy as np
import subprocess
import matplotlib.pyplot as plt

# mass = np.genfromtxt('Data_Test.csv', skip_header = 9,max_rows = 1)
F_block = np.loadtxt('Data_Test.csv', delimiter = ';', skiprows = 17, max_rows = 1, usecols = 1) # [lbs]
mass_pl = np.loadtxt('Data_Test.csv', delimiter = ';', skiprows = 7, max_rows = 9, usecols = 6) # read mass of all persons [kg]
CLCD = np.loadtxt('Data_Test.csv', delimiter = ';', skiprows = 27, max_rows = 6, usecols = (1,2,3,4,5,6,7,8,9)) # read data CL-CD measurements
trim = np.loadtxt('Data_Test.csv', delimiter = ';', skiprows = 57, max_rows = 7, usecols = (1,2,3,4,5,6,7,8,9,10,11,12)) # read data trim measurements
cg = np.loadtxt('Data_Test.csv', delimiter = ';', skiprows = 73, max_rows = 2, usecols = (1,2,3,4,5,6,7,8,9,10,11,12)) # read data cg measurements
# Format CLCD: time_min, time_sec, hp [ft], IAS [kts], a [deg], FFl [kg/h], FFr [kg/h], F_used [lbs], TAT [celsius]

# ---- General data ----
F_block = F_block * 0.45359237

# ---- CL-CD measurements ----
t_CLCD = CLCD[:,0] * 60 + CLCD[:,1] # total time in [s] for CL-CD measurements
hp_CLCD = CLCD[:,2] * 0.3048 # hp [m] for CL-CD measurements
IAS_CLCD = CLCD[:,3] * 0.514444444 # indicated airspeed [m/s] for CL-CD measurements
alpha_CLCD = np.deg2rad(CLCD[:,4]) # AoA in [rad] for CL-CD measurements
FFl_CLCD = CLCD[:,5] * 0.45359237 / 3600 # fuel flow left engine [kg/s] for CL-CD measurements
FFr_CLCD = CLCD[:,6] * 0.45359237 / 3600 # fuel flow right engine [kg/s] for CL-CD measurements
F_used_CLCD = CLCD[:,7] * 0.45359237 # fuel used by both engines in [kg] for CL-CD measurements
TAT_CLCD = CLCD[:,8] + 273.15 # Total ambient temperature [K] for CL-CD measurements

# ---- Trim measurements ----
t_trim = trim[:,0] * 60 + trim[:,1] # total time in [s] for trim measurements
hp_trim = trim[:,2] * 0.3048 # hp [m] for trim measurements
IAS_trim = trim[:,3] * 0.514444444 # indicated airspeed [m/s] for trim measurements
alpha_trim = np.deg2rad(trim[:,4]) # AoA in [rad] for trim measurements
FFl_trim = trim[:,8] * 0.45359237 / 3600 # fuel flow left engine [kg/s] for trim measurements
FFr_trim = trim[:,9] * 0.45359237 / 3600 # fuel flow right engine [kg/s] for trim measurements
F_used_trim = trim[:,10] * 0.45359237 # fuel used by both engines in [kg] for trim measurements
TAT_trim = trim[:,11] + 273.15 # Total ambient temperature [K] for trim measurements

de_trim = np.deg2rad(trim[:,5]) # elevator deflection angle in [rad] for trim measurements
detr_trim = np.deg2rad(trim[:,6]) # trim tab angle in [rad] for trim measurements
Fe_trim = trim[:,7] # Stick force in [N] for trim measurements

# ---- CG trim measurements ----
t_cg = cg[:,0] * 60 + cg[:,1] # total time in [s] for trim measurements
hp_cg = cg[:,2] * 0.3048 # hp [m] for trim measurements
IAS_cg = cg[:,3] * 0.514444444 # indicated airspeed [m/s] for trim measurements
alpha_cg = np.deg2rad(cg[:,4]) # AoA in [rad] for trim measurements
FFl_cg = cg[:,8] * 0.45359237 / 3600 # fuel flow left engine [kg/s] for trim measurements
FFr_cg = cg[:,9] * 0.45359237 / 3600 # fuel flow right engine [kg/s] for trim measurements
F_used_cg = cg[:,10] * 0.45359237 # fuel used by both engines in [kg] for trim measurements
TAT_cg = cg[:,11] + 273.15 # Total ambient temperature [K] for trim measurements

de_cg = np.deg2rad(cg[:,5]) # elevator deflection angle in [rad] for trim measurements
detr_cg = np.deg2rad(cg[:,6]) # trim tab angle in [rad] for trim measurements
Fe_cg = cg[:,7] # Stick force in [N] for trim measurements

# ---- CL and CD curves ----
OEW = 9165 * 0.45359237 # Operational empty weight [kg]
W_ramp = np.sum(mass_pl) + F_block + OEW # [kg]                       !!!!!!!!!!!!!!!!
rho_0 = 1.225 # [kg/m^3]
S = 30.0 # [m^2]
# CD_0 = 0.04 # zero-lift drag []
c = 2.0569 # mean chord [m]
b = 15.911 # Span [m]
A = b/c # Aspect ratio []
labda = -0.0065 # Temp gradient
T_0  = 288.15
p_0 = 101325
g = 9.81
R = 287
gamma = 1.4
nu_air = 13.28*10**-6 # Kinematic vicosity [m^2/s]


hpp_CLCD = np.reshape(hp_CLCD, (6,1))

p = p_0*((1+(labda*hp_CLCD/T_0))**(-g/(labda*R)))

# Mach number
M = np.sqrt((2/(gamma-1))*((1+(p_0/p)*(((1+(((gamma-1)/(2*gamma))*(rho_0/p_0)*IAS_CLCD*IAS_CLCD))**(gamma/(gamma-1)))-1))**((gamma-1)/gamma)-1))
T = np.reshape(TAT_CLCD/(1+((gamma-1)/2)*M*M),(6,1))
a = np.sqrt(gamma*R*T)

MM = np.reshape(M,(6,1))
V_t = MM*a

pp = np.reshape(p, (6,1))
rho = pp/(R*T)

V_e = V_t * np.sqrt(rho/rho_0)

TISA = np.reshape(T_0 + labda* hp_CLCD,(6,1))

DTT = TISA - T

Re = c * V_t / nu_air # Reynolds range:  10E6 - 21E6

FFll_CLCD = np.reshape(FFl_CLCD, (6,1))
FFrr_CLCD = np.reshape(FFr_CLCD, (6,1))
np.savetxt("matlab.dat",np.hstack((hpp_CLCD, MM, DTT, FFll_CLCD, FFrr_CLCD)), delimiter = ' ')
subprocess.run("thrust(1).exe")
Thrust = np.sum(np.loadtxt('thrust.dat'),1)
Thrustt = np.reshape(Thrust,(6,1))
# CL = (W_ramp-F_used_CLCD)*2*9.81/(rho_0*(np.reshape(V_e,(6,))**2)*S)
CL = (W_ramp-F_used_CLCD)*2*9.81/(rho_0*(IAS_CLCD**2)*S)
# CD1 = Thrust * 2 / (rho_0 * IAS_CLCD**2 * S)
CD1 = np.reshape(Thrustt * 2 / (rho * V_t**2 * S),(6,))
slope = np.mean(((CL**2)[1:]-(CL**2)[:-1])/(CD1[1:]-CD1[:-1]))
e = np.pi*A/slope
CD_0 = np.mean(CD1-CL**2/(np.pi*A*e))
CD = CD_0 + CL**2/(A*e*np.pi)
# CL[0]/slope
# CD_0 = 0.0203
# e = CL*CL/(np.pi*A*(CD1-CD_0))
CL_alpha = (CL[-1]-CL[0])/(alpha_CLCD[-1]-alpha_CLCD[0]) # [rad]

# ----- Trim curve measurements -----
# cg locations for cg shift measurements. Values taken from mass balance report by Wessel
x_cg1 = 7.1489 # [m]
x_cg2 = 7.0952 # [m]

hpp_cg = np.reshape(hp_cg, (2,1))
p_cg = p_0*((1+(labda*hp_cg/T_0))**(-g/(labda*R)))
# Mach number
M_cg = np.sqrt((2/(gamma-1))*((1+(p_0/p_cg)*(((1+(((gamma-1)/(2*gamma))*(rho_0/p_0)*IAS_cg*IAS_cg))**(gamma/(gamma-1)))-1))**((gamma-1)/gamma)-1)) 
T_cg = np.reshape(TAT_cg/(1+((gamma-1)/2)*M_cg*M_cg),(2,1))
a_cg = np.sqrt(gamma*R*T_cg)
MM_cg = np.reshape(M_cg,(2,1))
V_t_cg = MM_cg*a_cg
pp_cg = np.reshape(p_cg, (2,1))
rho_cg = pp_cg/(R*T_cg)

V_e_cg = V_t_cg * np.sqrt(rho_cg/rho_0)

CN_cg = np.mean((W_ramp-F_used_cg)*2*9.81/(rho_0*(np.reshape(V_e_cg,(2,))**2)*S))
Cm_delta = -1/(de_cg[-1]-de_cg[0])*CN_cg*(x_cg2-x_cg1)/c

# Trim: find d(delta)/d(alpha)
hpp_trim = np.reshape(hp_trim, (7,1))
p_trim = p_0*((1+(labda*hp_trim/T_0))**(-g/(labda*R)))
# Mach number
M_trim = np.sqrt((2/(gamma-1))*((1+(p_0/p_trim)*(((1+(((gamma-1)/(2*gamma))*(rho_0/p_0)*IAS_trim*IAS_trim))**(gamma/(gamma-1)))-1))**((gamma-1)/gamma)-1)) 
T_trim = np.reshape(TAT_trim/(1+((gamma-1)/2)*M_trim*M_trim),(7,1))
a_trim = np.sqrt(gamma*R*T_trim)
MM_trim = np.reshape(M_trim,(7,1))
V_t_trim = MM_trim*a_trim
pp_trim = np.reshape(p_trim, (7,1))
rho_trim = pp_trim/(R*T_trim)
V_e_trim = V_t_trim * np.sqrt(rho_trim/rho_0)
TISA_trim = np.reshape(T_0 + labda* hp_trim,(7,1))
DTT_trim = TISA_trim - T_trim

CN_trim = (W_ramp-F_used_trim)*2*9.81/(rho_0*(np.reshape(V_e_trim,(7,))**2)*S)
deda = (np.min(de_trim)-np.max(de_trim))/(np.max(alpha_trim)-np.min(alpha_trim)) # d(delta)/d(alpha)
Cm_a = -deda*Cm_delta

# ----- Reduction of measurement data -------
Ws = 60500 # [N] Standard aircraft mass
V_epin = np.reshape(np.reshape(V_e_trim,(7,)) * np.sqrt(Ws/((W_ramp-F_used_trim)*9.81)),(7,1))

# TCs computation
mfs = np.ones((7,1))*0.048 # [kg/s] standard fuel flow
d_eng = 0.69 # Characteristic diameter engine [m]
CmTC = -0.0064 # dCm/dTC []
FFll_trim = np.reshape(FFl_trim, (7,1))
FFrr_trim = np.reshape(FFr_trim, (7,1))
np.savetxt("matlab.dat",np.hstack((hpp_trim, MM_trim, DTT_trim, FFll_trim, FFrr_trim)), delimiter = ' ')
subprocess.run("thrust(1).exe")
Thrust_trim = np.sum(np.loadtxt('thrust.dat'),1)
Thrustt_trim = np.reshape(Thrust_trim,(7,1))
TC = Thrustt_trim/(0.5*rho_trim*V_t_trim**2*d_eng**2)

np.savetxt("matlab.dat",np.hstack((hpp_trim, MM_trim, DTT_trim, mfs, mfs)), delimiter = ' ')
subprocess.run("thrust(1).exe")
Ts_trim = np.sum(np.loadtxt('thrust.dat'),1)
Tss_trim = np.reshape(Ts_trim,(7,1))
TCs = Tss_trim/(0.5*rho_trim*V_t_trim**2*d_eng**2)

dee_trim = np.reshape(de_trim,(7,1))
W_trim = np.reshape((W_ramp-F_used_trim)*9.81,(7,1))
Cm0 = np.mean(-dee_trim*Cm_delta-CmTC*TC-(Cm_a*W_trim/(CL_alpha*0.5*rho_0*V_e_trim**2*S)))

# delta*
de_red = -1/Cm_delta*(Cm0+Cm_a*W_trim/(CL_alpha*0.5*rho_0*V_epin**2*S)+CmTC*TCs)
# Computation of stick force F_e*
Fe_red = np.reshape(Fe_trim,(7,1))*Ws/W_trim

plt.figure()
plt.subplot(1,2,1)
plt.title("Reduced elevator deflection vs reduced equivalent airspeed")
plt.grid()
plt.xlabel('V_e reduced [m/s]')
plt.ylabel('delta_e reduced [deg]')
plt.plot(np.sort(V_epin,0), np.sort(np.rad2deg(de_red),0))
plt.subplot(1,2,2)
plt.title("Reduced stick force vs reduced equivalent airspeed")
plt.grid()
plt.gca().invert_yaxis()
plt.xlabel('V_e reduced [m/s]')
plt.ylabel('F_e reduced [N]')
plt.plot(np.sort(V_epin,0), np.sort(Fe_red,0))
plt.show()

print("Oswald factor:", e)
print("CD_0:", CD_0)
print("CL_alpha:", CL_alpha)
print("Cm_alpha:", Cm_a)
print("Cm_delta:", Cm_delta)
print("Cm_0:", Cm0)

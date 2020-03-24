# Citation 550 - Linear simulation

#from Verification import geteigenvalues
import numpy as np
from math import pi, sin, cos
from matplotlib import pyplot as plt
from Verification import get_short_period_eigenvalues,get_phugoid_eigenvalues,get_aperiodic_roll_eigenvalue,get_spiral_motion_eigenvalue,get_dutch_roll_with_aperiodic_roll_eigenvalues,get_eigenvalues

# xcg = 0.25 * c

# Stationary flight condition

hp0    = 1500      	     # pressure altitude in the stationary flight condition [m]
V0     = 82.3            # true airspeed in the stationary flight condition [m/sec]
alpha0 =  0.04           # angle of attack in the stationary flight condition [rad]
th0    =  0.08           # pitch angle in the stationary flight condition [rad]

# Aircraft mass
m      = 10000            # mass [kg]

#Definition returns the state space matrix for both the symmetric as asymmetric case
def statespacematrix(hp0,V0,alpha0,th0,m):

    # aerodynamic properties (from table C)
    e      = 0.9502171503627244             # Oswald factor [ ]
    CD0    = 0.021055067079004692           # Zero lift drag coefficient [ ]
    CLa    = 4.492916965171044              # Slope of CL-alpha curve [ ]

    # Longitudinal stability (from table C)
    Cma    = -0.6824789003603976            # longitudinal stability [ ]
    Cmde   = -1.5450563993322224            # elevator effectiveness [ ]


    # Aircraft geometry

    S      = 30.00	            # wing area [m^2]
    Sh     = 0.2 * S            # stabiliser area [m^2]
    Sh_S   = Sh / S	            # [ ]
    lh     = 0.71 * 5.968       # tail length [m]
    c      = 2.0569	            # mean aerodynamic cord [m]
    lh_c   = lh / c	            # [ ]
    b      = 15.911	            # wing span [m]
    bh     = 5.791	            # stabiliser span [m]
    A      = b ** 2 / S         # wing aspect ratio [ ]
    Ah     = bh ** 2 / Sh       # stabiliser aspect ratio [ ]
    Vh_V   = 1	                # [ ]
    ih     = -2 * pi / 180      # stabiliser angle of incidence [rad]

    # Constant values concerning atmosphere and gravity

    rho0   = 1.2250          # air density at sea level [kg/m^3]
    lambdaa = -0.0065        # temperature gradient in ISA [K/m]
    Temp0  = 288.15          # temperature at sea level in ISA [K]
    R      = 287.05          # specific gas constant [m^2/sec^2K]
    g      = 9.81            # [m/sec^2] (gravity constant)

    # air density [kg/m^3]
    rho    = rho0 * np.power( ((1+(lambdaa * hp0 / Temp0))), (-((g / (lambdaa*R)) + 1)))
    W      = m * g            # [N]       (aircraft weight)

    # Constant values concerning aircraft inertia

    muc    = m / (rho * S * c)
    mub    = m / (rho * S * b)
    KX2    = 0.019
    KZ2    = 0.042
    KXZ    = 0.002
    KY2    = 1.25 * 1.114

    # Aerodynamic constants

    Cmac   = 0                      # Moment coefficient about the aerodynamic centre [ ]
    CNwa   = CLa                    # Wing normal force slope [ ]
    CNha   = 2 * pi * Ah / (Ah + 2) # Stabiliser normal force slope [ ]
    depsda = 4 / (A + 2)            # Downwash gradient [ ]

    # Lift and drag coefficient

    CL = 2 * W / (rho * V0 ** 2 * S)              # Lift coefficient [ ]
    CD = CD0 + (CLa * alpha0) ** 2 / (pi * A * e) # Drag coefficient [ ]

    # Stability derivatives

    CX0    = W * sin(th0) / (0.5 * rho * V0 ** 2 * S)
    CXu    = -0.095
    #CXu = -0.125            #Change in stability coefficient for phugoid 1/2
    CXa    = 0.47966		# Positive! (has been erroneously negative since 1993)
    CXadot = +0.08330
    CXq    = -0.28170
    CXde   = -0.03728

    CZ0    = -W * cos(th0) / (0.5 * rho * V0 ** 2 * S)
    CZu    = -0.37616
    #CZu = -0.63              #Change in stability coefficient for phugoid 2/2
    CZa    = -5.74340
    #CZa = -5.0                #Change in stability coefficient for short period 1/2
    CZadot = -0.00350
    CZq    = -5.66290
    CZde   = -0.69612

    Cmu    = +0.06990
    Cmadot = +0.17800
    Cmq    = -8.79415
    #Cmq = -12.3                 #Change in stability coefficient for short period 2/2

    CYb    = -0.7500
    CYbdot =  0
    CYp    = -0.0304
    CYr    = +0.8495
    CYda   = -0.0400
    CYdr   = +0.2300

    Clb    = -0.10260
    #Clb = -0.26                #Change in stability coefficient for spiral motion 1/2
    Clp    = -0.71085
    #Clp = -0.6              #Change in stability coefficient for aperiod roll
    Clr    = +0.23760
    Clda   = -0.23088
    Cldr   = +0.03440

    Cnb    =  +0.1348
    Cnbdot =   0
    Cnp    =  -0.0602
    Cnr    =  -0.2061
    #Cnr = -0.3                 #Change in stability coefficient for spiral motion 2/2
    Cnda   =  -0.0120
    Cndr   =  -0.0939

    '''
    # ----------------------------------Eigenvalue Verification ---------------------------------------
    # print("simplified short period",
    sp = get_short_period_eigenvalues(CZa, CZadot, muc, CZq, Cma, Cmadot, Cmq, KY2)

    # print("simplified phugoid",d
    phug = get_phugoid_eigenvalues(CXu, CZu, muc, CZ0)

    # print("simplified aperiodic roll",
    aperiodic = [get_aperiodic_roll_eigenvalue(Clp, mub, KX2)]
    # print("simplified dutch roll",get_dutch_roll_eigenvalues(mu_b, kz, cn_r, cy_beta, cn_beta))
    # print("simplified spiral",
    spiral = [get_spiral_motion_eigenvalue(CL, Clb, Cnb, Clr, Cnr, Clp, CYb, Cnp, mub)]
    # print("simplified Dutch roll",

    dutch_roll = get_dutch_roll_with_aperiodic_roll_eigenvalues(mub, KX2, KZ2, KXZ, Clr, Cnp, Cnr, Clp, Clb, Cnb)
    print(sp)
    print("sp: ", sp[0]*(V0/c))
    print("phug:", phug[0]*V0/c)
    print("aperiodic: ,",aperiodic[0]*V0/b)
    print("spiral: ", spiral[0]*V0/b)
    print("dr: ", dutch_roll[0]*V0/b, dutch_roll[1]*V0/b, dutch_roll[2]*V0/b)

    if plot:
        plt.figure(figsize=(10, 5))
        plt.subplot(121)
        plt.scatter([x.real*V0/c for x in sp], [x.imag*V0/c for x in sp], color='r', marker='x')
        plt.scatter([x.real*V0/c for x in phug], [x.imag*V0/c for x in phug], color='b', marker='x')
        plt.scatter([x.real*V0/b for x in aperiodic], [x.imag*V0/b for x in aperiodic], color='y', marker='x')
        plt.scatter([x.real*V0/b for x in spiral], [x.imag*V0/b for x in spiral], color='g', marker='x')
        plt.scatter([x.real*V0/b for x in dutch_roll], [x.imag*V0/b for x in dutch_roll], color='c', marker='x')
        plt.title("Eigenvalues of simplified models")
        plt.grid()
        plt.legend(labels=['short period', 'phugoid', 'aperiodic roll', 'spiral', '"special" dutch roll '])
        plt.ylabel("imaginary")
        plt.xlabel("real")
        #plt.show()
    '''
    # --------Symmetric equations of motion in the form of:  C1 * xdot + C2 * x + C3 * u     --------

    C1 = np.array([[(-2*muc*c)/(V0**2),0                    ,0    ,0],
                   [0                 ,(CZadot-2*muc)*(c/V0),0    ,0],
                   [0                 ,0                    ,-c/V0,0],
                   [0                 ,(Cmadot*c)/V0        ,0    ,(-2*muc*c**2 * KY2)/(V0**2)]])

    C2 = np.array([[CXu/V0,CXa,CZ0 ,CXq*(c/V0)],
                   [CZu/V0,CZa,-CX0,(CZq+2*muc)*(c/V0)],
                   [0     ,0  ,0   ,c/V0],
                   [Cmu/V0,Cma,0   ,(Cmq*c)/V0]])

    C3 = np.array([[CXde],
                   [CZde],
                   [0],
                   [Cmde]])

    ### State Space Symmetric for state vector [u,alpha,theta,q]
    A_sym = -np.dot(np.linalg.inv(C1),C2)
    B_sym = -np.dot(np.linalg.inv(C1),C3)

    #---------Asymmetric equations of motion in the form of: D1 * xdot + D2 * x + D3 * u ----------------

    D1 = np.array([[(CYbdot-2*mub)*(b/V0),0,0,0],
                   [0,-b/(2*V0),0,0],
                   [0,0,-2*mub*KX2*(b**2 / V0**2),2*mub*KXZ*(b**2 / V0**2)],
                   [Cnbdot*(b/V0),0,2*mub*KXZ*(b**2/V0**2),-2*mub*KZ2*(b**2 / V0**2)]])

    D2 = np.array([[CYb,CL,(CYp*b)/(2*V0),(CYr-4*mub)*(b/(2*V0))],
                   [0,0,b/(2*V0),0],
                   [Clb,0,(Clp*b)/(2*V0),Clr*(b/(2*V0))],
                   [Cnb,0,(Cnp*b)/(2*V0),(Cnr*b)/(2*V0)]])

    D3 = np.array([[CYda,CYdr],
                   [0,0],
                   [Clda,Cldr],
                   [Cnda,Cndr]])

    ### State Space Asymmetric for state vector [beta, phi, p, r]

    A_asym = -np.dot(np.linalg.inv(D1),D2)
    B_asym = -np.dot(np.linalg.inv(D1),D3)

    return A_sym,B_sym,A_asym,B_asym

'''
plot = True
if plot:
    a = statespacematrix(hp0,V0,alpha0,th0,m)

    eigenvalues_sym = get_eigenvalues(a[0])
    eigenvalues_asym = get_eigenvalues(a[2])

    plt.subplot(122)
    plt.scatter([x.real for x in eigenvalues_sym[0:2]], [x.imag for x in eigenvalues_sym[0:2]], color='r', marker='x')
    plt.scatter([x.real for x in eigenvalues_sym[2:4]], [x.imag for x in eigenvalues_sym[2:4]], color='b', marker='x')
    plt.scatter([x.real for x in eigenvalues_asym[0:1]], [x.imag for x in eigenvalues_asym[0:1]], color='y', marker='x')
    plt.scatter([x.real for x in eigenvalues_asym[3:4]], [x.imag for x in eigenvalues_asym[3:4]], color='g', marker='x')
    plt.scatter([x.real for x in eigenvalues_asym[1:3]], [x.imag for x in eigenvalues_asym[1:3]], color='c', marker='x')
    plt.title("Eigenvalues of total models")
    plt.grid()
    plt.legend(labels=['short period', 'phugoid', 'aperiodic roll', 'spiral', 'dutch roll '])
    plt.ylabel("imaginary")
    plt.xlabel("real")
    plt.show()


A_ = 4*muc**2*KY2*(CZadot-2*muc)
B_ = Cmadot*2*muc*(CZq+2*muc) - Cmq*2*muc*(CZadot-2*muc)- 2*muc*KY2*(CXu*(CZadot-2*muc)-2*muc*CZa)
C_ = Cma*2*muc*(CZq+2*muc)-Cmadot*(2*muc*CX0+CXu*(CZq+2*muc))+Cmq*(CXu*(CZadot-2*muc)-2*muc*CZa)+2*muc*KY2*(CXa*CZu-CZa*CXu)
D_ = Cmu*(CXa*(CZq+2*muc)-CZ0*(CZadot-2*muc))-Cma*(2*muc*CX0+CXu*(CZq+2*muc))+Cmadot*(CX0*CXu-CZ0*CZu)+Cmq*(CXu*CZa-CZu*CXa)
E_ = -Cmu*(CX0*CXa+CZ0*CZa)+Cma*(CX0*CXu+CZ0*CZu)
R_ = B_*C_*D_-A_*D_**2-B_**2*E_

# print(A_, B_, C_, D_, E_, R_)


# print("simplified short period",get_short_period_eigenvalues(cz_alpha, cza_dot, mu_c, cz_q, cm_alpha, cm_alpha_dot, cm_q, ky))
# print("simplified phugoid",get_phugoid_eigenvalues(cx_u, cz_u, mu_c, cz_0, ky))
# print("simplified aperiodic roll",get_aperiodic_roll_eigenvalue(cl_p, mu_b, kx))
# print("simplified dutch roll",get_dutch_roll_eigenvalues(mu_b, kz, cn_r, cy_beta, cn_beta))
# print("simplified spiral",get_spiral_motion_eigenvalue(cl, cl_beta, cn_beta, cl_r, cn_r, cl_p, cy_beta, cn_p, mu_b))
# print("simplified Dutch roll",get_dutch_roll_with_aperiodic_roll_eigenvalues(mu_b, kx, kz, kxz, cl_r, cn_p, cn_r, cl_p, cl_beta, cn_beta))



'''
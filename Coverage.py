import numpy as np
import control as ctrl
import scipy.io as sio
import matplotlib.pyplot as plt
from math import pi, sin, cos
import cmath

from Cit_par import statespacematrix

statespacematrix(hp0,V0,alpha0,th0,m)


from IAS_to_CAS import IAS_to_CAS

IAS_to_CAS(V_IAS)


from Verification import get_eigenvalues, get_coefficients_3x3_matrix,abcd_formula, abc_formula, get_short_period_eigenvalues,get_phugoid_eigenvalues,get_aperiodic_roll_eigenvalue,get_dutch_roll_eigenvalues, get_spiral_motion_eigenvalue, get_dutch_roll_with_aperiodic_roll_eigenvalues

get_eigenvalues(A)
get_coefficients_3x3_matrix(B)
abcd_formula(a, b, c, d)
abc_formula(a, b, c)
get_short_period_eigenvalues(cz_alpha, cza_dot, mu_c, cz_q, cm_alpha, cm_alpha_dot, cm_q, ky)
get_phugoid_eigenvalues(cx_u, cz_u, mu_c, cz_0)
get_aperiodic_roll_eigenvalue(cl_p, mu_b, kx)
get_dutch_roll_eigenvalues(mu_b, kz, cn_r, cy_beta, cn_beta)
get_spiral_motion_eigenvalue(cl, cl_beta, cn_beta, cl_r, cn_r, cl_p, cy_beta, cn_p, mu_b)
get_dutch_roll_with_aperiodic_roll_eigenvalues(mu_b, kx, kz, kxz, cl_r, cn_p, cn_r, cl_p, cl_beta, cn_beta)

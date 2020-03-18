import numpy as np
import cmath


def get_eigenvalues(A):

    eigenvalues = np.linalg.eigvals(A)

    return eigenvalues


def get_coefficients_3x3_matrix(B):
    # det = lambda x: (B[0, 0] - x)*((B[1, 1] - x)*(B[2, 2] - x) - B[2, 1]*B[1, 2]) - (B[0, 1])*(B[1, 0]*(B[2, 2] - x)
    # - B[1, 2]*B[1, 0]) + (B[0, 2])*(B[1, 0]*B[2, 1] - (B[1, 1] - x)*B[2, 0])
    a = -1
    b = B[0, 0] + B[1, 1] + B[2, 2]
    c = - B[0, 0] * (B[1, 1] + B[2, 2]) - B[1, 1]*B[2, 2] + B[1, 2]*B[2, 1] + B[0, 1]*B[1, 0] + B[0, 2]*B[2, 0]
    d = B[0, 0]*B[1, 1]*B[2, 2] - B[0, 0]*B[1, 2]*B[2, 1] - B[0, 1]*B[1, 0]*B[2, 2] + B[0, 1]*B[1, 2]*B[2, 0] + B[0, 2]*B[1, 0]*B[2, 1] - B[0, 2]*B[1, 1]*B[2, 0]
    return a, b, c, d


def abcd_formula(a, b, c, d):
    # this function solves for x in the following equation
    # a*x^3 + b*x^2 + c*x + d = 0
    # complex solutions are also given

    # first divide everything by a (a is not zero)
    bb = b/a
    cc = c/a
    dd = d/a
    # equation is now: x^3 + bb * x^2 + cc * x + dd = 0
    # convert equation into:
    # y^3 + py + q = 0, with p and q as follows:
    p = cc - bb**2/3
    q = 2/27*bb**3 - 1/3*bb*cc+dd
    ##
    # solve 27 z^2 + 27 qz - p^3 = 0 to get z
    z = abc_formula(27, 27*q, -p**3)
    print(z)
    # get r and phi from general solution
    r, phi = cmath.polar(z[0])
    # convert back to m using a general formula
    m1 = r**(1/3)*(np.cos(phi/3) - 1j*np.sin(phi/3))
    m2 = r**(1/3)*(np.cos((2*np.pi + phi)/3) + 1j*np.sin((2*np.pi + phi)/3))
    m3 = r**(1/3)*(np.cos((4*np.pi + phi)/3) + 1j*np.sin((4*np.pi + phi)/3))
    print("m1: ", m1)
    print("m2: ", m2)
    print("m3: ", m3)
    n1 = -p/(3*m1)
    n2 = -p/(3*m2)
    n3 = -p/(3*m3)
    y1 = n1 + m1
    y2 = n2 + m2
    y3 = n3 + m3
    x1 = y1 - bb/3
    x2 = y2 - bb/3
    x3 = y3 - bb/3

    print(x1, x2, x3)
    return x1, x2, x3


def abc_formula(a, b, c):
    D = b**2 - 4*a*c
    sol = []

    if D == 0:
        sol.append(-b / (2*a))
    elif D < 0:
        sol.append((-b + cmath.sqrt(D))/(2*a))
        sol.append((-b - cmath.sqrt(D))/(2*a))
    else:
        sol.append((-b + np.sqrt(D))/(2*a))
        sol.append((-b - np.sqrt(D))/(2*a))

    return sol


def get_short_period_eigenvalues(cz_alpha, cza_dot, mu_c, cz_q, cm_alpha, cm_alpha_dot, cm_q, ky):
    a = 2*mu_c*ky*(2*mu_c - cza_dot)
    b = -2*mu_c*ky*cz_alpha - (2*mu_c + cz_q)*cm_alpha_dot - (2*mu_c-cza_dot)*cm_q
    c = cz_alpha*cm_q - (2*mu_c + cz_q)*cm_alpha
    eig_values = abc_formula(a, b, c)
    return eig_values


def get_phugoid_eigenvalues(cx_u, cz_u, mu_c, cz_0):
    a = -4*mu_c**2
    b = 2*mu_c*cx_u
    c = -cz_u*cz_0
    eig_values = abc_formula(a, b, c)
    return eig_values


def get_aperiodic_roll_eigenvalue(cl_p, mu_b, kx):
    # kx is assumed to be squared already
    lambda_1 = cl_p/(4*mu_b*kx)
    return lambda_1


def get_dutch_roll_eigenvalues(mu_b, kz, cn_r, cy_beta, cn_beta):
    # kz is assumed to be squared already
    a = 8*mu_b**2*kz
    b = -2*mu_b*(cn_r + 2*kz*cy_beta)
    c = 4*mu_b*cn_beta + cy_beta*cn_r
    eig_values = abc_formula(a, b, c)
    return eig_values


def get_spiral_motion_eigenvalue(cl, cl_beta, cn_beta, cl_r, cn_r, cl_p, cy_beta, cn_p, mu_b):
    ##
    lambda_2 = (2*cl*(cl_beta*cn_r - cn_beta*cl_r))/(cl_p*(cy_beta*cn_r + 4*mu_b*cn_beta) - cn_p*(cy_beta*cl_r + 4*mu_b*cl_beta))
    return lambda_2


def get_dutch_roll_with_aperiodic_roll_eigenvalues(mu_b, kx, kz, kxz, cl_r, cn_p, cn_r, cl_p, cl_beta, cn_beta):
    ##
    a = 4*mu_b**2*(kx*kz - kxz**2)
    b = -mu_b*((cl_r + cn_p)*kxz + cn_r*kx + cl_p*kz)
    c = 2*mu_b*(cl_beta*kxz + cn_beta*kx) + 0.25*(cl_p*cn_r - cn_p*cl_r)
    d = 0.5*(cl_beta*cn_p - cn_beta*cl_p)
    eig_values = abcd_formula(a, b, c, d)
    return eig_values

# print("simplified short period",get_short_period_eigenvalues(cz_alpha, cza_dot, mu_c, cz_q, cm_alpha, cm_alpha_dot, cm_q, ky))
# print("simplified phugoid",get_phugoid_eigenvalues(cx_u, cz_u, mu_c, cz_0, ky))
# print("simplified aperiodic roll",get_aperiodic_roll_eigenvalue(cl_p, mu_b, kx))
# print("simplified dutch roll",get_dutch_roll_eigenvalues(mu_b, kz, cn_r, cy_beta, cn_beta))
# print("simplified spiral",get_spiral_motion_eigenvalue(cl, cl_beta, cn_beta, cl_r, cn_r, cl_p, cy_beta, cn_p, mu_b))
# print("simplified Dutch roll",get_dutch_roll_with_aperiodic_roll_eigenvalues(mu_b, kx, kz, kxz, cl_r, cn_p, cn_r, cl_p, cl_beta, cn_beta))


"""
A = np.matrix([[1, 2, 3], [5, 6, 7], [2, -1, 3]])
# A = np.matrix([[1, 0, 0], [0, 2, 3], [0, 6, 5]])

print("using numpy function: ", get_eigenvalues(A))
print("using own function: ", abcd_formula(get_coefficients_3x3_matrix(A)[0], get_coefficients_3x3_matrix(A)[1], get_coefficients_3x3_matrix(A)[2], get_coefficients_3x3_matrix(A)[3]))
"""
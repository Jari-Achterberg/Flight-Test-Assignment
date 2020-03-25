from matplotlib import pyplot as plt

# NUMERICAL MODEL
sp = [-1.30959588+2.14008295j, -1.30959588-2.14008295j]
phugoid = [-0.00473832+0.13011474j, -0.00473832-0.13011474j]
spiral = [0.01096293+0.j]
ap = [-4.24164388+0.j]
dr = [-0.28474598+1.99940621j, -0.28474598-1.99940621j]

# SIMPLIFIED NUM MODEL
phug_simple = [(-0.010170937122260445-0.0939576268367825j), (-0.010170937122260445+0.0939576268367825j)]
spiral_simple = [0.011294029988588404]
dr_simple = [(-0.22134045943061767+1.9658058541899388j), (-4.003140853704277-2.7798466841257503e-16j), (-0.22134045943061767-1.9658058541899388j)] # check
ap_simple = [-4.147995610820421]  # checked
sp_simple = [(-1.3045855138173366+2.1388159093429007j), (-1.3045855138173366-2.1388159093429007j)]  # checked

# validation
phug_val = [-0.005635342931381669 + 0.1595121936323835j, -0.005635342931381669 - 0.1595121936323835j]
spiral_val = [0.009200254586673019]  # [0.00967002205022245]
dr_val = [-0.23028145533553002 + 1.9482745138541355j, -0.23028145533553002 - 1.9482745138541355j]
ap_val = [-4.332168]
sp_val = [-1.3004637533957697 + 2.001014429038085j, -1.3004637533957697 - 2.001014429038085j]


def plot_eigenvalues(sp, phugoid, spiral, ap, dr, phug_simple, spiral_simple, dr_simple, ap_simple, sp_simple, phug_val, spiral_val, dr_val, ap_val, sp_val):
    plt.subplot(131)
    plt.scatter([x.real for x in sp], [x.imag for x in sp], color='r', marker='x')
    plt.scatter([x.real for x in phugoid], [x.imag for x in phugoid], color='b', marker='x')
    plt.scatter([x.real for x in ap], [x.imag for x in ap], color='y', marker='x')
    plt.scatter([x.real for x in spiral], [x.imag for x in spiral], color='g', marker='x')
    plt.scatter([x.real for x in dr], [x.imag for x in dr], color='c', marker='x')
    plt.title("Eigenvalues of total models")
    plt.grid()
    plt.legend(labels=['short period', 'phugoid', 'aperiodic roll', 'spiral', 'dutch roll '])
    plt.ylabel("imaginary")
    plt.xlabel("real")

    plt.subplot(132)
    plt.scatter([x.real for x in sp_simple], [x.imag for x in sp_simple], color='r', marker='x')
    plt.scatter([x.real for x in phug_simple], [x.imag for x in phug_simple], color='b', marker='x')
    plt.scatter([x.real for x in ap_simple], [x.imag for x in ap_simple], color='y', marker='x')
    plt.scatter([x.real for x in spiral_simple], [x.imag for x in spiral_simple], color='g', marker='x')
    plt.scatter([x.real for x in dr_simple], [x.imag for x in dr_simple], color='c', marker='x')
    plt.title("Eigenvalues of SIMPLIFIED models")
    plt.grid()
    plt.legend(labels=['short period', 'phugoid', 'aperiodic roll', 'spiral', 'dutch roll '])
    plt.ylabel("imaginary")
    plt.xlabel("real")

    plt.subplot(133)
    plt.scatter([x.real for x in sp_val], [x.imag for x in sp_val], color='r', marker='x')
    plt.scatter([x.real for x in phug_val], [x.imag for x in phug_val], color='b', marker='x')
    plt.scatter([x.real for x in ap_val], [x.imag for x in ap_val], color='y', marker='x')
    plt.scatter([x.real for x in spiral_val], [x.imag for x in spiral_val], color='g', marker='x')
    plt.scatter([x.real for x in dr_val], [x.imag for x in dr_val], color='c', marker='x')
    plt.title("Validation Eigenvalues using plots")
    plt.grid()
    plt.legend(labels=['short period', 'phugoid', 'aperiodic roll', 'spiral', 'dutch roll '])
    plt.ylabel("imaginary")
    plt.xlabel("real")
    plt.ylim(-2.3, 2.3)
    plt.show()


plot_eigenvalues(sp, phugoid, spiral, ap, dr, phug_simple, spiral_simple, dr_simple, ap_simple, sp_simple, phug_val, spiral_val, dr_val, ap_val, sp_val)

import numpy as np
import scipy.sparse as sps
import matplotlib.pyplot as plt


def finite(n): 
    # Define the parameters
    u_b = 10
    h_0 = 60e-6
    l = 0.1
    eta_0 = 0.1

    s_h = h_0/(np.exp(-1))- h_0

    # Define the grid
    x = np.linspace(0, l, n)
    dx = abs(x[0] - x[1])

    h = (h_0 + s_h)*np.exp(-x/l)


    
    # Create the diagonals
    south = np.zeros(len(h))
    for i in range(len(h)-2):
        i = i+1
        south[i-1] = 3 * h[i]**2 * (h[i+1] - h[i-1]) / (2 * dx) * -1/(2*dx) + h[i]**3 * 1/(dx**2)

    # Create the diagonals
    middle = np.zeros(len(h))
    for i in range(len(h)):
        middle[i] = h[i]**3 * -2/(dx**2)

    # Bondary conditions
    middle[0] = 1
    middle[-1] = 1

    # Create the diagonals
    north = np.zeros(len(h))
    for i in range(len(h)-2):
        i = i+1
        north[i+1] = 3 * h[i]**2 * (h[i+1] - h[i-1]) / (2 * dx) * 1/(2*dx) + h[i]**3 * 1/(dx**2)


    # Create the matrix
    data = np.array([north, middle, south])
    A = sps.spdiags(data, [1, 0, -1], n, n, format = 'csr').T


    # rhs
    rhs = np.zeros(len(h))
    for i in range(len(h)-2):
        i=i+1
        rhs[i] = 6 * eta_0 * u_b * (h[i+1] - h[i-1]) / (2 * dx)

    # Solve the system
    p = sps.linalg.spsolve(A, rhs)

    return x, p


# Plot the solution
N = 2**np.arange(17)+2

x = []
p = []

for i in range(len(N)):
    X, P = finite(N[i])
    x.append(X)
    p.append(P)


plt.figure()
for i in range(len(N)):
    plt.plot(x[i], p[i], label='N = %d' % N[i])
plt.legend()
plt.xlabel('x')
plt.ylabel('p')
plt.title('Finite difference solution')
plt.grid()
plt.show()


# Calculate the error
baseline = 450668.4416
error = []

for i in range(len(N)):
    error.append(abs(np.trapezoid(p[i], x[i]) - baseline))

plt.figure()
for i in range(len(N)):
    plt.loglog(N, error, label='N = %d' % N[i])
plt.legend()
plt.xlabel('N')
plt.ylabel('Error')
plt.title('Error')
plt.grid()
plt.show()


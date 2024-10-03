import numpy as np
import scipy.sparse as sps

# Define the parameters
ub = 10
h_0 = 0.000001
l = 0.1
eta_0 = 0.1

s_h = 0.0001

# Define the grid
n = 5
x = np.linspace(0, l, n)
dx = abs(x[0] - x[1])

h = (h_0 + s_h)*np.exp(-x/l)


# Create the diagonals
south = np.zeros(len(h))
for i in range(len(h)-2):
    i = i+1
    south[i] = 3 * h[i]**2 * (h[i+1] - h[i-1]) / (2 * dx) * -1/(2*dx) + h[i]**3 * 1/(dx**2)

# Bondary conditions
south[-1] = 0

# Create the diagonals
middle = np.zeros(len(h))
for i in range(len(h)):
    middle[i] = h[i]**3 * -2/(dx**2)

# Bondary conditions
middle[0] = 1
middle[-1] = 1


# Create the diagonals
north = np.zeros(len(h)-1)

import numpy as np
import scipy.sparse as sps
import matplotlib.pyplot as plt

def finit_2d(nx: int = 50, ny: int = 50, l: float = 0.1, b:float = 0.1, u_b: float = 10, h_0: float = 60e-6, eta_0: float = 0.1, type: int = 0):
    '''
    Finite difference method to solve the Reynolds equation for a thrust bearing with linear inclined pads
    type 0 = linear inclined pads
    type 1 = step pads
    '''

    s_h = np.sqrt(2) * h_0

    X = np.linspace(0, l, nx)
    Y = np.linspace(0, b, ny)
    dx = abs(X[0] - X[1])
    dy = abs(Y[0] - Y[1])


    h = np.zeros((ny, nx))

    if type == 0: # Linear incline
        for i in range(nx):
            h[i,:] = h_0 + s_h * (1 - X/l)
    elif type == 1: # Step pads
        for i in range(nx):
            if Y[i] < l*0.7182:
                h[:,i] = h_0 + s_h
            else:
                h[:,i] = h_0

    M = sps.eye(nx * ny)
    M = M.tocsr()

    rhs = np.zeros(nx * ny)

    for i in range(1, nx-1):
        for j in range(1, ny-1):
            c = j + ny * (i)  # Current index
            n = c + 1  # North neighbor
            s = c - 1  # South neighbor
            e = j + ny * (i+1)  # East neighbor
            w = j + ny * (i - 1)  # West neighbor

            # Filling the M matrix
            M[c, c] = -2 * h[j, i]**3 * (1/(dx**2) + 1/(dy**2))
            M[c, n] = 3 * h.flatten('F')[c]**2 * (h.flatten('F')[n] - h.flatten('F')[s]) / (2 * dy**2) + h.flatten('F')[c]**3 / (dy**2)
            M[c, s] = -3 * h.flatten('F')[c]**2 * (h.flatten('F')[n] - h.flatten('F')[s]) / (2 * dy**2) + h.flatten('F')[c]**3 / (dy**2)
            M[c, e] = 3 * h.flatten('F')[c]**2 * (h.flatten('F')[e] - h.flatten('F')[w]) / (2 * dx**2) + h.flatten('F')[c]**3 / (dx**2)
            M[c, w] = -3 * h.flatten('F')[c]**2 * (h.flatten('F')[e] - h.flatten('F')[w]) / (2 * dx**2) + h.flatten('F')[c]**3 / (dx**2)

            # Filling the rhs vector
            rhs[c] = 6 * u_b * eta_0 * (h.flatten('F')[e] - h.flatten('F')[w]) / (2 * dx)

    p = sps.linalg.spsolve(M, rhs)

    p = p.reshape(ny, nx, order = 'F')


    return X, Y, p

if __name__ == '__main__':
    x, y, p1 = finit_2d()
    x, y, p2 = finit_2d(type = 1)


    X, Y = np.meshgrid(x, y)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, p1)
    
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111, projection='3d')
    ax1.plot_surface(X, Y, p2)
    plt.show()

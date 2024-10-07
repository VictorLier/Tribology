
# In Hamrock chapter 9.2 a numerical method using Fourier series expansion is outlined for
# the analysis of parallel-step thrust bearings. This has been implemented and is available
# at DTU Learn.

# We consider in the following a bearing with length l = 30 mm running in an oil with
# viscosity η = 0.05 Pa s and with an outlet oil film thickness h0 = 20 µm and step height
# sh = 30 µm located at ns = 0.7. The velocity of the slider is ub = 10 m s−1
# .At first the bearing is considered to have width b = l i.e. λ = 1.

# • Investigate the influence of changing the number of terms to include. E.g. plot the
# pressure profile for various jmax = {1, 3, 5..., 99} and investigate the convergence for
# the load Wz

# • Plot a contour of the pressure. Consider (7.38) and (7.39) and plot streamlines in
# the contour plot to visualize the flow. The leakage flow has been integrated over
# the boundaries and is available.

# • Inspect the ratio between leakage and inlet flow.
# The Fourier series expansion implementation might be used as a reference for further
# development of your own general 2D solution using finite differences.

# Write the following matlab code in Python:
# % Parallel step bearing 
# %   Solved using Fourier series expansion
# %   Follows Hamrock et.al. chapter 9
# %
# % Send comments to csan@dtu.dk
# % Casper Schousboe Andreasen 2023
# clear all;
# % parameters
# ns=0.7;     
# l=30e-3;
# b=l;
# sh=30e-6;
# h0=20e-6;
# eta=0.05;
# u=10;

# nx=30; % number of nodes in x direction in both inlet and outlet
# ny=30; % total number of nodes in y-direction
# jmax=19; % Number of terms in series

# % Below values are computed
# H0=h0/sh;
# lambda=l/b;

# % Write data to screen
# fprintf('This bearing have the following properties:\n');
# fprintf('H0: %1.3f\n',H0);
# fprintf('ns: %1.3f\n',ns);
# fprintf('lambda: %1.3f\n',lambda);
# % Setup a grid
# [Xi,Yi]=meshgrid(linspace(0,ns,nx),linspace(0,1,ny));
# [Xo,Yo]=meshgrid(linspace(ns,1,nx),linspace(0,1,ny));
# % Initialize vectors
# Pi=zeros(size(Xi));
# Po=zeros(size(Xo));
# for j=1:2:jmax
#     F=24/(j^2*pi^2*lambda*(H0^3*coth(j*pi*lambda*(1-ns))+(1+H0)^3*coth(j*pi*lambda*ns)));
#     Po=Po+F*sin(j*pi*Yo).*sinh(j*pi*lambda*(1-Xo))/sinh(j*pi*lambda*(1-ns));

#     Pi=Pi+ F*sin(j*pi*Yi).*sinh(j*pi*lambda*Xi)/sinh(j*pi*lambda*ns);
# end
# % Plot the pressure distribution
# figure(1);clf;
# surf(Xi*lambda,Yi,Pi)
# hold on
# surf(Xo*lambda,Yo,Po);

# P=[Pi Po(:,2:end)];
# X=[Xi Xo(:,2:end)];
# Y=[Yi Yo(:,2:end)];
# % Compute the pressure gradients
# [dpidx,dpidy]=gradient(Pi*eta*u*l/sh^2,l*ns/(nx-1),b/(ny-1));
# [dpodx,dpody]=gradient(Po*eta*u*l/sh^2,l*(1-ns)/(nx-1),b/(ny-1));
# % Compute the flows (inlet section and outlet section)
# qix=-(sh+h0)^3/(12*eta)*dpidx+u*(sh+h0)/2;
# qiy=-(sh+h0)^3/(12*eta)*dpidy;
# qox=-h0^3/(12*eta)*dpodx+u*h0/2;
# qoy=-h0^3/(12*eta)*dpody;

# figure(2);clf;
# quiver(Xo,Yo,qox,qoy),hold on
# quiver(Xi,Yi,qix,qiy)
# contour(X,Y,P);
# streamslice(X,Y,[qix qox(:,2:end)],[qiy qoy(:,2:end)]); % duplicate points removed 


# %% Compute the flow across the border
# xx=[Xi Xo(:,2:end)]*l;
# yy=[Yi Yo(:,2:end)]*l/lambda;
# qxx=[qix qox(:,2:end)];
# qyy=[qiy qoy(:,2:end)];

# figure(5);clf;
# plot3(xx(:,1),yy(:,1),qxx(:,1));hold on
# plot3(xx(:,end),yy(:,end),qxx(:,end));
# plot3(xx(1,:),yy(1,:),-qyy(1,:));
# plot3(xx(end,:),yy(end,:),qyy(end,:));
# title('Flow over boundaries')

# inflow=trapz(yy(:,1),qxx(:,1));
# outflow=trapz(yy(:,end),qxx(:,end));
# up=trapz(xx(1,:),-qyy(1,:));
# down=trapz(xx(end,:),qyy(end,:));

# fprintf('Inflow %e, outflow %e, up %e, down %e, leakage %e, leakage%% %f\n',inflow,outflow,up,down,up+down,(up+down)/inflow*100);
# % Check the conservation
# mismatch=inflow-outflow-up-down;
# fprintf('Conservation mismatch: %e mismatch%%: %f\n',mismatch,mismatch/inflow*100);

# start

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import griddata

def coth(x):
    return np.cosh(x) / np.sinh(x)

# parameters
ns=0.7
l=30e-3
b=l
sh=30e-6
h0=20e-6
eta=0.05
u=10

nx=30 # number of nodes in x direction in both inlet and outlet
ny=30 # total number of nodes in y-direction
jmax=19 # Number of terms in series

# Below values are computed
H0=h0/sh
lambda_=l/b

# Write data to screen
print('This bearing have the following properties:')
print(f'H0: {H0:.3f}')
print(f'ns: {ns:.3f}')

# Setup a grid
Xi,Yi=np.meshgrid(np.linspace(0,ns,nx),np.linspace(0,1,ny))
Xo,Yo=np.meshgrid(np.linspace(ns,1,nx),np.linspace(0,1,ny))

# Initialize vectors
Pi=np.zeros(Xi.shape)
Po=np.zeros(Xo.shape)

for j in range(1, jmax+1, 2):
    F=24/(j**2*np.pi**2*lambda_*(H0**3*coth(j*np.pi*lambda_*(1-ns))+(1+H0)**3*coth(j*np.pi*lambda_*ns)))
    Po=Po+F*np.sin(j*np.pi*Yo)*np.sinh(j*np.pi*lambda_*(1-Xo))/np.sinh(j*np.pi*lambda_*(1-ns))
    Pi=Pi+ F*np.sin(j*np.pi*Yi)*np.sinh(j*np.pi*lambda_*Xi)/np.sinh(j*np.pi*lambda_*ns)

# Plot the pressure distribution
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(Xi*lambda_, Yi, Pi, cmap='viridis')
ax.plot_surface(Xo*lambda_, Yo, Po, cmap='viridis')
plt.show()

P=np.concatenate((Pi, Po[:,1:]), axis=1)
X=np.concatenate((Xi, Xo[:,1:]), axis=1)
Y=np.concatenate((Yi, Yo[:,1:]), axis=1)

# Compute the pressure gradients
dpidx,dpidy=np.gradient(Pi*eta*u*l/sh**2, l*ns/(nx-1), b/(ny-1))
dpodx,dpody=np.gradient(Po*eta*u*l/sh**2, l*(1-ns)/(nx-1), b/(ny-1))

# Compute the flows (inlet section and outlet section)
qix=-(sh+h0)**3/(12*eta)*dpidx+u*(sh+h0)/2
qiy=-(sh+h0)**3/(12*eta)*dpidy
qox=-h0**3/(12*eta)*dpodx+u*h0/2
qoy=-h0**3/(12*eta)*dpody

fig, ax = plt.subplots()
ax.quiver(Xo,Yo,qox,qoy)
ax.quiver(Xi,Yi,qix,qiy)
ax.contour(X,Y,P)
# Ensure that X and Y are equally spaced for streamplot
X_stream, Y_stream = np.meshgrid(np.linspace(0, 1, X.shape[1]), np.linspace(0, 1, X.shape[0]))

# Interpolate qix, qiy, qox, qoy to match the new grid

points = np.array([X.flatten(), Y.flatten()]).T
qix_interp = griddata(points, qix.flatten(), (X_stream, Y_stream), method='cubic')
qiy_interp = griddata(points, qiy.flatten(), (X_stream, Y_stream), method='cubic')
qox_interp = griddata(points, qox.flatten(), (X_stream, Y_stream), method='cubic')
qoy_interp = griddata(points, qoy.flatten(), (X_stream, Y_stream), method='cubic')

# Combine interpolated flows
qx_interp = np.concatenate((qix_interp, qox_interp[:, 1:]), axis=1)
qy_interp = np.concatenate((qiy_interp, qoy_interp[:, 1:]), axis=1)

ax.streamplot(X_stream, Y_stream, qx_interp, qy_interp)

# Compute the flow across the border
xx=np.concatenate((Xi, Xo[:,1:]))*l
yy=np.concatenate((Yi, Yo[:,1:]))*l/lambda_
qxx=np.concatenate((qix, qox[:,1:]))
qyy=np.concatenate((qiy, qoy[:,1:]))
fig, ax = plt.subplots()
ax.plot(xx[:,0],yy[:,0],qxx[:,0])
ax.plot(xx[:,-1],yy[:,-1],qxx[:,-1])
ax.plot(xx[0,:],yy[0,:],-qyy[0,:])
ax.plot(xx[-1,:],yy[-1,:],qyy[-1,:])
plt.show()

inflow=np.trapz(yy[:,0],qxx[:,0])
outflow=np.trapz(yy[:,-1],qxx[:,-1])
up=np.trapz(xx[0,:],-qyy[0,:])
down=np.trapz(xx[-1,:],qyy[-1,:])

print(f'Inflow {inflow:e}, outflow {outflow:e}, up {up:e}, down {down:e}, leakage {up+down:e}, leakage% {(up+down)/inflow*100:.2f}')

# Check the conservation
mismatch=inflow-outflow-up-down
print(f'Conservation mismatch: {mismatch:e} mismatch% {mismatch/inflow*100:.2f}')
# end
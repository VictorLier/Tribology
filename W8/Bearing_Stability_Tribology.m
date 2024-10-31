%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRIBOLOGY OF MACHINE ELEMENTS  (41524)                %
% MEK - DEPARTMENT OFCIVIL & MECHANICAL ENGINEERING     %
% DTU - TECHNICAL UNIVERSITY OF DENMARK                 %
%                                                       %
%                 Copenhagen, October 26th, 2023        %
%                                                       %
%                           Ilmar Ferreira Santos       %
%                                                       %                   
%                                                       %
% JOURNAL BEARING DYNAMICS  &  STABILITY ANALYSIS       % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Definition of Nondimensional Bearing Parameters
% S     = eta*N*L*D/W*(R/C)^2;   Sommerfeld Number
% E     = e/C;                  Excentricity
% Phi                            Atittude Angle
% Kij   = (C/W)*kij (i,j=x,y)   dimensionless stiffness coefficients
% Bij   = (C*w/W)*kij (i,j=x,y) dimensionless damping coefficients
%
% Definition of Journal Bearing Parameters
% D                              bearing inner diameter
% L                              bearing width
% R = D/2;                       rotor radius
% C                              bearing clearance
% W                              external load [N] 
% eta                            oil viscosity  [N.s/m2]
% N                              rotor angular velocity [1/s]
% w     = 2*pi*N                 rotor angular velocity [rad/s]

clear all;
close all;

% Bearing Properties 
% Table 1a : Two-axial-groove bearing, L/D = 0.5
%
%        S     E     Phi   Q     P     T   Kxx   Kxy   Kyx   Kyy  Bxx   Bxy  Byx  Byy
%
Table=[6.430 0.071 81.89 0.121 0.860  5.7  1.55 14.41 -6.60 1.88 28.75 1.89 1.89 13.31
       3.937 0.114 77.32 0.192 0.846  5.9  1.57  9.27 -4.20 1.89 18.44 1.93 1.93  8.58
       2.634 0.165 72.36 0.271 0.833  6.2  1.61  6.74 -3.01 1.91 13.36 2.00 2.00  6.28
       2.030 0.207 68.75 0.332 0.835  6.6  1.65  5.67 -2.50 1.93 11.18 2.07 2.07  5.33
       1.656 0.244 65.85 0.383 0.835  7.0  1.69  5.06 -2.20 1.95  9.93 2.15 2.15  4.80
       0.917 0.372 57.45 0.540 0.850  8.5  2.12  4.01 -1.30 1.85  7.70 2.06 2.06  3.23     
       0.580 0.477 51.01 0.651 0.900 10.5  2.67  3.70 -0.78 1.75  6.96 1.94 1.94  2.40 
       0.378 0.570 45.43 0.737 0.977 13.4  3.33  3.64 -0.43 1.68  6.76 1.87 1.87  1.89
       0.244 0.655 40.25 0.804 1.096 17.9  4.21  3.74 -0.13 1.64  6.87 1.82 1.82  1.54
       0.194 0.695 37.72 0.833 1.156 21.3  4.78  3.84  0.01 1.62  7.03 1.80 1.80  1.40
       0.151 0.734 35.20 0.858 1.240 25.8  5.48  3.98  0.15 1.61  7.26 1.79 1.79  1.27
       0.133 0.753 33.93 0.870 1.289 28.7  5.89  4.07  0.22 1.60  7.41 1.79 1.79  1.20
       0.126 0.761 33.42 0.875 1.310 30.0  6.07  4.11  0.25 1.60  7.48 1.79 1.79  1.18
       0.116 0.772 32.65 0.881 1.343 32.2  6.36  4.17  0.30 1.60  7.59 1.79 1.79  1.15
       0.086 0.809 30.04 0.902 1.473 41.4  7.51  4.42  0.47 1.59  8.03 1.79 1.79  1.03
       0.042 0.879 24.41 0.936 1.881 80.9 11.45  5.23  0.92 1.60  9.48 1.80 1.80  0.82 ];
 
% Journal Bearing -- Static Properties   
figure(1)
subplot(2,2,1), plot(Table(:,1),Table(:,2),'*-b','LineWidth',1.5)
title('Excentricity','FontSize',14)
xlabel('S=\eta*N*L*D/W*(R/C)^2','FontSize',14)
ylabel('\epsilon=e/C','FontSize',14)
grid    
subplot(2,2,2), plot(Table(:,1),Table(:,3),'*-b','LineWidth',1.5)
title('Atittude Angle','FontSize',14)
xlabel('S=\eta*N*L*D/W*(R/C)^2','FontSize',14)
ylabel('\phi [^o]','FontSize',14)
grid
subplot(2,2,3.5),
polar(3*pi/2+Table(:,3)*pi/180,Table(:,2))
title('Journal Center Locus','FontSize',14)
grid

% Journal Bearing -- Dynamic Properties (Stiffness)
figure(2)
subplot(2,2,1), plot(Table(:,1),Table(:,7),'*-b','LineWidth',1.5)
xlabel('S=\eta*N*L*D/W*(R/C)^2','FontSize',14)
ylabel('Kxx','FontSize',14)
grid    
subplot(2,2,2), plot(Table(:,1),Table(:,8),'*-b','LineWidth',1.5)
xlabel('S=\eta*N*L*D/W*(R/C)^2','FontSize',14)
ylabel('Kxy')
grid
subplot(2,2,3), plot(Table(:,1),Table(:,9),'*-b','LineWidth',1.5)
xlabel('S=\eta*N*L*D/W*(R/C)^2','FontSize',14)
ylabel('Kyx','FontSize',14)
grid
subplot(2,2,4), plot(Table(:,1),Table(:,10),'*-b','LineWidth',1.5)
xlabel('S=\eta*N*L*D/W*(R/C)^2','FontSize',14)
ylabel('Kyy','FontSize',14)
grid

% Journal Bearing -- Dynamic Properties (Damping)
figure(3)
subplot(2,2,1), plot(Table(:,1),Table(:,11),'*-b','LineWidth',1.5)
xlabel('S=\eta*N*L*D/W*(R/C)^2','FontSize',14)
ylabel('Bxx','FontSize',14)
grid
subplot(2,2,2), plot(Table(:,1),Table(:,12),'*-b','LineWidth',1.5)
xlabel('S=\eta*N*L*D/W*(R/C)^2','FontSize',14)
ylabel('Bxy','FontSize',14)
grid
subplot(2,2,3), plot(Table(:,1),Table(:,13),'*-b','LineWidth',1.5)
xlabel('S=\eta*N*L*D/W*(R/C)^2','FontSize',14)
ylabel('Byx','FontSize',14)
grid
subplot(2,2,4), plot(Table(:,1),Table(:,14),'*-b','LineWidth',1.5)
xlabel('S=\eta*N*L*D/W*(R/C)^2','FontSize',14)
ylabel('Byy','FontSize',14)
grid


% bearing geometry
 D     = 0.100 ;    % [m]
 L     = 0.050 ;    % [m]
 R     = D/2;       % [m]
 r     = 0.0499 ;   % [m]
 C    = 100e-6 ;   % [m]
 C    = R-r;       % [m]
% oil properties  
 eta    = 0.0271;    % [N.s/m2]
% acceleration of gravity 
 g     = 9.8;       % [m/s^2]
% rotor mass and weight 
 mass  = 500;       % [kg]
 W     = mass*g/2;  % [N] 
%  N = rotational velocity [1/s] 
%  N = Table(1,1)*(W)/(eta*L*D)/(R/C)^2
%  N = 3.0376;   line = 16;
%   N = 6.2199;   line  = 15;
%  N = 8.3897;   line = 14;
%  N = 9.1129;   line = 13;
%  N = 9.6192;   line = 12;
%  N = 10.9210;  line = 11;
%  N = 14.0310;  line = 10;
%  N = 17.6472;  line = 9;
%   N = 27.3387;  line = 8 ;
%   N = 41.9483;  line = 7;
   N = 66.3218;  line = 6; % STABLE
%  N = 119.7697; line = 5; % UNSTABLE
%  N = 146.8192; line = 4;
%  N = 190.5033; line = 3;
%  N = 284.7424; line = 2;
%  N = 465.0480; line = 1;

 w     = 2*pi*N;   % [rad/s]
 S     = eta*N*L*D/W*(R/C)^2 % Sommerfeld Number
    
% Coefficients of the Stiffness Matrix  
kxx =  Table(line,7)*W/C  ; % [N/m]  
kxy =  Table(line,8)*W/C  ; % [N/m]
kyx =  Table(line,9)*W/C  ; % [N/m]
kyy =  Table(line,10)*W/C ; % [N/m]

% Coefficients of the Damping Matrix 
bxx=  Table(line,11)*W/(w*C)  ; % [N/(m/s)]
bxy=  Table(line,12)*W/(w*C)  ; % [N/(m/s)]
byx=  Table(line,13)*W/(w*C)  ; % [N/(m/s)]
byy=  Table(line,14)*W/(w*C)  ; % [N/(m/s)]

x_st   = Table(line,2)*cos(Table(line,3)*pi/180)*C; % equilibrium position [m]
y_st   = Table(line,2)*sin(Table(line,3)*pi/180)*C; % equilibrium position [m]


%Mass Matrix
M= [mass 0; 0 mass];

%Damping Matrix
B=[bxx bxy; byx byy];

%Stiffness Matrix
K= [kxx kxy; kyx kyy];

%State Matrix
A1= [        M        zeros(size(M))    ; 
    zeros(size(M))   M  ] ;

A2= [       B         K         ; 
           -M        zeros(size(M))];

%Dynamical Properties of the Mass-Spring System
[u,s]=eig(-A2,A1);          % eigenvectors [ ] & eigenvalues [rad/s]
diag(s)                     % eigenvalues [rad/s] 
xi=diag(-real(s)/abs(s))    % damping factors [ ]
wd_hz=diag(imag(s))/(2*pi)  % damped natural frequencies [Hz]
 

%_____________________________________________________
%Inicial Condition (Perturbation)
 x_ini    =  0.000;       % initial displacement [m]
 y_ini    =  0.000;       % initial displacement [m]
 vx_ini   =  0.03;        % initial velocity [m/s]
 vy_ini   =  0.000;       % initial velocity [m/s]
 time_max =  0.26;        % integration time [s]
 
%_____________________________________________________

%_____________________________________________________
%EXACT SOLUTION 
n=2000;                % number of points for plotting
j=sqrt(-1);            % complex number 

% initial conditions of velocity and displacement
z_ini = [vx_ini vy_ini x_ini y_ini]';

% poles s or eigenvalues 
s1=s(1,1);
s2=s(2,2);
s3=s(3,3);
s4=s(4,4);

% eigenvectors
u1=u(1:4,1);
u2=u(1:4,2);
u3=u(1:4,3);
u4=u(1:4,4);

% calculation of the constants based on the initial conditions of movement

C_ini=inv(u)*(z_ini);

c1=C_ini(1);
c2=C_ini(2);
c3=C_ini(3);
c4=C_ini(4);

% Homogenous solutions based on modal superposition

for i=1:n,
    t(i)=(i-1)/n*time_max;
    z_exact = c1*u1*exp(s1*t(i)) + ...
              c2*u2*exp(s2*t(i)) + ...
              c3*u3*exp(s3*t(i)) + ...
              c4*u4*exp(s4*t(i));
           
 % cartesian coordinates          
    x_exact(i)     = x_st + real(z_exact(3));
    y_exact(i)     = y_st + real(z_exact(4));
    
 % polar coordinates   
    polar_rho(i)   = sqrt(x_exact(i)*x_exact(i)+y_exact(i)*y_exact(i));
    polar_angle(i) = acos(x_exact(i)/polar_rho(i));
end

% Journal Bearing Stability - Time Domain Analysis
figure(4)
subplot(2,2,1), plot(t,x_exact,'b','LineWidth',1.5)
title('Vertical Direction ','FontSize',14)
xlabel('time [s]','FontSize',14)
ylabel('x(t) [m]','FontSize',14)
grid
subplot(2,2,2), plot(t,y_exact,'b','LineWidth',1.5)
title('Horizontal Direction ','FontSize',14)
xlabel('time [s]','FontSize',14)
ylabel('y(t) [m]','FontSize',14)
grid
subplot(2,2,3), plot(y_exact,x_exact,'b','LineWidth',1.5)
title('Orbit of Journal Center','FontSize',14)
xlabel('y(t) [m]','FontSize',14)
ylabel('x(t) [m]','FontSize',14)
grid

subplot(2,2,4),
polar(polar_angle*0,polar_rho/C*0)
hold on
polar(3*pi/2+polar_angle,polar_rho/C)
title('Journal Center Orbit','FontSize',14)
grid

figure(5)
polar(polar_angle*0,polar_rho/C*0)
hold on
polar(3*pi/2+polar_angle,polar_rho/C)
title('Journal Center Orbit','FontSize',14)
grid








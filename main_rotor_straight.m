%%Main file for the structural calculations

clc; clear; close all

%% Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 100;                                %Number of lengthwise elements
 
%geometric parameters
%rotor geometry
h = 60;                                 %rotor height in m
D = 30;                                 %rotor diameter in m
phase = 0;                              %Phase shift angle in degrees


%airfoil properties
airfoil.filename = 'NACAM23.txt';      	%Airfoil file
airfoil.t = 3*10^-3;                 	%Thickness in m
airfoil.C = 1.3;                      	%Chord length in m

theta_p = 0;                            %pitch angle 

%Airodynamic properties
lambda = 3;                           	%tip speed ratio
Vw = 12;                             	%Wind speed
qy = 1e3;
qx = 1e3;

%material properties
E = 1e5;                                %youngsmodulus
rho_mat = 2702;                        	%material density in kg/m3

%general properties
g = 9.81;                              	%gravitational acceleration





%% Airfoil moments of area, include in function
[Ixx, Iyy, Ixy, area] = airfoilInertia(airfoil,theta_p);
[Ixx0, Iyy0, ~, ~] = airfoilInertia(airfoil,0);
%% CHECK IF IXY > 0 


%% Centrifugal force
Mtot = h*area*rho_mat;                  %total weight
Wtot = Mtot*g;                          %total weight

omega  = Vw*lambda/(0.5*D);                               %Calculates rotational speed in rad/s
Fc     = Mtot*omega^2*(0.5*D);

% COORDINATE SYSTEM??????????????????

% for i = 1:length(r)
%     Fc(i) = M(i)*r(i)*omega^2;                            %Calculates magnitude centrifugal force in N
%     xFc(i) = (r(i)+Fc(i)/(10000./N))*cos(theta(i));       %Make array of position points of the centrifugal force
%     zFc(i) = (r(i)+Fc(i)/(10000./N))*sin(theta(i));
% end



%% call matrix method
%The first column indicates the distance from the start of the rotor
%(bottom up). The second column is the degre of freedom. 
% DoF 1 : positive y direction (Force in y)
% Dof 2 : positve angle about x (Moment about x)
% DoF 3 : positive x direction (Force in x)
% DoF 4 : positive angle about y (moment about y)
% The third column indicates the set value of the boundary condtion (0
% means no displacement) in bcs. In matrix R, the third column represents
% the reaction forces.

bcs = [[0.0, 1, 0];
       [0.0, 2, 0];
       [0.0, 3, 0];
       [0.0, 4, 0];
       [h, 1, 0];
       [h, 2, 0];
       [h, 3, 0];
       [h, 4, 0]];

[R] = MatrixMethod(h,theta_p,qy,qx,Ixx0,Iyy0,E,bcs)


%% Reaction force calculations
% Fy = 0.5*sum(W);                                                            %reaction force in axis direction
% Frup = 0;
% for i = 1:length(r)                                                         %
%     Frup = Frup + (Fc(i)*y(i)+W(i)*r(i))/h;
% end
% Frupx = Frup*cos(mean(theta));
% Frupz = Frup*sin(mean(theta));
% Frlow = sum(Fc)-Frup;
% Frlowx = Frlow*cos(mean(theta));
% Frlowz = Frlow*sin(mean(theta));
% 
% plot3([x(length(y)),x(length(y))-Frupx/(10000/N)],[z(length(y)),z(length(y))-Frupz/(10000/N)],[y(length(y)),y(length(y))],'black:')
% plot3([x(1),x(1)-Frlowx/(10000/N)],[z(1),z(1)-Frlowz/(10000/N)],[y(1),y(1)],'black:')
% 
% %____Internal Forces_______________________________________________________
% 
% j = -Fy;
% jx = -Frlowx;
% jz = -Frlowz; 
% for i = 1:length(y)
%     j = j + W(i);
%     jx = jx + Fc(i)*cos(theta(i));
%     jz = jz + Fc(i)*sin(theta(i));
%     intfn(i) = j;
%     intfvx(i) = jx;
%     intfvz(i) = jz;
% end
% 
% %____Stresses______________________________________________________________
% tauxyV = intfvx/Aair;
% tauzyV = intfvz/Aair;
% sigmaN = intfn/Aair;


%%  Plot forces on beam
% %____Plot weight___________________________________________________________
% hold on
% plot3(x, z, yW,'r--')                        	%plots the line of the weights
% for i = 1:length(y)
%    plot3([x(i),x(i)],[z(i),z(i)],[y(i),yW(i)],'r')
% end

%____Plot centrifugal force________________________________________________
% plot3(xFc, zFc, y, 'y--')
% for i = 1:length(y)
%    plot3([x(i),xFc(i)],[z(i),zFc(i)],[y(i),y(i)],'y')
% end

% fig = plt.figure(4)
% ax = fig.gca(projection='3d')
% ax.scatter(tauxyV, tauzyV, sigmaN)
% #ax.legend()
% #plt.grid(True)
% #plt.title("Stress envalope")
% plt.show()
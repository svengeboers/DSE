function [Cl, Cd] = liftdrag(alpha)

% Cl = 2*pi*sin(alpha)*1.1;
% Cd = 0;%0.023;%*ones(1,length(alpha));
<<<<<<< HEAD

filename = '/Users/Emile/Desktop/DSE/Polars/NACA0015_RE265000.txt';
delimiterIn = ' ';
headerlinesIn = 12;

A = importdata(filename,delimiterIn, headerlinesIn);

alphas = deg2rad(A.data(:,1));
Cls = A.data(:,2);
Cds = A.data(:,3);

Cl = interp1(alphas',Cls,alpha);
Cd = interp1(alphas',Cds,alpha);

end
=======
% clc; clear all; close all;
% alpha = [-3:3];

filename = 'NACA0021_Re3.000_M0.02_N9.0_360_V.txt';
delimiterIn = ' ';
headerlinesIn = 14;

A = importdata(filename, delimiterIn, headerlinesIn);

alphas = A.data(:,1);
Cls = A.data(:,2);
Cds = A.data(:,3);

Cl = interp1(alphas',Cls,rad2deg(alpha));
Cd = interp1(alphas',Cds,rad2deg(alpha));
end

% load('D:\damdetavernier\Desktop\PhD\7 - Matlab\7.2 - VAWT Codes\Actuator Cylinder - Fixed Pitch\Actuator Cylinder - Fixed Pitch V2\AIRFOILS\NACA0018_RE5.1e+06_pol.mat')
% Cl = interp1(deg2rad(POLAR(1).alpha),POLAR(1).CL,alpha);
% Cd = interp1(deg2rad(POLAR(1).alpha),POLAR(1).CD,alpha);
>>>>>>> 08c3dc90e4ecbb0435634a0b2fceb83c694e3fd1

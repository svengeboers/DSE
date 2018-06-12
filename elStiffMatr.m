function Kel = elStiffMatr(L,E,Izz,Iyy,theta)
	
% 	Description: generates the element stiffness matrix for
% 				 a beam element, 8 DoF per element (4 per node)
% 	
% 	Inputs:
% 		L       =  length of beam element [m]
% 		E       =  Young's Modulus [Pa]
% 		Izz,Iyy =  Moments of inertia [m^4]
% 		theta   =  Rotation to global coordinate system [deg]
% 	
% 	Outputs:
% 		Kel   =  Element stiffness matrix in global coordinates
% 
% %testing
% clear all; close all;
% theta = 0;
% Izz = 2;
% Iyy = 1;
% L = 0.1;
% E = 100000;
% %end testing	

th = deg2rad(theta);

Kel = zeros(8,8);

% Load element stifness matrix in local coordinates
P1 = [[12,6*L];[6*L,4*L*L]];
P2 = [[-12,6*L];[-6*L,2*L*L]]; P3 = [[12,-6*L];[-6*L,4*L*L]];

Kel(1:2,1:2) = P1*Izz;
Kel(3:4,3:4) = P1*Iyy;
Kel(5:6,1:2) = (P2*Izz)';
Kel(7:8,3:4) = (P2*Iyy)';
Kel(1:2,5:6) = P2*Izz;
Kel(3:4,7:8) = P2*Iyy;
Kel(5:6,5:6) = P3*Izz;
Kel(7:8,7:8) = P3*Iyy;

Kel = Kel*(E/(L^3));


%%%%% IMPORTANT NOTE %%%%%%
% The curent method of loading the element stiffness matrix, and then rotatating
% the matrix to account for a roted beam only works for beams that have Izy = 0
% when not rotated. To do the analysis for beams with Izy nonzero, the equations to
% load the element matrix should be adjusted to include for Izy.
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define Rotation matrix T
T = zeros(8,8);
P1 = [[cos(th), 0, -sin(th) , 0];
    [0,cos(th),0,sin(th)];
    [sin(th),0,cos(th),0];
    [0,-sin(th),0,cos(th)]];
T(1:4,1:4) = P1;
T(5:end,5:end) = P1;

Kel =  T*Kel*(T');


end
% In this program the structural loads will be calucalated. The paper
% used for this is "structural analysis of a curved beam element defined
% in global coordinates". The idea of the matrix method used here
% is to integrate the distributed loads on the structure which will result
% in the shear forces on the structrure, then the moments, twist angles
% and displacements can be calculated subsequentially. The input of the 
% matrix will thus be the distributed loads of gravity, lift and drag
% in x, y and z direction.

load('/home/sgeboers/Dropbox/DSE/DSE/QX.mat')
load('/home/sgeboers/Dropbox/DSE/DSE/QY.mat')
load('/home/sgeboers/Dropbox/DSE/DSE/S.csv')

% First attempt at trying to find an expression for the forces on the
% wind turbine blade
% We try to interpolate the load equation for each 
% angle among the curve of angles
Z = 15:10:55;

% Calculate the s values for certain the Z values at which the forces are
% calculated
Sspline = spline(S(1:end,1),S(1:end,2),Z);

% Find a spline that relates the QX values to the spline values and graph
% them out
QX(1:42:end);
pp = spline(Sspline,QX(1:42:end),1:1:5);
size(QX(1:42:end));

% What follows are placeholder expressions for the distributed forces
% in x, y and z direction

syms s sb s1 s2 st
qx = 3*s^2+2;
qy = 3*s^2+2;
qz = 3*s^2+2;

% Definite integrals for previously mentioned expressions
int(qx);
int(qy);
int(qz);

S = 0:0.1:50;

%% Defining constants for in THE MATRIX
a = 5;
b = 3;
c = 3.5;
E = 1000;
A = 1700;
G = 3000;
It = 2000;
In = 1500;
Ib = 2100;


% Transformation matrix
vtx = -a/c*sin(s/c);
vty = a/c*cos(s/c);
vtz = b/c;
vnx = -c/c*cos(s/c);
vny = -c/c*cos(s/c);
vnz = 0;
vbx = b/c*sin(s/c);
vby = -b/c*cos(s/c);
vbz = a/c;

% Other coefficients
exx = a^2*sin(s/c)^2/(c^2*E*A);
eyy = a^2*cos(s/c)^2/(c^2*E*A);
ezz = b^2/(c^2*E*A);
exy = -a^2*sin(2*s/c)/(2*c^2*E*A);
exz = -a*b*sin(s/c)/(c^2*E*A);
eyz = a*b*cos(s/c)/(c^2*E*A);

yxx = a^2*sin(s/c)^2/(c^2*G*It)+cos(s/c)^2/(E*In)+b^2*sin(s/c)^2/(c^2*E*Ib);
yyy = a^2*cos(s/c)^2/(c^2*G*It)+sin(s/c)^2/(E*In)+b^2*cos(s/c)^2/(c^2*E*Ib);
yzz = b^2/(c^2*G*It)+a^2/(c^2*E*Ib);
yxy = -a^2*sin(2*s/c)/(2*c^2*G*It)+sin(2*s/c)^2/(2*E*In)-b^2*sin(2*s/c)^2/(c^2*E*Ib);
yxz = -a*b*sin(s/c)/(c^2*G*It)+a*b*sin(s/c)/(c^2*E*Ib);
yyz = a*b*cos(s/c)/(c^2*G*It)-a*b*cos(s/c)/(c^2*E*Ib);

% Load step function
syms Fxb Fx1 Fx2 Fxt Fyb Fy1 Fy2 Fyt Fzb Fz1 Fz2 Fzt
syms c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 c12

Vx = Fxb + Fx1*heaviside(s-s1) + Fx2*heaviside(s-s2) - qx*s + c1;
Vy = Fyb + Fy1*heaviside(s-s1) + Fy2*heaviside(s-s2) - qy*s + c2;
Vz = Fzb + Fz1*heaviside(s-s1) + Fz2*heaviside(s-s2) - qz*s + c3;

% Moment equations
Mx = int(vtz*Vy-vty*Vz) + c4;
My = int(vtx*Vz-vtz*Vx) + c5;
Mz = int(vty*Vx-vtx*Vy) + c6;

% Angular displacements
thetax = int(yxx*Mx+yyx*My+yzx*Mz) + c7;
thetay = int(yxy*Mx+yyy*My+yzy*Mz) + c8;
thetaz = int(yxz*Mx+yyz*My+yzz*Mz) + c9;

% Displacements
deltax = vtz*thetay-vty*thetaz + c10;
deltay = vtx*thetaz-vtz*thetax + c11;
deltaz = vty*thetax-vtx*thetay + c12;






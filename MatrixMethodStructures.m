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
% in x, y and z direction. For the aerodynamic forces, for now a constant
% distribution is used instead of a variable distributed load

syms s 
sb = 0; 
s1 = 0; %Placeholder values for s1 and s2
s2 = 0; 
qx = 0*s+300; %Placeholder values
qy = 0*s+50;
qz = 0*s+9.80665;

% Definite integrals for previously mentioned expressions
int(qx);
int(qy);
int(qz);

%% Defining constants for THE MATRIX
a = 15; %this is the radius of the helix
H = 45; %Height of the rotor
b = H*3/pi; %Constant that defines the helix
c = sqrt(a^2+b^2); %Constant that defines the helix
S = sb:0.1:st; %Array which will be used to plot things along the length of the beam later on
st = c*pi/3; %Length of the curve
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

% Reaction forces
Fxb = 200;
Fx1 = 0;
Fx2 = 0;
Fxt = 500;
Fyb = 300;
Fy1 = 0;
Fy2 = 0;
Fyt = 600;
Fzb = 400;
Fz1 = 0;
Fz2 = 0;
Fzt = 700;

syms c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 c12

Vx = Fxb + 0*Fx1*heaviside(s-s1) + 0*Fx2*heaviside(s-s2) + 0*Fxt*heaviside(s-st) - qx*s + c1;
Vy = Fyb + 0*Fy1*heaviside(s-s1) + 0*Fy2*heaviside(s-s2) + 0*Fyt*heaviside(s-st) - qy*s + c2;
Vz = Fzb + 0*Fz1*heaviside(s-s1) + 0*Fz2*heaviside(s-s2) + 0*Fzt*heaviside(s-st) - qz*s + c3;

% Moment equations
Mx = int(vtz*Vy-vty*Vz) + c4;
My = int(vtx*Vz-vtz*Vx) + c5;
Mz = int(vty*Vx-vtx*Vy) + c6;

% Angular displacements
thetax = int(yxx*Mx+yxy*My+yxz*Mz) + c7;
thetay = int(yxy*Mx+yyy*My+yyz*Mz) + c8;
thetaz = int(yxz*Mx+yyz*My+yzz*Mz) + c9;

% Displacements
deltax = vtz*thetay-vty*thetaz + c10;
deltay = vtx*thetaz-vtz*thetax + c11;
deltaz = vty*thetax-vtx*thetay + c12;

% Calculation of constant c1, c2 and c3
% Shear force at the end should equal the reaction force

eqn1 = subs(Vx,s,st)==-Fxt;
eqn2 = subs(Vy,s,st)==-Fyt;
eqn3 = subs(Vz,s,st)==-Fzt;

% Calculation of constant c4, c5 and c6
% Moments yet to be obtained from Louis K

eqn4 = subs(Mx,s,st)==0;
eqn5 = subs(My,s,st)==0;
eqn6 = subs(Mz,s,st)==0;

% Calculation of constant c7, c8 and c9
% angular displacements at the end should be 0

eqn7 = subs(thetax,s,st)==0;
eqn8 = subs(thetay,s,st)==0;
eqn9 = subs(thetaz,s,st)==0;

% Calculation of constant c10, c11 and c12
% Displacement at the end should be 0 (clamped)
eqn10 = subs(deltax,s,st)==0;
eqn11 = subs(deltay,s,st)==0;
eqn12 = subs(deltaz,s,st)==0;


[A,B] = equationsToMatrix([eqn1 eqn2 eqn3 eqn4 eqn5 eqn6 eqn7 eqn8 eqn9 eqn10 eqn11 eqn12],[c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 c12]);

X = linsolve(A,B);

% Now that the solution for the constants is known, they can be used as an
% input to rewrite the equations. What will follow is just a substitution
% of the constants into the equations.

Vx = Fxb + 0*Fx1*heaviside(s-s1) + 0*Fx2*heaviside(s-s2) + 0*Fxt*heaviside(s-st) - qx*s + X(1);
Vy = Fyb + 0*Fy1*heaviside(s-s1) + 0*Fy2*heaviside(s-s2) + 0*Fyt*heaviside(s-st) - qy*s + X(2);
Vz = Fzb + 0*Fz1*heaviside(s-s1) + 0*Fz2*heaviside(s-s2) + 0*Fzt*heaviside(s-st) - qz*s + X(3);

% Moment equations
Mx = int(vtz*Vy-vty*Vz) + X(4);
My = int(vtx*Vz-vtz*Vx) + X(5);
Mz = int(vty*Vx-vtx*Vy) + X(6);

% Angular displacements
thetax = int(yxx*Mx+yxy*My+yxz*Mz) + X(7);
thetay = int(yxy*Mx+yyy*My+yyz*Mz) + X(8);
thetaz = int(yxz*Mx+yyz*My+yzz*Mz) + X(9);

% Displacements
deltax = vtz*thetay-vty*thetaz + X(10);
deltay = vtx*thetaz-vtz*thetax + X(11);
deltaz = vty*thetax-vtx*thetay + X(12);


Vx = double(subs(Vx,s,S));
Vy = double(subs(Vy,s,S));
Vz = double(subs(Vz,s,S));

Mx = double(subs(Mx,s,S));
My = double(subs(My,s,S));
Mz = double(subs(Mz,s,S));

thetax = double(subs(thetax,s,S));
thetay = double(subs(thetay,s,S));
thetaz = double(subs(thetaz,s,S));

deltax = double(subs(deltax,s,S));
deltay = double(subs(deltay,s,S));
deltaz = double(subs(deltaz,s,S));


%plot(S,Vx)
%plot(S,Mx)
%plot(S,thetax)
%plot(S,deltax)
























clc; clear all; close all;
tic

%% Input parameters

% Blade
B = 3;                      % Number of blades [-]        
c = 1.6;            % Chord length [m]
twist = 0;                  % Twist angle [rad] [keep 0]

% Flow conditions
Vinf = 6;                   % Incoming wind speed, from west to east [m/s]
rho = 1.225;                % Air density [kg/m^3]
omega = lambda*Vinf/R;      % Rotational speed, counter-clockwise(+) [rad/s]
lambda = 3;                 % Tip speed ratio [-]
theta_p = deg2rad(0);       % Pitch angle [rad]


% Discretisation
H = 60;         %Height [m]
dZ = 10;       %Spacing [m]
N = H/dZ       %Number of sections

%% Discretisation

Ptot = 0;
Ttot = 0;

for i = 1:1:N
    h = i*dZ;                   % height
    R = h*0.5;                  % Radius
    sigma = B*c/(2*R);          % Solidity
    

    [Cp,CT] = actuatorcylinder2(sigma,B,R,Vinf,lambda);
    Ptot = Ptot + Cp*0.5*rho^2*Vinf^3*2*R;
    Ttot = Ttot + CT*0.5*rho^2*Vinf^3*2*R;
end

Ttot
Ptot

toc
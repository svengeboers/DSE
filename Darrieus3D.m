clc; clear all; close all;
tic

% 

%% Input parameters

% Blade
B = 3;                      % Number of blades [-]        
c = 1.6;                    % Chord length [m]
twist = 0;                  % Twist angle [rad] [keep 0]

% Flow conditions
%Vinf = 6;                   % Incoming wind speed, from west to east [m/s]
rho = 1.225;                % Air density [kg/m^3]
omega = 2;                  % Rotational speed, counter-clockwise(+) [rad/s]
theta_p = deg2rad(0);       % Pitch angle [rad]
gamma = deg2rad(10);         % sweep angle [rad] (change in azamuthal angle per meter)

% Discretisation
H = 60;                             %Total height [m]
dZ = 10;                            %Section height [m]
Htower = 5;                         % Tower height [m]
N = H/dZ;                           %Number of sections
f = @(x) -(1/60)*(x-30)^2+15;       %Radius function, x = height
df = @(x) -(2/60)*x+(30/60);        

%% Discretisation

Ptot = 0;
Ttot = 0;

for i = 0:1:N
    h = i*dZ;                   % height
    R = f(h);                   % Radius
    htot = h+Htower;            % Total height incl tower [m]
    Vinf = Weibull2(htot);      % Incoming wind speed, from west to east [m/s]
    sigma = B*c/(2*R);          % Solidity
    lambda = omega*R/Vinf;      % Tip speed ratio [-]
    delta = atan(1/df(h));      % Angle of curvature [rad]
    Dtheta = h*gamma;           % Azimuthal offset [rad]
    
    if lambda >= 2.07
        [Cp,CT] = actuatorcylinder2(B,R,Vinf,lambda,delta,c,Dtheta);
        Ptot = Ptot + dZ*Cp*0.5*rho^2*Vinf^3*2*R;
        Ttot = Ttot + dZ*CT*0.5*rho^2*Vinf^3*2*R;
    end
    
end

Ttot
Ptot

toc
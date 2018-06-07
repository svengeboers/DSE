clc; clear all; close all;
tic

% This program calculates the thrust, power, torque and forces on the
% blades of the 3D turbine. It discretises the turbine over the hight in N
% sections. At each section the actuator cylinder is applied.

%% Input parameters

% Blade
B = 3;                      % Number of blades [-] Make sure N in actuatorcylinder2 can be divided by B
beta = deg2rad(360/B);      % Angle between the blades [rad]
c = 1.6;                    % Chord length [m]
twist = 0;                  % Twist angle [rad] [keep 0]

% Flow conditions
rho = 1.225;                % Air density [kg/m^3]
omega = 2;                  % Rotational speed, counter-clockwise(+) [rad/s]
theta_p = deg2rad(0);       % Pitch angle [rad]
gamma = deg2rad(10);        % Sweep angle [rad] (change in azamuthal angle of the blade per meter increase in height)

% Discretisation
H = 60;                             %Total height [m]
dZ = 10;                            %Section height [m]
Htower = 5;                         %Tower height [m]
M = H/dZ;                           %Number of sections
f = @(x) -(1/60)*(x-30)^2+15;       %Radius function, x = height, the function of the curvature of the blade
df = @(x) -(2/60)*x+1;        

%% Discretisation

Ptot = 0;
Ttot = 0;
Atot = 0;
Flst = [];
CTlst = [];
Cplst = [];
k = 0;

for i = 0:1:M
    h = i*dZ;                   % height
    R = f(h);                   % Radius
    htot = h+Htower;            % Total height incl tower [m]
    Vinf = Weibull2(htot);      % Incoming wind speed, from west to east [m/s]
    sigma = B*c/(2*R);          % Solidity
    lambda = omega*R/Vinf;      % Tip speed ratio [-]
    delta = -atan(df(h));       % Angle of curvature [rad]
    Dtheta = h*gamma;           % Azimuthal offset [rad]
    
    if lambda >= 2.07
        [Cp,CT,theta,Qn,Qz,Qt,~] = actuatorcylinder2(B,R,Vinf,lambda,delta,c);
        Ptot = Ptot + dZ*Cp*0.5*rho^2*Vinf^3*2*R;       
        Ttot = Ttot + dZ*CT*0.5*rho^2*Vinf^3*2*R;
        
        Qx = Qt.*cos(theta)+Qn.*sin(theta);        % transform the loads to x,y,z coordinates
        Qy = Qt.*sin(theta)-Qn.*cos(theta);
        
        loadlst = [];
        Cplst = [Cplst;[h,Vinf,Cp]];
        CTlst = [CTlst;[h,Vinf,CT]];
        
        for i = 1:1:length(Qx)                      % create a matrix with the height, Radius, azimuthal angle and loads in x,y,z (and Qt for calculation of torque)
            loadlst = [loadlst;[h,R,theta(i),Qt(i),Qx(i),Qy(i),Qz(i)]];
        end
        Flst = [Flst;loadlst];                      % Flst = [[h,R,theta,Qt,Qx,Qy,Qz];[h,R,theta,Qt,Qx,Qy,Qz]...[h,R,theta,Qt,Qx,Qy,Qz];[h,R,theta,Qt,Qx,Qy,Qz]]
        k = k+1;
    end
    
end

%% Calculate the torque

[~,~,~,~,~,~,N] = actuatorcylinder2(3,1,4,3,0,1);

Torquelst = [];

for i = 1:1:N               %Calculate the torque over the rotor for each position
    Torque = 0;
    for j = 1:1:k           %Calculate the torque over the whole rotor for one position
        n = k*N-N+1;
        B1 = Flst(n,4)*Flst(n,2);
        B2 = Flst((n+N/B),4*Flst(n,2);
        B3 = Flst((n+2*N/B),4*Flst(n,2);
        Torque = Torque+B1+B2+B3;
    end
    
end

%% Print thrust, power

Ttot        %total thrust of the rotor
Ptot        %total power of the rotor

toc
function [Cptot,CTtot] = fDarrieus3D(c,omega,B,gamma)

%% Input parameters

% Blade
B = B;                      % Number of blades [-] Make sure N in actuatorcylinder2 can be divided by B
beta = deg2rad(360/B);      % Angle between the blades [rad]
c = c;                      % Chord length [m]
twist = 0;                  % Twist angle [rad] [keep 0]

% Flow conditions
rho = 1.225;                % Air density [kg/m^3]
omega = omega;              % Rotational speed, counter-clockwise(+) [rad/s]
theta_p = deg2rad(0);       % Pitch angle [rad]
gamma = gamma;              % Sweep angle [rad] (change in azamuthal angle of the blade per meter increase in height)

% Discretisation
H = 60;                             %Total height [m]
dZ = 10;                            %Section height [m]
Htower = 5;                         %Tower height [m]
M = H/dZ;                           %Number of sections
f = @(x) -(1/H)*(x-H/2)^2+H/4;      %Radius function, x = height, the function of the curvature of the blade
df = @(x) -(2/H)*x+1;      

%% Discretisation

Ptot = 0;
Ttot = 0;
Cptot = 0;
CTtot = 0;
Atot = 0;
Flst = [];
CTlst = [];
Cplst = [];
k = 0;

for i = 1:1:M
    h = i*dZ-0.5*dZ;            % height
    R = f(h);                   % Radius
    htot = h+Htower;            % Total height incl tower [m]
    Vinf = Weibull2(htot);      % Incoming wind speed, from west to east [m/s]
    sigma = B*c/(2*R);          % Solidity
    lambda = omega*R/Vinf;      % Tip speed ratio [-]
    delta = -atan(df(h));       % Angle of curvature [rad]
    
    [Cp,CT,theta,Qn,Qz,Qt,~] = actuatorcylinder2(B,R,Vinf,lambda,delta,c);
    
    Ptot = Ptot + dZ*Cp*0.5*rho^2*Vinf^3*2*R;       
    Ttot = Ttot + dZ*CT*0.5*rho^2*Vinf^3*2*R;
    Cptot = Cptot + Cp*2*R*dZ;
    CTtot = CTtot + CT*2*R*dZ;
    Atot = Atot + 2*R*dZ;

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

Cptot = Cptot/Atot;   %total Cp of the rotor
CTtot = CTtot/Atot;  %total Ct of the rotor
Ttot;                 %total thrust of the rotor
Ptot;                 %total power of the rotor
Torque_estimated = Ptot/omega;
lambda;

end
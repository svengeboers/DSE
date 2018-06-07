function [Cp,CT,theta,Qn,Qz,Qt,N] = actuatorcylinder2(B,R,Vinf,lambda,delta,c)
%% Input parameters
% Blade
sigma = B*c/(2*R);              % Solidity [-]
B = B;                      % Number of blades [-]
R = R;                      % Radius [m]
c = c;                      % Chord length [m]
twist = 0;                  % Twist angle [rad] [keep 0]

% Flow conditions
Vinf = Vinf;                % Incoming wind speed, from west to east [m/s]
rho = 1.225;                % Air density [kg/m^3]
lambda = lambda;            % Tip speed ratio [-]
theta_p = deg2rad(0);       % Pitch angle [rad]
omega = lambda*Vinf/R;      % Rotational speed, counter-clockwise(+) [rad/s]

% Additional parameters
N = 42;                     % Number of cells in cylinder [-]   [ > 36 ]
Ngrid = 16;                 % Number of cells in grid [-]       [dont change]
f = 1.001;                  % Factor for evaluation points [-]  [dont change]
tol = 1e-7;                 % Tolerance for while loop [-]      [dont change]
relax = 0.1;                % Relaxation paramter [-]           [dont change]

%% Define control points
% Discretise control points
dtheta = 2*pi/N;
theta = (1:N)*dtheta-1/2*dtheta;

xc = -R*sin(theta);
yc = R*cos(theta);

% Discretise evaluation points
xe = f*xc/R;
ye = f*yc/R;


%% Influence matrix
% Prelocation
Rwx = zeros(N,N);
Rwy = zeros(N,N);

for i=1:N           % Influence of control points
    for j=1:N       % On evaluation points
        theta_int = linspace((theta(i)-0.5*dtheta),(theta(i)+0.5*dtheta),1000);
        Rx_int = (-(xe(j)+sin(theta_int)).*sin(theta_int)+(ye(j)-cos(theta_int)).*cos(theta_int))./((xe(j)+sin(theta_int)).^2+(ye(j)-cos(theta_int)).^2);
        Ry_int = (-(xe(j)+sin(theta_int)).*cos(theta_int)-(ye(j)-cos(theta_int)).*sin(theta_int))./((xe(j)+sin(theta_int)).^2+(ye(j)-cos(theta_int)).^2);
        
        Rwx(i,j) = -trapz(theta_int,Rx_int);
        Rwy(i,j) = -trapz(theta_int,Ry_int);
    end
end


%% Compute induced velocities

% Prealocate
wx = zeros(1,N);
wy = zeros(1,N);
wxold = ones(1,N);
wyold = ones(1,N);
wx_i = zeros(N,N);
wy_i = zeros(N,N);

iter = 0;

while (abs(wx-wxold)>tol) | (abs(wy-wyold)>tol)
    % Initialise
    iter = iter+1;
    wxold = wx;
    wyold = wy;
    
    % Calculate velocities
    Vx = Vinf+omega*R.*cos(theta)+Vinf.*wx;    % X-direction, same as the wind speed
    Vy = omega*R.*sin(theta)+Vinf.*wy;         % Y-direction, pointing north
    Vt = Vx.*cos(theta)+Vy.*sin(theta);        % Tangential, from leading to trailing edge
    Vn = Vx.*sin(theta)-Vy.*cos(theta);        % Radial, pointing inwards
    Vrel = sqrt(Vx.^2+Vy.^2);                  % Relative velocity
    
    % Calculate angle of attack
    phi = atan(Vn./Vt);                        % Flow angle [rad]
    alpha = phi-theta_p-twist;                 % Angle of attack [rad]
    
    % Airfoil lift and drag
    [Cl, Cd] = liftdrag(alpha);
    
    % Blade loading
    Cn = Cl.*cos(alpha)+Cd.*sin(alpha);
    Ct = Cl.*sin(alpha)-Cd.*cos(alpha);
    
    Fn = 0.5*rho*Vrel.^2*c.*Cn;              % Normal, pointing outwards
    Ft = 0.5*rho*Vrel.^2*c.*Ct;              % Tangential, from trailing to leading edge

    Qn = cos(delta)*B*(Fn.*cos(theta_p)-Ft.*sin(theta_p))./(2*pi*R*rho*Vinf^2);
    Qz = sin(delta)*B*(Fn.*cos(theta_p)-Ft.*sin(theta_p))./(2*pi*R*rho*Vinf^2);
    Qt = -B*(Ft.*cos(theta_p)+Fn.*sin(theta_p))./(2*pi*R*rho*Vinf^2);
    

   % Influence of control point(i) on other locations(j)
    wx = 1/(2*pi)*Qn*Rwx + 1/(2*pi)*Qt*Rwy;
    wy = 1/(2*pi)*Qn*Rwy - 1/(2*pi)*Qt*Rwx;
    
    % Add effect of cylinder
    if f >= 1                            % In the wake
        for j = N/2+1:N
            wx(j) = wx(j)+Qn(j)-Qn(N+1-j)-Qt(j)*(yc(j)/sqrt(R^2-yc(j)^2))-Qt(N+1-j)*(yc(j)/sqrt(R^2-yc(j)^2));
        end    
    elseif f < 1                        % Inside the cylinder
        for j = 1:N/2
            wx(j) = wx(j)-Qn(j)-Qt(j)*(ye(j)/sqrt(1-ye(j)^2));
        end
        for j = N/2+1:N
            wx(j) = wx(j)-Qn(N+1-j)-Qt(N+1-j)*(ye(j)/sqrt(1-ye(j)^2));
        end
    end
    
    % Mod-Lin model
    CT = dtheta*sum(Qn.*sin(theta) + Qt.*cos(theta));
    
    %%% method 1
    a = 0.0892074*CT^3 +0.0544955 * CT^2 + 0.251163 * CT -0.0017077;
    
  %  if a<=0.15
        ka = abs(1/(1-a));
   % else
  %      ka = 1/(1-a)*(0.65+0.35*exp(-4.5*(a-0.15)));
  %  end
   
     wxnew = ka.*wx;
     wynew = ka.*wy;
     
     % Relaxation
      wx = ((1-relax).*wxold + relax.*wxnew);
      wy = ((1-relax).*wyold+relax.*wynew);
      
%       if imag(sum(wx))~=0
%           fprintf('complex wx')
%           iter
%       end
     
if iter>200
        warning('More than 200 iterations, code is stopped')
        break
    end

end


%% Power coefficient
Cp_ideal=trapz(theta,Qn.*Vn./(Vinf)); 
Cp_real=1/(2*pi)*trapz(theta,B.*((Ft.*cos(theta_p)+Fn.*sin(theta_p)).*omega.*R))/(0.5*rho*Vinf^3*2*R);

Cp_local = -Qt.*lambda;
Cp = (sum(Cp_local*dtheta));

CT = dtheta*sum(Qn.*sin(theta) + Qt.*cos(theta));
CT2 = trapz(theta,Qn.*sin(theta)+Qt.*cos(theta));
end
clc; clear all; close all;
tic
%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Actuator Cylinder %%%
%%%%%%%%%%%%%%%%%%%%%%%%%

%% Input parameters
% Blade
sigma = 0.16;               % Solidity [-]
B = 3;                      % Number of blades [-]
R = 15;                     % Radius [m]
c = sigma*2*R/B;            % Chord length [m]
twist = 0;                  % Twist angle [rad] [keep 0]

% Flow conditions
Vinf = 6;                   % Incoming wind speed, from west to east [m/s]
rho = 1.225;                % Air density [kg/m^3]
lambda = 3.1;               % Tip speed ratio [-]
theta_p = deg2rad(0);       % Pitch angle [rad]
omega = lambda*Vinf/R;      % Rotational speed, counter-clockwise(+) [rad/s]

% Additional parameters
N = 40;                     % Number of cells in cylinder [-]   [ > 36 ]
Ngrid = 16;                 % Number of cells in grid [-]       [dont change]
f = 1.001;                  % Factor for evaluation points [-]  [dont change]
tol = 1e-7;                 % Tolerance for while loop [-]      [dont change]
relax = 0.1;                % Relaxation paramter [-]           [dont change]

% Output
flowfield = 'On';


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

    Qn = B*(Fn.*cos(theta_p)-Ft.*sin(theta_p))./(2*pi*R*rho*Vinf^2);
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
Cp = (sum(Cp_local*dtheta))

CT = dtheta*sum(Qn.*sin(theta) + Qt.*cos(theta))
CT2 = trapz(theta,Qn.*sin(theta)+Qt.*cos(theta));
%% Figures
figure(1)
    subplot(2,3,1)
    grid on
    hold on
    plot(theta*180/pi,alpha*180/pi,'-x','Color',[22 61 90]/255)
    xlabel('theta [deg]')
    ylabel('alpha [deg]')
    
    subplot(2,3,2)
    grid on
    hold on
    plot(theta*180/pi,Vrel,'-x','Color',[22 61 90]/255)
    xlabel('theta [deg]')
    ylabel('Vrel [m/s]')
    
    subplot(2,3,3)
    grid on
    hold on
    plot(theta*180/pi,Qn,'-x','Color',[22 61 90]/255)
    xlabel('theta [deg]')
    ylabel('Qn [-]')
    
    subplot(2,3,4)
    grid on
    hold on
    plot(theta*180/pi,Qt,'-x','Color',[22 61 90]/255)
    xlabel('theta [deg]')
    ylabel('Qt [-]')    
    
    subplot(2,3,5)
    grid on
    hold on
    plot(theta*180/pi,wx,'-x','Color',[22 61 90]/255)
    xlabel('theta [deg]')
    ylabel('wx [-]')      
    
    subplot(2,3,6)
    grid on
    hold on
    plot(theta*180/pi,wy,'-x','Color',[22 61 90]/255)
    xlabel('theta [deg]')
    ylabel('wy [-]')     
    
    
    
%% Flow field 

if flowfield == 'On'

    % Define grid
    xg = repmat(linspace(-3.5*R,3.5*R,Ngrid),1,Ngrid);
    yg = repmat(linspace(-2*R,2*R,Ngrid),Ngrid,1);
    yg = yg(:)';

    % Influence matrix
    Rwx_grid = zeros(N,Ngrid^2);
    Rwy_grid = zeros(N,Ngrid^2);
    wx_i_grid = zeros(N,Ngrid^2);
    wy_i_grid = zeros(N,Ngrid^2);

    for i=1:N          % Influence of control points
        for j=1:Ngrid^2       % On evaluation points
            theta_int = linspace((theta(i)-0.5*dtheta),(theta(i)+0.5*dtheta),1000);
            Rx_int = (-(xg(j)+sin(theta_int)).*sin(theta_int)+(yg(j)-cos(theta_int)).*cos(theta_int))./((xg(j)+sin(theta_int)).^2+(yg(j)-cos(theta_int)).^2);
            Ry_int = (-(xg(j)+sin(theta_int)).*cos(theta_int)-(yg(j)-cos(theta_int)).*sin(theta_int))./((xg(j)+sin(theta_int)).^2+(yg(j)-cos(theta_int)).^2);

            Rwx_grid(i,j) = -trapz(theta_int,Rx_int);
            Rwy_grid(i,j) = -trapz(theta_int,Ry_int);
        end
    end

    % Influence of control point(i) on other locations(j)
    wx_grid = 1/(2*pi)* Qn*Rwx_grid;
    wy_grid = 1/(2*pi)* Qn*Rwy_grid;

    % Add effect of cylinder
    for j = 1:Ngrid^2
        if (xg(j)^2+yg(j)^2>R^2) & (xg(j)> 0) & (yg(j)>-R) & (yg(j)<R) % Outside the cylinder, in the wake
            index = N/2+find(abs(yc(N/2+1:N)-yg(j))==min(abs(yc(N/2+1:N)-yg(j))));    
            wx_grid(j) = wx_grid(j)+Qn(index)-Qn(N+1-index);
        elseif xg(j)^2+yg(j)^2 < R^2                                % Inside the cylinder
            index = find(abs(yc(1:N/2)-yg(j))==min(abs(yc(1:N/2)-yg(j))));
            wx_grid(j) = wx_grid(j)-Qn(index);
        end
    end

    % Lin-Mod model
    wxnew_grid = ka.*wx_grid;
    wynew_grid = ka.*wy_grid;

    % Including free stream wind speed
    Vx_grid = Vinf+wxnew_grid*Vinf;
    Vy_grid = wynew_grid*Vinf;

    % Plot velocity field
    starty1 = -2*R:0.2*R:2*R;
    startx1 = 3.5*R*ones(size(starty1));
    starty2 = -2*R:0.2*R:2*R;
    startx2 = -3.5*R*ones(size(starty2));

    figure(3)
        subplot(2,2,1)
        hold on
        grid on
        axis equal
        quiver(xg,yg,wxnew_grid,wynew_grid)
        plot(xc,yc)
        xlabel('x/R [-]')
        xlabel('y/R [-]')
        title('Induced velocity field - vector field')

        subplot(2,2,2)
        hold on
        grid on
        axis equal
        streamline(reshape(xg,[Ngrid,Ngrid])',reshape(yg,[Ngrid,Ngrid])',reshape(wxnew_grid,[Ngrid,Ngrid])',reshape(wynew_grid,[Ngrid,Ngrid])',startx1,starty1)
        plot(xc,yc)
        xlabel('x/R [-]')
        xlabel('y/R [-]')
        title('Induced velocity field - streamlines')

        subplot(2,2,3)
        hold on
        grid on
        axis equal
        quiver(xg,yg,Vx_grid,Vy_grid)
        plot(xc,yc)
        xlabel('x/R [-]')
        xlabel('y/R [-]')
        title('Real velocity field - vector field')

        subplot(2,2,4)
        hold on
        grid on
        axis equal
        streamline(reshape(xg,[Ngrid,Ngrid])',reshape(yg,[Ngrid,Ngrid])',reshape(Vx_grid,[Ngrid,Ngrid])',reshape(Vy_grid,[Ngrid,Ngrid])',startx2,starty2)
        plot(xc,yc)
        xlabel('x/R [-]')
        xlabel('y/R [-]')
        title('Real velocity field - streamlines')

hfig = figure();
set(hfig,'Position',[100 100 1500 420])
f1 = subplot('Position',[0.06 0.16 0.42 0.81]);
        hold on
        grid on
        axis equal
        f = quiver(xg,yg,Vx_grid,Vy_grid,'Color',[22 61 90]/255);
        set(f,'LineWidth',1,'Color',[22 61 90]/255)
        plot(R*cos(0:0.01:2*pi),R*sin(0:0.01:2*pi), 'LineWidth',3,'Color', [83 166 157]/255)
        xlabel('X-coordinate, x/R [-]')
        ylabel('Y-coordinate, y/R [-]')
        %title('Real velocity field - vector field')
        xlim(R*[-3.5 3.5])
        set(gca,'fontsize',16)

f2 = subplot('Position',[0.56 0.16 0.42 0.81]);
        hold on
        grid on
        axis equal
        h = streamline(reshape(xg,[Ngrid,Ngrid])',reshape(yg,[Ngrid,Ngrid])',reshape(Vx_grid,[Ngrid,Ngrid])',reshape(Vy_grid,[Ngrid,Ngrid])',startx2,starty2);
        set(h,'LineWidth',1,'Color',[22 61 90]/255)
        plot(R*cos(0:0.01:2*pi),R*sin(0:0.01:2*pi), 'LineWidth',3,'Color', [83 166 157]/255)
        xlabel('X-coordinate, x/R [-]')
        ylabel('Y-coordinate, y/R [-]')
        %title('Real velocity field - streamlines')
        xlim(R*[-3.5 3.5])
        set(gca,'fontsize',16)
end    

toc
   
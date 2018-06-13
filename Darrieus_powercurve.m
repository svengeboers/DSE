clc; clear all; close all;
tic

%% Input parameters

% Blade
B = 3;                      % Number of blades [-] Make sure N in actuatorcylinder2 can be divided by B
c = 1.1;                    % Chord length [m]
R = 15;                     % Radius [m]

% Flow conditions
omegamax = 3.2933333;       % Rotational speed, counter-clockwise(+) [rad/s]
lambda = 3.8;
Vinfs = 0:1:20;             % range of wind speeds [m/s]
Vrated = 13;                % rated wind speed [m/s]

%% Make list with powers

Power = [];
lambdas = [];

for Vinf = Vinfs
    Vinf
    if Vinf <= Vrated
        [Ptot,~] = fDarrieus3D2(c,lambda,B,Vinf);    
        Power = [Power;Ptot];
        lambdas = [lambdas;lambda];
    end
    
    if Vinf > Vrated
        lambda = omegamax*R/Vinf;      % Tip speed ratio [-]
        [Ptot,~] = fDarrieus3D2(c,lambda,B,Vinf);
        Power = [Power;Ptot];
        lambdas = [lambdas;lambda];
    end
end

%% Plot 

figure(1)
plot(Vinfs,Power)
xlabel('Vinf [m/s]')
ylabel('P [W]')
title('Chord 1.1m')

figure(2)
plot(Vinfs,lambdas)
xlabel('Vinf [m/s]')
ylabel('P [W]')
title('Chord 1.1m')

toc
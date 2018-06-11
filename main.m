clc; clear all; close all;
tic

%% Input variables
lambdas=2.5:0.1:5.5;        % Tip speed ratio [-]
sigmas=0.01:0.01:0.2;        % Solidity [-]
c3 = 3;                     % Number of blades [-]
c4 = 15;                    % Radius [m]
c5 = 6;                     % Incoming wind speed, from west to east [m/s]
c6 = deg2rad(0);            % Pitch angle [rad]
filename = '/Users/Emile/Desktop/DSE/Polars/NACA2321_Re2.200_M0.02_N9.0_360_M.txt';

%%
Cp = zeros(length(lambdas),length(sigmas));
CT = zeros(length(lambdas),length(sigmas));
i=0;


%% Loop
for c1 = lambdas  
    i = i + 1;
    j = 0;
    c1
    for c2 = sigmas
        j = j + 1;
        [cp,ct] = actuatorcylinder(c1,c2,c3,c4,c5,c6,filename); 
        Cp(i,j) = cp; 
        CT(i,j) = ct;
    end
end


%% Merge
[lambdas,sigmas] = meshgrid(lambdas,sigmas);

%% Print
[I_rowCp, I_colCp] = find(Cp == max(max(Cp)));
[I_rowCT, I_colCT] = find(CT == max(max(CT)));
fprintf('The max Cp is %.2f with a lambda of %.1f and a sigma of %.2f. \n', max(max(Cp)), lambdas(1,I_rowCp), sigmas(I_colCp,1));
fprintf('The max CT is %.2f with a lambda of %.1f and a sigma of %.2f. \n', max(max(CT)), lambdas(1,I_rowCT), sigmas(I_colCT,1));

%% Plot
figure('position', [100, 200, 1500, 600])
subplot(1,2,1)
contourf(lambdas',sigmas',Cp,0.0:0.01:1.0,'ShowText','on')
xlabel('TSR')
ylabel('Solidity')
title('Cp')

subplot(1,2,2)
contourf(lambdas',sigmas',CT,0.0:0.01:1.5,'ShowText','on')
xlabel('TSR')
ylabel('Solidity')
title('CT')

toc
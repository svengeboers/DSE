clc; clear all; close all;
tic

%% input variables

chords = 0.5:0.1:2.5;       % Chord length [m]
omegas = 1:0.1:2;           % Rotational speed, counter-clockwise(+) [rad/s]

B = 3;                      % Number of blades [-] Make sure N in actuatorcylinder2 can be divided by B
gamma = deg2rad(10);        % Sweep angle [rad] (change in azamuthal angle of the blade per meter increase in height)



%%
Cp = zeros(length(chords),length(omegas));
CT = zeros(length(chords),length(omegas));
i=0;


%% Loop
for c = chords  
    i = i + 1;
    j = 0;
    c
    for omega = omegas
        j = j + 1;
        [Cptot,CTtot] = fDarrieus3D(c,omega,B,gamma); 
        Cp(i,j) = Cptot; 
        CT(i,j) = CTtot;
    end
end

%% Merge
[chords,omegas] = meshgrid(chords,omegas);

%% Print
[I_rowCp, I_colCp] = find(Cp == max(max(Cp)));
[I_rowCT, I_colCT] = find(CT == max(max(CT)));
fprintf('The max Cp is %.2f with a chord of %.1f and a omega of %.2f. \n', max(max(Cp)), chords(1,I_rowCp), omegas(I_colCp,1));
fprintf('The max CT is %.2f with a chord of %.1f and a omega of %.2f. \n', max(max(CT)), chords(1,I_rowCT), omegas(I_colCT,1));

%% Plot
figure('position', [100, 200, 1500, 600])
subplot(1,2,1)
contourf(chords',omegas',Cp,0.0:0.01:1.0,'ShowText','on')
xlabel('chord [m]')
ylabel('rotational speed [rad/s]')
title('Cp')

subplot(1,2,2)
contourf(chords',omegas',CT,0.0:0.01:1.5,'ShowText','on')
xlabel('chord [m]')
ylabel('rotational speed [rad/s]')
title('CT')

toc
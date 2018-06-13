%%Main file for the structural calculations

clc; clear; close all

%% Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 100;                                %Number of lengthwise elements
 
%geometric parameters
%rotor geometry
h = 60;                                 %rotor height in m
D = 30;                                 %rotor diameter in m


%airfoil properties
airfoil.filename = 'NACA2321.txt';     %Airfoil file
airfoil.t = 6*10^-3;                 	%Thickness in m
airfoil.C = 1.1;                      	%Chord length in m
theta_p = 0;                            %pitch angle 

airfoil.spar1 = true;                       	%spar 1 option
airfoil.s1t = 4*10^-3;                       	%spar thickness in m
airfoil.s1c = 0.2;                            	%cordwise position 

airfoil.spar2 = true;                           %spar 2 option
airfoil.s2t = 4*10^-3;                        	%spar thickness
airfoil.s2c = 0.5;                          	%cordwise position 

%Airodynamic properties
lambda = 3.4;                           %tip speed ratio
Vw = 13;                             	%Wind speed

Qn = 2399;                            	%N/m, the maximum aerodynamic forces on the blade
Qt = -454;                           	%N/m

qy = -Qn;
qx = Qt;

%material properties: APA6
E = 16.0e9;                             %Pa, young's modulus 
rho_mat = 2096;                         %material density in kg/m3
yieldstrength = 1229;                    %MPa


%general properties
g = 9.81;                              	%gravitational acceleration
a = 0.0;                                %strut location 1 
b = 1-a;                                %strut location 2


%% Airfoil moments of area, include in function
[Ixx0, Iyy0, ~, area, x_stress, y_stress] = airfoilInertia(airfoil,0,1);
[Ixx, Iyy, Ixy, ~, ~, ~] = airfoilInertia(airfoil,theta_p,0);


%% Centrifugal force
Mtot = h*area*rho_mat;                  %total mass
Wtot = Mtot*g;                          %total weight

omega  = Vw*lambda/(0.5*D);            	%rotational speed rad/s
Fc     = Mtot*omega^2*(0.5*D)/h;        %N/m
qy     = qy+Fc;

%% call matrix method
%The first column indicates the distance from the start of the rotor
%(bottom up). The second column is the degre of freedom. 
% DoF 1 : positive y direction (Force in y)
% Dof 2 : positve angle about x (Moment about x)
% DoF 3 : positive x direction (Force in x)
% DoF 4 : positive angle about y (moment about y)
% The third column indicates the set value of the boundary condtion (0
% means no displacement) in bcs. In matrix R, the third column represents
% the reaction forces.

bcs = [[a*h, 1, 0];
       [a*h, 2, 0];
       [a*h, 3, 0];
       [a*h, 4, 0];
       [b*h, 1, 0];
       [b*h, 2, 0];
       [b*h, 3, 0];
       [b*h, 4, 0]];

[R,znodes] = MatrixMethod(h,theta_p,qy,qx,Ixx0,Iyy0,E,bcs);

%% compute internal bending moment
Mx = zeros(1,length(znodes));
My = zeros(1,length(znodes));

i = 0;
for z = znodes
    i = i+1;
    Vx(i) = -qy*z;
    Vy(i) = -qx*z;
    Mx(i) = -qy*0.5*z^2;
    My(i) = -qx*0.5*z^2;
    
    for ii = 1:length(R(:,1));
        if z >= R(ii,1)
            %check DoF
            if R(ii,2) == 1
                Mx(i)=Mx(i) - R(ii,3)*(z - R(ii,1));
                Vy(i)=Vy(i) - R(ii,3);
            elseif R(ii,2) == 2
                Mx(i)=Mx(i) + R(ii,3);
            elseif R(ii,2) == 3
                My(i)=My(i) - R(ii,3)*(z - R(ii,1));
                Vx(i)=Vx(i) - R(ii,3);
            else
                My(i)=My(i) - R(ii,3);
            end
        end
            
    end
    
    %compute the normal stress in points 1-4
 	sigma1(i) = ((My(i)-Mx(i)*(Ixy/Ixx))./(Iyy - (Ixy/Ixx)*Ixy))*x_stress(1) +...
        ((Mx(i) - My(i).*(Ixy/Iyy))/(Ixx - (Ixy/Iyy).*Ixy))*y_stress(1);
    sigma2(i) = ((My(i)-Mx(i)*(Ixy/Ixx))./(Iyy - (Ixy/Ixx)*Ixy))*x_stress(2) +...
        ((Mx(i) - My(i).*(Ixy/Iyy))/(Ixx - (Ixy/Iyy).*Ixy))*y_stress(2);
    sigma3(i) = ((My(i)-Mx(i)*(Ixy/Ixx))./(Iyy - (Ixy/Ixx)*Ixy))*x_stress(3) +...
        ((Mx(i) - My(i).*(Ixy/Iyy))/(Ixx - (Ixy/Iyy).*Ixy))*y_stress(2);
    sigma4(i) = ((My(i)-Mx(i)*(Ixy/Ixx))./(Iyy - (Ixy/Ixx)*Ixy))*x_stress(4) +...
        ((Mx(i) - My(i).*(Ixy/Iyy))/(Ixx - (Ixy/Iyy).*Ixy))*y_stress(2);
    tau(i) = sqrt(Vx(i)^2+Vy(i)^2)/area;
    sigma1(i) = sqrt(sigma1(i)^2+3*tau(i)^2);
    sigma2(i) = sqrt(sigma2(i)^2+3*tau(i)^2);
    sigma3(i) = sqrt(sigma3(i)^2+3*tau(i)^2);
    sigma4(i) = sqrt(sigma4(i)^2+3*tau(i)^2);
end

%Add the normal stress to points 1-4
sigma1(find(znodes>R(1,1) & znodes<R(5,1))) = sigma1(find(znodes>R(1,1) & znodes<R(5,1))) - ...
    Wtot*0.5/area;
sigma2(find(znodes>R(1,1) & znodes<R(5,1))) = sigma2(find(znodes>R(1,1) & znodes<R(5,1))) - ...
    Wtot*0.5/area;
sigma3(find(znodes>R(1,1) & znodes<R(5,1))) = sigma3(find(znodes>R(1,1) & znodes<R(5,1))) - ...
    Wtot*0.5/area;
sigma4(find(znodes>R(1,1) & znodes<R(5,1))) = sigma4(find(znodes>R(1,1) & znodes<R(5,1))) - ...
    Wtot*0.5/area;

[~,idxm1] = max(abs(sigma1));[~,idxm2] = max(abs(sigma2));
[~,idxm3] = max(abs(sigma3));[~,idxm4] = max(abs(sigma4));
sigmam1 = sigma1(idxm1);sigmam2 = sigma2(idxm2);
sigmam3 = sigma3(idxm3);sigmam4 = sigma4(idxm4);

margin = ((yieldstrength*10^6/max(abs([sigmam1,sigmam2,sigmam3,sigmam4]))) - 1)*100;


%% Print the maximum stress
fprintf('Maximum stresses in the structure: \n')
fprintf('Maximum stress in cross-sectional point 1: %f MPa \n',sigmam1*10^-6)
fprintf('\t at z = %f m \n',znodes(idxm1))
fprintf('Maximum stress in cross-sectional point 2: %f MPa \n',sigmam2*10^-6)
fprintf('\t at z = %f m \n',znodes(idxm2))
fprintf('Maximum stress in cross-sectional point 3: %f MPa \n',sigmam3*10^-6)
fprintf('\t at z = %f m \n',znodes(idxm3))
fprintf('Maximum stress in cross-sectional point 4: %f MPa \n',sigmam4*10^-6)
fprintf('\t at z = %f m \n \n',znodes(idxm4))
fprintf('Yield Strength: %f MPa \n \n',yieldstrength)
fprintf('Margin: %f percent \n',margin)


%% Plot internal bending moment
figure(2)
subplot(3,1,1)
plot(znodes,Mx)
ylabel('M_y [Nm]')
grid on
subplot(3,1,3)
plot(znodes,[sigma1;sigma2;sigma3;sigma4]*10^-6)
legend('location 1','location 2','location 3', 'location 4')
ylabel('\sigma [MPa]')
xlabel('z [m]')
grid on
subplot(3,1,2)
plot(znodes,My)
ylabel('M_x [Nm]')
grid on


clc
clear
close all

%____Inputs________________________________________________________________

dL = 100;                                                                   %Number of lengthwise elements
%geometric parameters
%Rotor
h = 55;                                                                     %Rotor height in m
D = 40;                                                                     %Diameter in m
Phase = 60;                                                                 %Phase shift angle in degrees
%airfoil
t = 3*10^-3;                                                                %Thickness in m
C = 1.3;                                                                    %Chord length in m
filename = 'NACAM23.txt';                                                   %Airfoil file
theta_p = 0;                                                                %pitch angle 
%Airodynamic shit
Rtip = 3;                                                                   %tip speed ratio
%general
Vw = 12;                                                                    %Wind speed
rho = 2702;                                                                 %material density in kg/m3
g = 9.81;                                                                   %gravitational acceleration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%____Node position_________________________________________________________

rmax = 0.5*D;                                                               %maximum radius
ymax = 0.5*h;                                                               %height of maximum radius
a = rmax/(ymax*(h-ymax));                                                   %shape parameter a (parabolic shape)
b = a*h;                                                                    %shape parameter b (parabolic shape)

for i = 1:dL                                                                
    y(i) = i*h/dL;                                                          %List of vertical positions
    r(i) = -a*y(i)^2+b*y(i);                                                %List of radial positions
end

Phase = Phase/360.*2*pi;                                                    %Convert phase shift to radians
if Phase == 0                                                               %list of tangential position in radians
    theta = zeros(1,length(y));
    for i = 1:dL
        x(i) = r(i)*cos(theta(i));                                         
        z(i) = r(i)*sin(theta(i)); 
    end
else
    for i = 1:dL
        theta(i) = i*Phase/dL;
        x(i) = r(i)*cos(theta(i));                                         
        z(i) = r(i)*sin(theta(i)); 
    end
end
plot3(x,z,y)                                                                %Plot :D
xlim([-h/2,h/2])
ylim([-h/2,h/2])

%__________________________________________________________________________
%____________Inertias______________________________________________________
%__________________________________________________________________________

delimiterIn = ' ';
headerlinesIn = 2;

A = importdata(filename,delimiterIn, headerlinesIn);                        %Opens aifoil profile

x1 = A.data(1:17,1)*C;                                                      %x position points on top
y1 = A.data(1:17,2)*C;                                                      %x position points on bottom
x2 = A.data(18:34,1)*C;                                                     %y position points on top
y2 = A.data(18:34,2)*C;                                                     %x position points on bottom

for i = 1:length(x1)-1
    eleLT(i) = sqrt((x1(i+1)-x1(i))^2+(y1(i+1)-y1(i))^2);                   %calculate length of elements on top
    elePxT(i) = x1(i) + (x1(i+1)-x1(i))/2;                                  %calculate x coordinate middle of elements on top
    elePyT(i) = y1(i) + (y1(i+1)-y1(i))/2;                                  %calculate y coordinate middle of elements on top
    eleLB(i) = sqrt((x2(i+1)-x2(i))^2+(y2(i+1)-y2(i))^2);                   %same for bottom
    elePxB(i) = x2(i) + (x2(i+1)-x2(i))/2;
    elePyB(i) = y2(i) + (y2(i+1)-y2(i))/2;    
end

TopAreas = t*eleLT;                                                         %Areas elements on top
BottomAreas = t*eleLB;                                                      %Areas elements on bottom

%____Centroid location_____________________________________________________

xbar = 0;
ybar = 0;
for i = 1:length(TopAreas)
    xbar = xbar + TopAreas(i)*elePxT(i) + BottomAreas(i)*elePxB(i);
    ybar = ybar + TopAreas(i)*elePyT(i) + BottomAreas(i)*elePyB(i);
end
xbar = xbar/(sum(TopAreas)+sum(BottomAreas));                               %centroid x position
ybar = ybar/(sum(TopAreas)+sum(BottomAreas));                               %centroid y position

% figure(1)
% plot(x1,y1,x2,y2,[0,C],[ybar,ybar],[xbar,xbar],[-0.3*C,0.3*C])                  %Plot :D
% ylim([-0.3*C,0.3*C])
% xlim([0,C])

elePxT = elePxT - xbar;                                                     %redefine mid positions with respect to centroid
elePyT = elePyT - ybar;                                 
elePxB = elePxB - xbar;                                  
elePyB = elePyB - ybar;     

% figure(2)
% plot(elePxT,elePyT,'.',elePxB,elePyB,'.')
% ylim([-0.3*C,0.3*C])
% xlim([-0.5*C,0.5*C])

%____Rotation______________________________________________________________
theta_p = 2*pi*theta_p/360;
for i = 1:length(elePxT)
   rt(i) = sqrt(elePxT(i)^2+elePyT(i)^2); 
   rb(i) = sqrt(elePxB(i)^2+elePyB(i)^2);  
   phit(i) = atan(elePyT(i)/elePxT(i))-theta_p;
   phib(i) = atan(elePyB(i)/elePxB(i))-theta_p;
   if elePxT(i)<0
       elePxT(i) = -rt(i)*cos(phit(i));
       elePyT(i) = -rt(i)*sin(phit(i));
   else
       elePxT(i) = rt(i)*cos(phit(i));
       elePyT(i) = rt(i)*sin(phit(i));
   end    
   if elePxB(i)<0
       elePxB(i) = -rb(i)*cos(phib(i));
       elePyB(i) = -rb(i)*sin(phib(i));
   else
       elePxB(i) = rb(i)*cos(phib(i));
       elePyB(i) = rb(i)*sin(phib(i));
   end
end

% figure(3)
% plot(elePxT,elePyT,'-',elePxB,elePyB,'-')
% ylim([-0.3*C,0.3*C])
% xlim([-0.5*C,0.5*C])


%____element inertias______________________________________________________

for i = 1:length(eleLT)
    IxeleT(i) = 1/12*eleLT(i)*t^3 + elePyT(i)^2*TopAreas(i);
    IxeleB(i) = 1/12*eleLB(i)*t^3 + elePyB(i)^2*BottomAreas(i);
    IyeleT(i) = 1/12*eleLT(i)^3*t + elePxT(i)^2*TopAreas(i);
    IyeleB(i) = 1/12*eleLB(i)^3*t + elePxB(i)^2*BottomAreas(i);
end

%____Total inertias________________________________________________________

Ixx = sum(IxeleT)+sum(IxeleB);
Iyy = sum(IyeleT)+sum(IyeleB);

for i = length(elePyT)
   IxyeleT(i) = elePxT(i)*elePyT(i)*TopAreas(i);
   IxyeleB(i) = elePxB(i)*elePyB(i)*BottomAreas(i);
end

Ixy = sum(IxyeleT) + sum(IxyeleB);
Aair = sum(TopAreas)+sum(BottomAreas);
%__________________________________________________________________________
%__________________________________________________________________________
%__________________________________________________________________________



%____Weight________________________________________________________________

                                      
for i = 2:length(y)                                                         %calculate section length between nodes
     L(i-1) = sqrt(sqrt((x(i)-x(i-1))^2+(y(i)-y(i-1))^2)^2+(z(i)-z(i-1))^2);
end 
for i = 1:length(y)                                                         %Calculate mass of segment around node in kg
    if i == 1
        M(i) = rho*Aair*L(i)*0.5;
    else
        if i == length(y)
            M(i) = rho*Aair*L(i-1)*0.5;
        else
            M(i) = rho*Aair*(L(i-1)*0.5+L(i)*0.5);
        end
    end
end
        

W = g*M;                                                                    %Calculate weight of segment around node in N
yW = y-W/(10000./dL);                                                       %Create array of segment weight underneat blade shape


%____Plot weight___________________________________________________________
hold on
plot3(x, z, yW,'r--')                                                       %plots the line of the weights
for i = 1:length(y)
   plot3([x(i),x(i)],[z(i),z(i)],[y(i),yW(i)],'r')
end

%____Centrifugal force_____________________________________________________

omega  = Vw*Rtip/(0.5*D);                                                   %Calculates rotational speed in rad/s
for i = 1:length(r)
    Fc(i) = M(i)*r(i)*omega^2;                                              %Calculates magnitude centrifugal force in N
    xFc(i) = (r(i)+Fc(i)/(10000./dL))*cos(theta(i));                        %Make array of position points of the centrifugal force
    zFc(i) = (r(i)+Fc(i)/(10000./dL))*sin(theta(i));
end

%____Plot centrifugal force________________________________________________

plot3(xFc, zFc, y, 'y--')
for i = 1:length(y)
   plot3([x(i),xFc(i)],[z(i),zFc(i)],[y(i),y(i)],'y')
end


% _________________________________________________________________________
% Simple beam aproximation_________________________________________________
% This uses the forces calculated before, but applies them to a straight beam
%__________________________________________________________________________

%____Reaction forces_______________________________________________________

Fy = 0.5*sum(W);                                                            %reaction force in axis direction
Frup = 0;
for i = 1:length(r)                                                         %
    Frup = Frup + (Fc(i)*y(i)+W(i)*r(i))/h;
end
Frupx = Frup*cos(mean(theta));
Frupz = Frup*sin(mean(theta));
Frlow = sum(Fc)-Frup;
Frlowx = Frlow*cos(mean(theta));
Frlowz = Frlow*sin(mean(theta));

plot3([x(length(y)),x(length(y))-Frupx/(10000/dL)],[z(length(y)),z(length(y))-Frupz/(10000/dL)],[y(length(y)),y(length(y))],'black:')
plot3([x(1),x(1)-Frlowx/(10000/dL)],[z(1),z(1)-Frlowz/(10000/dL)],[y(1),y(1)],'black:')

%____Internal Forces_______________________________________________________

j = -Fy;
jx = -Frlowx;
jz = -Frlowz; 
for i = 1:length(y)
    j = j + W(i);
    jx = jx + Fc(i)*cos(theta(i));
    jz = jz + Fc(i)*sin(theta(i));
    intfn(i) = j;
    intfvx(i) = jx;
    intfvz(i) = jz;
end

%____Stresses______________________________________________________________

tauxyV = intfvx/Aair;
tauzyV = intfvz/Aair;
sigmaN = intfn/Aair;
% 
% fig = plt.figure(4)
% ax = fig.gca(projection='3d')
% ax.scatter(tauxyV, tauzyV, sigmaN)
% #ax.legend()
% #plt.grid(True)
% #plt.title("Stress envalope")
% plt.show()
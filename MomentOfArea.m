clc
clear
close all

%____Inputs________________________________________________________________
t = 3*10^-3;                                                                %Thickness in m
C = 1.3;                                                                    %Chord length in m
filename = 'NACA23021.txt';                                                   %Airfoil file
theta_p = 15;                                                                %pitch angle

spar1 = true;                                                               %is there a spar???
s1t = 4*10^-3;                                                              %spar thickness in m
s1c = 0.2;                                                                  %cordwise position 
spar2 = false;                                                               %is there a second spar???
s2t = 4*10^-3;                                                              %spar thickness
s2c = 0.6;                                                                  %cordwise position 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

delimiterIn = ' ';
headerlinesIn = 2;

A = importdata(filename,delimiterIn, headerlinesIn);                        %Opens aifoil profile

x1 = A.data(1:18,1)*C;                                                      %x position points on top
y1 = A.data(1:18,2)*C;                                                      %x position points on bottom
x2 = A.data(18:35,1)*C;                                                     %y position points on top
y2 = A.data(18:35,2)*C;                                                     %x position points on bottom

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
%____Spars_________________________________________________________________

if spar1 == true
    s1x = C*s1c;                                                            %x location front spar
    for i = 1:length(x1)
       if s1x < x2(i)                                                       %determine top and bottom of spar
          i
          s1y(1) = y1(20-i);
          s1y(2) = y2(i-1);
          break
       end
       if s1c == x1(i)
           s1y(1) = y1(20-i);
           s1y(2) = y2(i-1);
           break
       end
    end
    s1height = s1y(1)-s1y(2);                                               %length front spar
    s1area = s1height*s1t;                                                  %area front spar
    s1mid = s1y(2)+0.5*s1height;                                            %midpoint
    figure(1)
    plot([s1x,s1x],s1y)                                                     %Plot :D
    ylim([-0.3*C,0.3*C])
    xlim([0,C])
end

if spar2 == true                                                            %same for rear spar
    s2x = C*s2c;
    for i = 1:length(x1)
       if s2c == x1(i)
           s2y(1) = y1(20-i);
           s2y(2) = y2(i);
           break
       end
       if s2x < x2(i)
          i
          s2y(1) = y1(20-i);
          s2y(2) = y2(i-1);
          break
       end
    end
    s2height = s2y(1)-s2y(2);
    s2area = s2height*s2t;
    s2mid = s2y(2)+0.5*s2height;
    hold on
    plot([s2x,s2x],s2y)                                                     %Plot :D
end

%____Centroid location_____________________________________________________

xbar = 0;
ybar = 0;
for i = 1:length(TopAreas)
    xbar = xbar + TopAreas(i)*elePxT(i) + BottomAreas(i)*elePxB(i);
    ybar = ybar + TopAreas(i)*elePyT(i) + BottomAreas(i)*elePyB(i);
end

if spar1 == true
    if spar2 == true
        xbar = xbar + s1area*s1mid  + s2area*s2mid;
        ybar = ybar + s1area*s1x  + s2area*s2x;
        Aair = sum(TopAreas)+sum(BottomAreas)+s1area+s2area;
    else
        xbar = xbar + s1area*s1mid;
        ybar = ybar + s1area*s1x;
        Aair = sum(TopAreas)+sum(BottomAreas)+s1area;
    end
else
    Aair = sum(TopAreas)+sum(BottomAreas);
end

xbar = xbar/Aair;                                                           %centroid x position
ybar = ybar/Aair;                                                           %centroid y position

if spar1 == true
    hold on
end
if spar1 == false
    figure(1)
    hold on
end
plot(x1,y1,x2,y2,[0,C],[ybar,ybar],[xbar,xbar],[-0.3*C,0.3*C])              %Plot :D
scatter(x1,y1)
scatter(x2,y2)
ylim([-0.3*C,0.3*C])
xlim([0,C])

elePxT = elePxT - xbar;                                                     %redefine mid positions with respect to centroid
elePyT = elePyT - ybar;                                 
elePxB = elePxB - xbar;                                  
elePyB = elePyB - ybar;     

% figure(2)
% plot(elePxT,elePyT,'.',elePxB,elePyB,'.')
% ylim([-0.3*C,0.3*C])
% xlim([-0.5*C,0.5*C])

%Rotation__________________________________________________________________
theta_p = 2*pi*theta_p/360;
if spar1 == true
   s1r = sqrt((s1mid-ybar)^2 + (s1x-xbar)^2);
   s1phi = atan((s1mid-ybar)/(s1x-xbar))-theta_p;
   if s1x-xbar < 0
       s1x = -s1r*cos(s1phi);
       s1mid = -s1r*sin(s1phi);
   else
       s1x = s1r*cos(s1phi);
       s1mid = s1r*sin(s1phi);
   end
   figure(3)
   scatter(s1x,s1mid)
   if spar2 == true
       s2r = sqrt((s2mid-ybar)^2 + (s2x-xbar)^2);
       s2phi = atan((s2mid-ybar)/(s2x-xbar))-theta_p;
          if s2x-xbar < 0
              s2x = -s2r*cos(s2phi);
              s2mid = -s2r*sin(s2phi);
          else
              s2x = s2r*cos(s2phi);
              s2mid = s2r*sin(s2phi);
          end
   hold on       
   scatter(s2x,s2mid)
   end
end

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

if spar1 == true
    hold on
else
    figure(3)
end
plot(elePxT,elePyT,'-',elePxB,elePyB,'-')
ylim([-0.3*C,0.3*C])
xlim([-0.5*C,0.5*C])


%____Element inertias______________________________________________________

for i = 1:length(eleLT)
    IxeleT(i) = 1/12*eleLT(i)*t^3 + elePyT(i)^2*TopAreas(i);
    IxeleB(i) = 1/12*eleLB(i)*t^3 + elePyB(i)^2*BottomAreas(i);
    IyeleT(i) = 1/12*eleLT(i)^3*t + elePxT(i)^2*TopAreas(i);
    IyeleB(i) = 1/12*eleLB(i)^3*t + elePxB(i)^2*BottomAreas(i);
end

%____Total inertias________________________________________________________
phi = 2*pi*(90-theta_p)/360;

if spar1 == true
    blub1x = s1height^2*cos(phi)^2 + s1t^2*sin(phi)^2;
    blub1y = s1t^2*cos(phi)^2 + s1height^2*sin(phi)^2;
    if spar2 == true
        sIx = s1t*s1height^3/12*blub1x + (s1mid-ybar)^2*s1area + s2t*s2height^3/12 + (s2mid-ybar)^2*s2area;
        sIy = s1t^3*s1height/12 + (s1mid-xbar)^2*s1area + s2t^3*s2height/12 + (s2mid-xbar)^2*s2area;
        Ixx = sum(IxeleT)+sum(IxeleB) + sIx;
        Iyy = sum(IyeleT)+sum(IyeleB) + sIy;
    else
        sIx = s1t*s1height^3/12*blub1x + (s1mid-ybar)^2*s1area;
        sIy = s1t^3*s1height/12 + (s1mid-xbar)^2*s1area;
        Ixx = sum(IxeleT)+sum(IxeleB) + sIx;
        Iyy = sum(IyeleT)+sum(IyeleB) + sIy;       
    end
else
    Ixx = sum(IxeleT)+sum(IxeleB);
    Iyy = sum(IyeleT)+sum(IyeleB);
end

for i = length(elePyT)
   IxyeleT(i) = elePxT(i)*elePyT(i)*TopAreas(i);
   IxyeleB(i) = elePxB(i)*elePyB(i)*BottomAreas(i);
end

Ixy = sum(IxyeleT) + sum(IxyeleB);

%____Outputs_______________________________________________________________

Ixx, Iyy, Ixy, Aair









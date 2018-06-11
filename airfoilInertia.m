function [Ixx, Iyy, Ixy, area] = airfoilInertia(airfoil,theta_p)
%INPUTS:
% theta_p               - Airfoil pitch angle [deg]
% elem.xT                - x positions of top airfoil elements                                                     
% elem.yT                - y positions of top airfoil elements               
% elem.xB                - x positions of bottom airfoil elements               
% elem.yB                - y positions of bottom airfoil elements   
%
%OUTPUTS:
% Ixx, Iyy, Ixy          - moments of inertia 
    
    %% Airfoil moments of area
    t = airfoil.t;
    C = airfoil.C;
    
    %Opens aifoil profile
    delimiterIn = ' ';
    headerlinesIn = 2;
    airfoil = importdata(airfoil.filename,delimiterIn, headerlinesIn);  


    x1 = airfoil.data(1:17,1)*C;                   	%x position points on top
    y1 = airfoil.data(1:17,2)*C;                    %y position points on top
    x2 = airfoil.data(18:34,1)*C;                   %x position points on bottom
    y2 = airfoil.data(18:34,2)*C;                   %y position points on bottom
    clear airfoil

    for i = 1:length(x1)-1
        LT(i) = sqrt((x1(i+1)-x1(i))^2+(y1(i+1)-y1(i))^2);        	%calculate length of elements on top
        xT(i) = x1(i) + (x1(i+1)-x1(i))/2;                        	%calculate x coordinate middle of elements on top
        yT(i) = y1(i) + (y1(i+1)-y1(i))/2;                          %calculate y coordinate middle of elements on top
        LB(i) = sqrt((x2(i+1)-x2(i))^2+(y2(i+1)-y2(i))^2);           %same for bottom
        xB(i) = x2(i) + (x2(i+1)-x2(i))/2;
        yB(i) = y2(i) + (y2(i+1)-y2(i))/2;    
    end

    TopAreas = t*LT;                     	%Areas elements on top
    BottomAreas = t*LB;                     	%Areas elements on bottom
    TopAreas = TopAreas; BottomAreas = BottomAreas;
    area = sum(TopAreas)+sum(BottomAreas);         

    %____Centroid location_____________________________________________________

    xbar = 0;
    ybar = 0;
    for i = 1:length(TopAreas)
        xbar = xbar + TopAreas(i)*xT(i) + BottomAreas(i)*xB(i);
        ybar = ybar + TopAreas(i)*yT(i) + BottomAreas(i)*yB(i);
    end
    xbar = xbar/(sum(TopAreas)+sum(BottomAreas));                               %centroid x position
    ybar = ybar/(sum(TopAreas)+sum(BottomAreas));                               %centroid y position

    % figure(1)
    % plot(x1,y1,x2,y2,[0,C],[ybar,ybar],[xbar,xbar],[-0.3*C,0.3*C])                  %Plot :D
    % ylim([-0.3*C,0.3*C])
    % xlim([0,C])

    xT = xT - xbar;                                                     %redefine mid positions with respect to centroid
    yT = yT - ybar;                                 
    xB = xB - xbar;                                  
    yB = yB - ybar;     

    
    % apply pitch angle rotation
    theta_p = deg2rad(theta_p);
    for i = 1:length(xT)
       rt(i) = sqrt(xT(i)^2+yT(i)^2); 
       rb(i) = sqrt(xB(i)^2+yB(i)^2);  
       phit(i) = atan(yT(i)/xT(i))-theta_p;
       phib(i) = atan(yB(i)/xB(i))-theta_p;
       if xT(i)<0
           xT(i) = -rt(i)*cos(phit(i));
           yT(i) = -rt(i)*sin(phit(i));
       else
           xT(i) = rt(i)*cos(phit(i));
           yT(i) = rt(i)*sin(phit(i));
       end    
       if xB(i)<0
           xB(i) = -rb(i)*cos(phib(i));
           yB(i) = -rb(i)*sin(phib(i));
       else
           xB(i) = rb(i)*cos(phib(i));
           yB(i) = rb(i)*sin(phib(i));
       end
    end

    %calculate element inertias
    for i = 1:length(LT)
        IxeleT(i) = 1/12*LT(i)*t^3 + yT(i)^2*TopAreas(i);
        IxeleB(i) = 1/12*LB(i)*t^3 + yB(i)^2*BottomAreas(i);
        IyeleT(i) = 1/12*LT(i)^3*t + xT(i)^2*TopAreas(i);
        IyeleB(i) = 1/12*LB(i)^3*t + xB(i)^2*BottomAreas(i);
    end

    %total inertias
    Ixx = sum(IxeleT)+sum(IxeleB);
    Iyy = sum(IyeleT)+sum(IyeleB);

    for i = length(yT)
       IxyeleT(i) = xT(i)*yT(i)*TopAreas(i);
       IxyeleB(i) = xB(i)*yB(i)*BottomAreas(i);
    end

    Ixy = sum(IxyeleT) + sum(IxyeleB);
end



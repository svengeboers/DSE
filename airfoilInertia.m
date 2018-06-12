function [Ixx, Iyy, Ixy, area, x_stress, y_stress] = airfoilInertia(airfoil,theta_p,plotting)
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
    
    spar1 = airfoil.spar1;                     	%second spar
    s1t = airfoil.s1t;                       	%spar thickness in m
    s1c = airfoil.s1c;                         	%cordwise position 
    spar2 = airfoil.spar2;                      %first spar
    s2t = airfoil.s2t;                        	%spar thickness
    s2c = airfoil.s2c;
    
    
    %Opens aifoil profile
    delimiterIn = ' ';
    headerlinesIn = 2;
    airfoil = importdata(airfoil.filename,delimiterIn, headerlinesIn);  


    x1 = airfoil.data(1:17,1)*C;                   	%x position points on top
    y1 = airfoil.data(1:17,2)*C;                    %y position points on top
    x2 = airfoil.data(18:34,1)*C;                   %x position points on bottom
    y2 = airfoil.data(18:34,2)*C;                   %y position points on bottom
    clear airfoil
    

    
    %% calculate points to compute stresses in 
    [~,idx] = max(y1);  x_stress(1)=x1(idx); y_stress(1)=y1(idx);%airfoil top
    [~,idx] = max(x1);  x_stress(2)=x1(idx); y_stress(2)=y1(idx);%airfoil TE
    [~,idx] = min(y2);  x_stress(3)=x2(idx); y_stress(3)=y2(idx);%airfoul bottom
    [~,idx] = min(x1);  x_stress(4)=x1(idx); y_stress(4)=y1(idx);%airfoul LE
    
%   plot airfoil
    if plotting
        figure
        plot(x_stress(1),y_stress(1),'bo')
        hold on
        plot(x_stress(2),y_stress(2),'ko')
        plot(x_stress(3),y_stress(3),'ro')
        plot(x_stress(4),y_stress(4),'go')
        legend('point 1','point 2','point 3','point 4')
        plot(x1,y1,'r')
        plot(x2,y2,'b')

    end
    
    %%
    for i = 1:length(x1)-1
        LT(i) = sqrt((x1(i+1)-x1(i))^2+(y1(i+1)-y1(i))^2);        	%calculate length of elements on top
        xT(i) = x1(i) + (x1(i+1)-x1(i))/2;                        	%calculate x coordinate middle of elements on top
        yT(i) = y1(i) + (y1(i+1)-y1(i))/2;                          %calculate y coordinate middle of elements on top
        LB(i) = sqrt((x2(i+1)-x2(i))^2+(y2(i+1)-y2(i))^2);          %same for bottom
        xB(i) = x2(i) + (x2(i+1)-x2(i))/2;
        yB(i) = y2(i) + (y2(i+1)-y2(i))/2;    
    end

    TopAreas = t*LT;                            %Areas elements on top
    BottomAreas = t*LB;                     	%Areas elements on bottom
    area = sum(TopAreas)+sum(BottomAreas);    
    
    
    %% spars
    if spar1 == true
        s1x = C*s1c;                           	%x location front spar
        for i = 1:length(x1)
           if s1x < x2(i)                      	%determine top and bottom of spar
              
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
        s1height = s1y(1)-s1y(2);           	%length front spar
        s1area = s1height*s1t;                	%area front spar
        s1mid = s1y(2)+0.5*s1height;          	%midpoint
        if plotting
            figure(1)
            hold on
            plot([s1x,s1x],s1y)                    	%Plot :D
            ylim([-0.3*C,0.3*C])
            xlim([0,C])
        end
    end

    if spar2 == true                           	%same for rear spar
        s2x = C*s2c;
        for i = 1:length(x1)
           if s2c == x1(i)
               s2y(1) = y1(20-i);
               s2y(2) = y2(i);
               break
           end
           if s2x < x2(i)
              
              s2y(1) = y1(20-i);
              s2y(2) = y2(i-1);
              break
           end
        end
        
        s2height = s2y(1)-s2y(2);
        s2area = s2height*s2t;
        s2mid = s2y(2)+0.5*s2height;
        if plotting
            hold on
            plot([s2x,s2x],s2y)    
        end
    end

    
    %% Centroid location
    xbar = 0;
    ybar = 0;
    for i = 1:length(TopAreas)
        xbar = xbar + TopAreas(i)*xT(i) + BottomAreas(i)*xB(i);
        ybar = ybar + TopAreas(i)*yT(i) + BottomAreas(i)*yB(i);
    end
    xbar = xbar/(sum(TopAreas)+sum(BottomAreas));                               %centroid x position
    ybar = ybar/(sum(TopAreas)+sum(BottomAreas));                               %centroid y position

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



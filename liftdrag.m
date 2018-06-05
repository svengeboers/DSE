function [Cl, Cd] = liftdrag(alpha)

% Cl = 2*pi*sin(alpha)*1.1;
% Cd = 0;%0.023;%*ones(1,length(alpha));

filename = '/Users/Emile/Desktop/DSE/Polars/NACA0015_RE265000.txt';
delimiterIn = ' ';
headerlinesIn = 12;

A = importdata(filename,delimiterIn, headerlinesIn);

alphas = deg2rad(A.data(:,1));
Cls = A.data(:,2);
Cds = A.data(:,3);

Cl = interp1(alphas',Cls,alpha);
Cd = interp1(alphas',Cds,alpha);

end

function [windspeed] = Weibull2(height2)

years=importdata('year.txt');
months=importdata('month.txt');
days = importdata('day.txt');
speeds_rott = importdata('speed rotterdam.txt'); 

%constants
h0 = 10;                          %reference height
z0 = 0.5;                         %roughness
height = linspace (1, 100, 1000); %discretisation of height. 

%%% Weibull fit CHANGE V TO WHATEVER YOU WANT TO MAKE THE WEIBULL FIT FOR

v=speeds_rott;

%%% PROCESS DATA

% Remove nil speed data (to avoid infeasible solutions in the following)
v(find(v==0)) = [];

% Extract the unique values occuring in the series
uniqueVals = unique(v);

% Get the number of unique values
nbUniqueVals = length(uniqueVals);

% Find the number of occurences of each unique wind speed value
for i=1:nbUniqueVals
    nbOcc = v(find(v==uniqueVals(i)));
    N(i) = length(nbOcc);
end

% Get the total number of measurements
nbMeas = sum(N);

% To take into account the measurement resolution
% (i.e., a measured wind speed of 2.0 m/s may actually correspond to a
% real wind speed of 2.05 or 1.98 m/s), compute the delta vector which
% contains the difference between two consecutive unique values
delta(1) = uniqueVals(1);
for i=2:(nbUniqueVals)
    delta(i) = uniqueVals(i) - uniqueVals(i-1);
end

% Get the frequency of occurence of each unique value
for i=1:nbUniqueVals
    prob(i) = N(i)/(nbMeas*delta(i));
end

% Get the cumulated frequency
freq = 0;
for i=1:nbUniqueVals
    freq = prob(i)*delta(i) + freq;
    cumFreq(i) = freq;
end
% Linearize distributions (see papers)
ln = log(uniqueVals);
lnln = log(-log(1-cumFreq));

% Check wether the vectors contain inifinite values, if so, remove them
test = isinf(lnln);
for i=1:nbUniqueVals
    if (test(i)==1)
        ln(i)= [];
        lnln(i)= [];
    end
end

% Extract the line parameters (y=ax+b) using the polyfit function
params = polyfit(ln,lnln',1);
a = params(1);
b = params(2);
y=a*ln+b;

% Extract the Weibull parameters c and k
k = a;
height2 = height2;
c = (log(height2/z0)/log(h0/z0))*exp(-b/a);

% Define the cumulative Weibull probability density function
% F(V) = 1-exp(-((v/c)^k)) = 1-exp(-a2), with a1 = v/c, a2 = (v/c)^k
a1 = uniqueVals/c;
a2 = a1.^k;
cumDensityFunc = 1-exp(-a2); 

% Define the Weibull probability density function
%f(v)=k/c*(v/c)^(k-1)*exp(-((v/c)^k))=k2*a3.*exp(-a2), 
% with  k2 = k/c, a3 = (v/c)^(k-1)
k1 = k-1;
a3 = a1.^k1;
k2 = k/c;
densityFunc = k2*a3.*exp(-a2); 

[maximumY,maxX] = max(densityFunc);
windspeed = uniqueVals(maxX);


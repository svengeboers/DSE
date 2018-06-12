%In this program the effects of fatigue will be assessed. The plan of
%attack is  as follows:
%1: Find the probability of a certain rotational speed using the Weibull
%distribution and the tip speed ratio and divide this into bins that can be
%related to a certain load.
%2: Sum up these loads in a specific way (explained later). This will give
%us the specific damage.
%3: Calculate the number of cycles the structure can hold
%4: Assess and repeat if necessary




%%1 Find the probability of a certain rotational speed using the Weibull
%%distribution.
syms v ;
a = 4.66*3.8;
b = 1.97;

f = b/a*(v./a)^(b-1)*exp(-(v./a).^b);

V = 0:0.1:40;
F = double(subs(f,v,V));
plot(V,F)

Prob = makedist('Weibull','a',4.66,'b',1.97);

%Define damages hereafter. A load that is caused by a wind speed between 5
%and 6 m/s might for example make the structure fail after 10 million
%cycles. Therefore a weight will be added to the loads that is proportional
%to the chances they occur.

%Sum = int(f,v,5,6)*n/Nf56 + int(f,v,7,10)*n/Nf710


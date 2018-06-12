function [x_loc, N] = discretize(length, locations, N)
    % 	Description: discretizes beam into beam elements and returns an array
    % 				 with the x-location of where each element begins.
    % 				 There will be more nodes than the value given to 
    %                correct for non-int division.
    % 	Inputs:
    % 		length      =  length of beam [m]
    % 		locations   =  array with x-coord of where nodes should be
    %                       placed as a requirement
    % 		N           =  total number of nodes along beam [int]
    % 	
    % 	Outputs:
    % 		x_loc   =  x-location of nodes
    % 		sec_loc =  index of first node in each section

%function code verification: test using test data. Assess whether the
%output nodes are spaced properly and whether at least N nodes are placed
%in total.
% clear; close all; clc
% %load test data
% load('testdiscretize.mat')

%make bcs a row vector for the loop later


if iscolumn(locations)
    locations = locations';
end

%add 0 to the vector as the first element and remove from bcs    
x_loc = [];
locations(find(locations==0)) = [];

%add L to the bcs vector if it is not initially present
if ~any(locations==length)
    locations = [locations length];
end

%temporarily decrease N nodes by 1 so adding L to x_loc later makes N

locations = sort(unique(locations));

x0 = 0;
ii = 0;
for bc = locations
    ii = ii + 1;

    l_ratio = (bc - x0) / length;
    nod_sec = ceil(N * l_ratio) + 1;
    
    %space the points, eclude end unless final node
    if bc == length
        x_sec = linspace(x0, bc, nod_sec-1);  
    else
        x_sec = linspace(x0, bc, nod_sec); x_sec(end)=[];       
    end
    
    %add x locations to vector xloc
    x_loc = [x_loc x_sec];
    
    %update the value x0
    x0 = bc;
end
N = size(x_loc); N = N(2);
end



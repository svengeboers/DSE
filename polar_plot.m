clc; clear all; close all;
tic

filename = 'NACA2321_Re2.200_M0.02_N9.0_360_M.txt';
delimiterIn = ' ';
headerlinesIn = 14;

A = importdata(filename, delimiterIn, headerlinesIn);

alphas = A.data(:,1);
Cls = A.data(:,2);
Cds = A.data(:,3);

figure(1)
    subplot(2,1,1)
    %grid on
    hold on
    plot(alphas,Cls,'-x','Color',[22 61 90]/255)
    xlabel('alpha [deg]')
    ylabel('Cl [-]')
    
    subplot(2,1,2)
    %grid on
    hold on
    plot(alphas,Cds,'-x','Color',[22 61 90]/255)
    xlabel('alpha [deg]')
    ylabel('Cd [-]')

toc
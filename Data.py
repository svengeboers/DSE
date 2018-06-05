#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun May 28 14:31:58 2017

@author: sgeboers
"""

#Note for self: maybe replace the outlier by the mean instead of removing it
import numpy as np
import matplotlib.pyplot as plt
from math import pi,e
import scipy as sc
from matplotlib import cm
TimeLength = 600 #s

def Outlier(component,data,m,mean,std,removed):    #m is amount of sigma from the mean, data is the data
    
    for i in range(len(component)):
        if component[i]>mean+3*std or component[i]<mean-3*std:
            data = np.delete(data,(i-removed),axis=0)
            removed = removed + 1
    return data

def TurbulenceIntensity(data,intervalwidth,samplingfrequency):
    lst = []
    for i in range(len(data)-intervalwidth*samplingfrequency):
        IntervalData = data[i:i+intervalwidth*samplingfrequency]
        TI = np.std(IntervalData)/np.mean(IntervalData)
        lst.append(TI)
    array = np.array(lst)
    return array

IntensityRatiolst = []
IntensityBefore = []
IntensityAfter = []
lstLoss = []
lst = [1,4,5,6,7,8,9,10]
#All the different WM's:
for number in lst:
    string = 'DATA/WM'+str(number)
    data= np.genfromtxt(string,delimiter = ',',skip_header=5,skip_footer=2,usecols = [1,2,3])
    
    #NaN removal
    data = data[~np.isnan(data).any(1)]
    X = data[:,0]
    Y = data[:,1]
    Z = data[:,2]
    
    #Outlier removal (3-sigma method for each component (then delete entire row (so x, y and z)))
    Xmean = np.mean(X)
    Ymean = np.mean(Y)
    Zmean = np.mean(Z)
    
    Xstd = np.std(X)
    Ystd = np.std(Y)
    Zstd = np.std(Z)
    
    #Outlier detection and removal for the three different parts:
    OriginalLength = len(data)
    data = Outlier(X,data,3,Xmean,Xstd,OriginalLength-len(data))
    data = Outlier(Y,data,3,Ymean,Ystd,OriginalLength-len(data))
    data = Outlier(Z,data,3,Zmean,Zstd,OriginalLength-len(data))
    #print OriginalLength - len(data)
    #Reimport data since the array changed dimensions
    X = data[:,0]
    Y = data[:,1]
    Z = data[:,2]
    
    #Basic calculations (other calculations will follow from these results):
    XY = np.sqrt(X**2+Y**2)
    XYZ = np.sqrt(X**2+Y**2+Z**2)
    VerticalInflow = np.arctan2(Z,XY)
    YawWall = np.arctan2(Y,X)
    MeanYaw = np.mean(YawWall) #This is the direction the wind turbine will be placed in (the mean wind direction so to say)
    Yaw = abs(YawWall-MeanYaw) #This is the actual yaw on the wind turbine, abs is technically not necessary
    turbulenceIntensity = TurbulenceIntensity(XYZ,TimeLength,4)
    turbulenceIntensity60 = TurbulenceIntensity(XYZ,60,4)
    Ratio = np.mean(turbulenceIntensity)/np.mean(turbulenceIntensity60)
    IntensityRatiolst.append(Ratio)
    if number == 7 or number == 9 or number ==8:
        IntensityBefore.append(np.mean(turbulenceIntensity))
    else:
        IntensityAfter.append(np.mean(turbulenceIntensity))
        
        
    
    #PLOTTING
    #Turbulence intensity
    plt.close('all')
    TurbulenceString = 'TurbulenceIntensityWM'+str(number) + '.png'
    xturbulenceintensity = np.linspace(0,2-TimeLength*4./3600,len(turbulenceIntensity))
    plt.plot(xturbulenceintensity,turbulenceIntensity)
    plt.xlabel('Time [h]')
    plt.ylabel('Turbulence intensity [-]')
    #plt.title('Turbulence intensity in function of time for WM' + str(number) + '\n with a time interval of ' + str(TimeLength)+' s',fontweight='bold')
    plt.minorticks_on()
    plt.grid(b=True, which='major',linewidth=2)
    plt.grid(b=True, which='minor') 
    plt.savefig(TurbulenceString)
    plt.close()
        
    #Rolled vertical inflow angle
    def Theta(data,intervalwidth,samplingfrequency):
         lst = []
         for i in range(len(data)-intervalwidth*samplingfrequency):
            IntervalData = data[i:i+intervalwidth*samplingfrequency]
            Mean = np.mean(IntervalData)
            lst.append(Mean)
            array = np.array(lst)
         return array
    
    ThetaString = 'RolledThetaWM'+str(number) + '.png'
    RolledTheta = Theta(VerticalInflow,600,4)
    plt.title('Theta')
    plt.plot(RolledTheta)
    plt.xlabel('Data point')
    plt.ylabel('Theta [rad]')
    plt.title('Rolled theta in function of every datapoint for WM'+str(number),fontweight='bold')
    plt.minorticks_on()
    plt.grid(b=True, which='major',linewidth=2)
    plt.grid(b=True, which='minor') 
    plt.savefig(ThetaString)
    plt.close()

    #ENERGY CALCULATIONS      
    #1   Using the moment method to estimate Weibull parameters
    MeanXY = np.mean(XY)
    StdXY = np.std(XY)
    gamma = lambda x: sc.special.gamma(x)
    Constant = (StdXY/MeanXY)**2
    kCalculation = lambda k: gamma(1+2/k)/((gamma(1+1/k))**2)-1-Constant

    k = sc.optimize.newton(kCalculation,2.0) #Use Rayleigh function as initial guess for k
    c = MeanXY/gamma(1+1/k) 
    
    #2   Defining Weibull
    WeibullV = np.linspace(0,14,500)
    V = np.linspace(3,14,500) #Velocity array for possible plotting
    
    def Weibull(V,k,c):
        Weibull = k/c*(V/c)**(k-1) * e**(-(V/c)**k)
        return Weibull
    #3   Plotting Weibull vs histograms
    WeibullvsHistogramString = 'WeibullvsHistogramWM'+str(number)+'.png'
    plt.plot(WeibullV,Weibull(WeibullV,k,c))
    binwidth = 0.5 #m/s
    plt.hist(XY,normed=True,bins=np.arange(0, max(XY) + binwidth, binwidth))
    plt.xlabel('Velocity [m/s]')
    plt.ylabel('Probability [-]')
    plt.title('Weibull vs histogram for WM'+str(number),fontweight='bold')
    plt.minorticks_on()
    plt.grid(b=True, which='major',linewidth=2)
    plt.grid(b=True, which='minor') 
    plt.savefig(WeibullvsHistogramString)
    plt.close()
    
    #4   Importing power curve
    PowerCurveArray = np.genfromtxt('DATA/PowerSpline',delimiter = ',')          
    PowerCurve = sc.interpolate.interp1d(PowerCurveArray[:,0],PowerCurveArray[:,1],kind='cubic')
    
    #5   Define product
    def Product(V):
        product = PowerCurve(V)*Weibull(V,k,c)
        return product
    
    Integration = sc.integrate.quad(Product,3.0,14,limit=200) 
    ExpectedPower = Integration[0]
    #Plotting product
    ProductString = 'ProductGraphWM'+str(number)+'.png'
    plt.plot(V,Product(V))
    plt.xlabel('Velocity [m/s]')
    plt.ylabel('Weighted Power [W]')
    plt.title('Product of the Weibull and power curve for WM'+str(number),fontweight='bold')
    plt.minorticks_on()
    plt.grid(b=True, which='major',linewidth=2)
    plt.grid(b=True, which='minor') 
    plt.savefig(ProductString)
    plt.close()
    
    #6   Calculating the loss due to yaw:
    YawLoss = np.mean((np.cos(Yaw))**3)
    lstLoss.append(YawLoss)
    #7   Checking for high/low turbulence to increase/decrease the energy output and calculating the energy output.
    AverageTurbulence = np.mean(turbulenceIntensity)
    factor = 0
    if AverageTurbulence < 0.14:
        factor = 0.96
    elif AverageTurbulence > 0.18:
        factor = 1.04
    else:
        factor = 1.0
    Energy = 2*ExpectedPower*factor*YawLoss #This is the 2-hour produced energy (according to calculations)
    #8   Printing out energy output for each WM:
    #OPTIMAL SPACING (note that I put a lot of things inside the for loop to keep this part seperated from other calculations)
    
    #1   Defining calculating functions as well as coth(x)
    
    def LevelToIntensity(NoiseLevelIndB):
        I0 = 10.**(-12)          #This is the treshold hearing intensity, matching 0 dB
        NoiseLevel = float(NoiseLevelIndB)
        Intensity = I0*10**(NoiseLevel/10)
        return Intensity
        
    def IntensityToLevel(Intensity):
        I0 = 10.**(-12)          #This is the treshold hearing intensity, matching 0 dB
        Intensity = Intensity
        NoiseLevelIndB = 10*np.log10(Intensity/I0)
        return NoiseLevelIndB
       
    coth = lambda x: np.cosh(x)/np.sinh(x)
    
    
    #2   Defining some constants:
    SLHighway10 = 53.5 #dB, this is the sound level of a highway at 10 m distance
    ResidentialDay = 45. #dB, sound level that is allowed during daytime in a residential area
    ResidentialNight = 35. #dB, sound level that is allowed at night in a residential area
    d1 = 10. #m, distance between the highway and the sound barrier
    
    #3   Creating data mesh    
    b = np.arange(0.1, 150, 0.5)
    d = np.arange(0.1, 150, 0.5)
    b, d = np.meshgrid(b, d)
    
    #4   Calculating maximum velocity and individual sound power
    Vmax = round(np.mean(XY)+3*np.std(XY),2)  #m/s
    IntensityTurbine40cm = lambda V: 4*10**(-6)*e**(0.2216*V)
    IntensityIndividualTurbine = IntensityTurbine40cm(Vmax)
    PowerIndividual = IntensityIndividualTurbine*pi*0.16*4
    SoundPowerHighway = LevelToIntensity(SLHighway10)*pi*d1**2*4
    
    #5   Calculating intensity and sound level
    Intensity = PowerIndividual/(4*b*d)*coth(d/b*pi)+SoundPowerHighway/(4*pi*(d+d1)**2)
    SL = IntensityToLevel(Intensity)
    
    #6   Plots contour curve    
    levels = [41.,47.] #Contour levels that will be shown (the legal limits for night and day respectively)
    ContourString = 'ContourPlotWM' + str(number) + '.png'
    fig = plt.figure()
    CS = plt.contourf(d, b, SL, levels,cmap=cm.Greys)
    cbar=plt.colorbar()
    cbar.set_label('Sound level [dB]', rotation=270)
    plt.xlabel('Distance [m]')
    plt.ylabel('Spacing [m]')
    #plt.title('Sound level in function of distance and spacing \n with a velocity of '+str(Vmax)+' m/s for WM'+str(number),fontweight='bold')
    plt.minorticks_on()
    plt.grid(b=True, which='major',linewidth=2)
    plt.grid(b=True, which='minor') 
    plt.savefig(ContourString)
    plt.close()
    
    #COST CALCULATIONS
    Length = 1250*1000
    #Spacing = 45
    InitialCost = 1428.83
    MaintenanceYear = 178.603
    Maintenance10Year = 893.017
    Energy2Hour = Energy/1000
    EnergyYearkWh = Energy2Hour/2*24*365
    Price = 0.18190 - 0.0619018 #sell profit - net cost
    FamilyUsage = 3300
    Production = Price*EnergyYearkWh
    Time = np.linspace(0,30,301)
    Profit = Time*Production-InitialCost
    
    for year in Time:
        if year % 3 == 0:
            Profit[int(year)*10::]=Profit[int(year)*10::]-MaintenanceYear
        if year % 10 == 0 and year != 0 and year % 20 != 0:
            Profit[int(year)*10::]=Profit[int(year)*10::]-Maintenance10Year
        if year % 20 == 0 and year != 0:
            Profit[int(year)*10::]=Profit[int(year)*10::]-InitialCost
    
    plt.close('all')
    CostString = 'NetRevenueWM'+str(number)+'.png'
    plt.plot(Time,Profit)
    plt.axhline(y=0,color = 'k')
    plt.xlabel('Time (years)')
    plt.ylabel('Net revenue (EUR)')
    plt.title('Net revenue of one turbine in function of the time',fontweight='bold')  
    #plt.grid(b=True, which='both')
    plt.minorticks_on()
    plt.grid(b=True, which='minor')
    plt.grid(b=True, which='major',linewidth=2) 
    plt.savefig(CostString)
    plt.close()
    #Total = EnergyYearkWh*Length/Spacing
    #Families = Total/FamilyUsage
    #print Families
    

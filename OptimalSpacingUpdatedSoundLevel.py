# -*- coding: utf-8 -*-
"""
Created on Tue Apr 25 13:34:46 2017

@author: Sven Geboers
"""

from math import pi,e
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

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
   
#Definine the mathematical function coth(x)
coth = lambda x: (e**(x)-e**(-x))/(e**(x)-e**(-x)) #np.cosh(x)/np.sinh(x)

#Closes all previous plots so that we don't have to click them away manually
plt.close('all')

#Defining some constants:
SLHighway10 = 53.5 #dB, this is the sound level of a highway at 10 m distance
d1 = 10. #m, distance between the highway and the sound barrier

#Creating data mesh    
b = np.arange(0.1, 150, 0.5)
d = np.arange(0.1, 150, 0.5)
b, d = np.meshgrid(b, d)

#Calculating maximum velocity and individual sound power
Vmax = 9.25  #m/s
IntensityTurbine40cm = lambda V: 4*10**(-6)*e**(0.2216*V)
IntensityIndividualTurbine = IntensityTurbine40cm(Vmax)
PowerIndividual = IntensityIndividualTurbine*pi*0.16 * 4
SoundPowerHighway = LevelToIntensity(SLHighway10)*pi*d1**2 * 4

#Calculating intensity and sound level
Intensity = PowerIndividual/(4*b*d)*coth(d/b*pi)+SoundPowerHighway/(4*pi*(d+d1)**2)
SL = IntensityToLevel(Intensity)

#Plots contour curve    
levels = [41.,47.] #Contour levels that will be shown
fig = plt.figure()
CS = plt.contourf(d, b, SL, levels,cmap=cm.Greys)
cbar=plt.colorbar()
cbar.set_label('Sound level in dB', rotation=270)
plt.xlabel('Distance (m)')
plt.ylabel('Spacing (m)')
plt.title('Sound level in function of distance and spacing \n with a velocity of 9.25 m/s for WM6',fontweight='bold')
plt.minorticks_on()
plt.grid(b=True, which='major',linewidth=2)
plt.grid(b=True, which='minor') 
plt.show()

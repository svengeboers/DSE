import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from math import pi
def Helix(z,H,D,alpha,phase):
    a = -2*D/(H**2)
    R = a*(z-H/2)**2+D/2
    x = R*np.cos(z*alpha+phase)
    y = R*np.sin(z*alpha+phase)
    
    # Calculate S in function of Z:
    dz = np.diff(z,n=1,axis=0)
    dx = np.diff(x,n=1,axis=0)
    dy = np.diff(y,n=1,axis=0)
    ds = np.sqrt(dz**2+dx**2+dy**2)
   
    S = np.cumsum(ds)
    S = np.insert(S,0,0) #insert first value of 0.0
    return x,y,z,S

#S = np.linspace(0,5,100)
#Array1=Helix(2.5,6.25,2,S,0)
z=np.linspace(0,60,1000)

Array1=Helix(z,60,30,pi/3/60,0)
Array2=Helix(z,60,30,pi/3/60,2*pi/3)
Array3=Helix(z,60,30,pi/3/60,4*pi/3)



Array = np.array(list(zip(z,Array1[3])))
np.savetxt('S.csv',Array,delimiter=',')

a = 15
b = 45*3/pi
c = np.sqrt(a**2+b**2)
s = np.linspace(0,c*pi/3,1000)
x = a*np.cos(s/c)
y = a*np.sin(s/c)
z = b*s/c

fig=plt.figure()
ax=fig.add_subplot(111,projection='3d')
plt.plot(x,y,z)
plt.show()

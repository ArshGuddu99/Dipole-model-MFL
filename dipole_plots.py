#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 25 15:45:01 2020

@author: arshad
"""

import sympy as sym
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('ggplot')

#%%
a=0.001 #Crack half width
b=0.0005 #Crack height
h=0.005 #Crack depth

z_0=0.005 #Measurement location
y_0=0.001 #Measurement location

x_1=-0.005 #Starting location
x_2=0.005 #Ending location

sigma=sym.Symbol('sigma')   #Stress density (area)
H_a=sym.Symbol('H_a')   #External magnetization
mu=1.25
x_0=sym.Symbol('x_0')
x=sym.Symbol('x')
z=sym.Symbol('z')
y1=h*(1-(x/a))
y2=h*(1+(x/a))
dH_1x=sym.Symbol('H_1x')
dH_2x=sym.Symbol('H_2x')
dH_1z=sym.Symbol('H_1z')
dH_2z=sym.Symbol('H_2z')

#%%
'define basic quantities'
sigma=(2.65/(2*sym.pi))*(1+h/(2*a))/(1+h/(2*a*mu))
dH_1x=(x_0-x)/((x-x_0)**2+(y2-y_0)**2+(z-z_0)**2)*sym.sqrt(1+(h/a)**2)
dH_2x=(x-x_0)/((x-x_0)**2+(y1-y_0)**2+(z-z_0)**2)*sym.sqrt(1+(h/a)**2)
dH_1z=(z_0-z)/((x-x_0)**2+(y2-y_0)**2+(z-z_0)**2)*sym.sqrt(1+(h/a)**2)
dH_2z=(z-z_0)/((x-x_0)**2+(y1-y_0)**2+(z-z_0)**2)*sym.sqrt(1+(h/a)**2)

#%%
'Set up the numerical integration'
nx=100 #number of steps in x
nz=100 #number of steps in z
dx=a/nx   #Integration step width in x
dz=b/nz   #Integration step width in z

no_x=20 #Number of discrete steps in the line of calculation
x_c=np.linspace(x_1,x_2,no_x)

H_x=[]
H_z=[]

for x0 in x_c:   
    xl_1=np.linspace(-a,0,nx)
    xl_2=np.linspace(0,a,nx)
    zl=np.linspace(0,b,nz)
    
    dH_1x0=dH_1x.subs(x_0,x0)
    dH_2x0=dH_2x.subs(x_0,x0)
    dH_1z0=dH_1z.subs(x_0,x0)
    dH_2z0=dH_2z.subs(x_0,x0)
    
    H_1xpp=[dH_1x0.subs(x,xl_1[i]) for i in range(nx)]
    H_1xp=np.sum(H_1xpp)*dx
    H_1xp=[H_1xp.subs(z,zl[i]) for i in range(nz)]
    H_1x=np.sum(H_1xp)*dz
    H_2xpp=[dH_2x0.subs(x,xl_2[i]) for i in range(nx)]
    H_2xp=np.sum(H_2xpp)*dx
    H_2xp=[H_2xp.subs(z,zl[i]) for i in range(nz)]
    H_2x=np.sum(H_2xp)*dz
    H_1zpp=[dH_1z0.subs(x,xl_1[i]) for i in range(nx)]
    H_1zp=np.sum(H_1zpp)*dx
    H_1zp=[H_1zp.subs(z,zl[i]) for i in range(nz)]
    H_1z=np.sum(H_1zp)*dz
    H_2zpp=[dH_2z0.subs(x,xl_2[i]) for i in range(nx)]
    H_2zp=np.sum(H_2zpp)*dx
    H_2zp=[H_2zp.subs(z,zl[i]) for i in range(nz)]
    H_2z=np.sum(H_2zp)*dz
    
    Hx=H_1x+H_2x
    Hz=H_1z+H_2z

    H_x.append(Hx)
    H_z.append(Hz)
    
#%%
'Plot'
fig=plt.figure(figsize=(10,10))
ax=fig.add_subplot()
ax.plot(x_c,H_x)
ax.plot(x_c,H_z)
plt.legend(('Tangential Component','Normal Component'))
ax.set_ylabel(r"$\frac{H}{H_a}$")
ax.set_xlabel('Distance along measurement axis')
plt.title("Modeled Magnetic field intensity response\n Triangular Crack 2mm wide(x), 0.5mm high (z) , 5mm deep (y), measured at an axis at y=1mm and z=5mm")
plt.show()
#%%
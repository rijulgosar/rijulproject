# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 19:05:38 2016

@author: rijulgosar
"""

'''
MICROBIAL GROWTH FOR BIOGAS PRODUCTION
THREE KINDS OF MICROBES-
1- HYDROLYTIC: Convert Feed to Glucose
2- ACETOGENIC: Convert Glucose to Acetic Acid
3- METHANOGENIC: Convert Acetic Acid to Methane

ASSUMPTIONS-
    - Feed and Glucose compositions are in ecxess, and do not vary much with progress of reaction
    -Concentration of acetic acid changes due to acetogenic microbe, which changes pH
    -Microbial growth follows Monod Kinetics
    -Biomass change with pH is modelled using an equation given in a paper by Tom, Wang, and Marshall
    -The microbes used are:
        1. Hydrolytic: Clostridia 
        2. Acetogens: Syntrophobacter wolinii
        3. Methanogens: Syntrophomonas wolfei
    -Growth of microbes stops at very low pH, i.e., when conc of acetic acid becomes too large
DATA
'''
import numpy
import scipy
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import matplotlib.animation as animation

#For hydrolytic microbe
uhmax=17.25         #max specific growth rate
kh1=1.39*10**(-9)   #constant Koh
kh2=9.19*10**(-7)   #constant Kh
Khs=5               #Monod's constant for hydrolytic microbe
yh=2                #yield coefficient for hydrolytic microbe
khd=3               #death constant
yh=2                #yield constant for hydrolytic microbe

#For acetogenic microbe
uamax=2.4           #max specific growth rate
ka1=3.09*10**(-7)   #constant Koh
ka2=2.42*10**(-7)   #constant Kh
Kas=4               #Monod's constant for acetogenic microbe
ya=1.67             #yield coefficient for acetogenic microbe
kad=2               #death constant
ya=0.3              #yield constant for acetogenic microbe

#For methanogenic microbe
ummax=0.208         #max specific growth rate
km1=7.05*10**(-10)  #constant Koh
km2=6.86*10**(-8)   #constant Kh
Kms=3               #Monod's constant for methanogenic microbe
ym=0.43             #yield coefficient for methanogenic microbe
kmd=0.5             #death constant
ym=1                #yield constant for methanogenic microbe

s=10    #substrate input-feed
ph=7    #initially
n=10    #grid points
v=10    #reactor volume, litres
'''
fig=plt.figure()
ax=fig.add_subplot(111)
'''
from transpose import transpose

t=numpy.linspace(0,1,10)

xainit=10
xminit=10

for t1 in t:
    
    #for hydrolytic microbe
    xh=numpy.exp(((uhmax*s)/(Khs+s)-khd)*t)
    
    #for acetogenic microbe
    def ODEA(xa,t):
        da=-(uamax*xa)/(1+(10**(-ph)/ka2)+ka1/(10**(-ph)))-kad*xa
        return da
    
    def ODEM(xm,t):
        dm=-(ummax*xm*xm)/(1+(10**(-ph)/km2)+km1/(10**(-ph)))-kmd*xm
        return dm
    
    def ODEM2(xm,t):
        dm=-(ummax*xm*xm)/(1+(10**(-ph)/km2)+km1/(10**(-ph)))-kmd*xm-kmd*xa
        return dm

    xa=odeint(ODEA,xainit,t)
    #xm=odeint(ODEM,10,t)
    
    ph=-numpy.log10(xa/v)
    s=-xh/yh
    
    print ph
    
    if xa<5:
        xm=odeint(ODEM,xminit,t)
    else:
        xm=odeint(ODEM2,xminit,t)
    
    print transpose (xh)
    print xa
    print xm

fig, ax = plt.subplots()
ax.set_ylim(0, 20)
ax.set_xlim(0, 10)
ax.grid()
xdata, ydata = [], []

def run(data):
    # update the data
    t,xa = data
    xdata.append(t)
    ydata.append(xa)
    xmin, xmax = ax.get_xlim()
'''   
    if t >= xmax:
        ax.set_xlim(xmin, 2*xmax)
        ax.figure.canvas.draw()
    line.set_data(xdata, ydata)

    return line,
'''

ani = animation.FuncAnimation(fig, run, ODEA(xa,t), blit=True, interval=1,
    repeat=False)
plt.show()

'''
ax.imshow(xa,cmap=plt.cm.RdBu_r)
plt.pause(0.2)
plt.show()
'''
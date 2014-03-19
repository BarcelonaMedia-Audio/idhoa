from constants import *
from functions import PlotOverSphere
import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

'''
* This is not original work, but more a mix of online resources.
'''

def SpherePlotting(title,var):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    x,y,z = PlotOverSphere(phiTest,thetaTest,var)

    ax.scatter(x,y,z)
    ax.set_title(title, fontsize=20)
    #ax.plot_wireframe(x,y,z, rstride=10, cstride=10)
    plt.show()


def SpeakersPlotting(phi,theta,rho):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    x,y,z = PlotOverSphere(phi,theta,rho)

    ax.scatter(x,y,z)
   # ax.plot_wireframe(x,y,z, rstride=10, cstride=10)
    plt.show()


def Polar(title,angle,*variables):
    plt.ion() 
    # radar black, solid grid lines
    plt.rc('grid', color='k', linewidth=1, linestyle='-')
    plt.rc('xtick', labelsize=15)
    plt.rc('ytick', labelsize=15)
    
    # force square figure and square axes looks better for polar, IMO
    width, height = plt.rcParams['figure.figsize']
    size = min(width, height)
    # make a square figure
    fig = plt.figure(figsize=(size, size))
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], polar=True, axisbg='1.0') # axisbg is the background colour
    # white canvas colour
    rect = fig.patch
    rect.set_facecolor('w')
    
    colorlist = ['k','r','g','c','m','y','b','w']
    stylelist = ['dashdot','dashed','solid']
    i=0
    for var in variables:
		if type(var)!=tuple: ax.plot(angle, var, color=colorlist[i], linestyle=stylelist[i], linewidth=4)
		if type(var)==tuple: leg=var
		i+=1
        
    ax.set_rmax(1.2)
    plt.grid(True)
    
    plt.legend( leg, loc = 'upper right', bbox_to_anchor = (1.125, 1.13))
    
    ax.set_title(title, fontsize=20)
    plt.show()
    fig.savefig(str(DEG)+"-"+str(DEC)+"-"+title+".eps")


def PlSpherePt(NP):
    # it is different from SpherePt because has some redundancy at 0 and 2*pi   
    # you can increase redundancy to ease the plot task...
    # (plots in vertical plane now are a bit messy)
    thetaPrev = [(np.pi/2.0-(float(i)/NP)*np.pi) for i in range(NP+1)]
    theta = []
    phi = []
    
    for i in range(len(thetaPrev)):
        n = max(int(2*NP*np.cos(thetaPrev[i])),1)
        phi.append( [(float(jj)/n*2)*np.pi for jj in range(n+1)] )
        temp = [thetaPrev[i]] *(n+1)
        theta.append(temp)

    phiok = [item for sublist in phi for item in sublist]
    thetaok = [item for sublist in theta for item in sublist]

    if len(phiok)!=len(thetaok): sys.exit("Died generating points on the sphere")
    return phiok, thetaok


def PlSpherePtRotat(NP):
    # it is different from SpherePt because has some redundancy at 0 and 2*pi   
    # you can increase redundancy to ease the plot task...
    # (plots in vertical plane now are a bit messy)
    thetaPrev = [((float(i)/(2*NP))*2.0*np.pi) for i in range(2*NP+1)]
    theta = []
    phi = []
    
    for i in range(len(thetaPrev)):
        n = max(int(2*NP*np.cos(thetaPrev[i])),1)
        phi.append( [(float(jj)/n*2)*np.pi for jj in range(n+1)] )
        temp = [thetaPrev[i]] *(n+1)
        theta.append(temp)

    phiok = [item for sublist in phi for item in sublist]
    thetaok = [item for sublist in theta for item in sublist]

    if len(phiok)!=len(thetaok): sys.exit("Died generating points on the sphere")
    return phiok, thetaok


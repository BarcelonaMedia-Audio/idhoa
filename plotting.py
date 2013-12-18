'''
 * This file is part of IDHOA software.
 *
 * Copyright (C) 2013 Barcelona Media - www.barcelonamedia.org
 * (Written by Davide Scaini <davide.scaini@barcelonamedia.org> for Barcelona Media)
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>,
 * or write to the Free Software Foundation, Inc.,
 * 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA
'''


def SpherePlotting(title,var):
    from constants import phiTest,thetaTest
    import matplotlib as cm
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    x,y,z = PlotOverSphere(phiTest,thetaTest,var)

    ax.scatter(x,y,z)
    ax.set_title(title, fontsize=20)
    #ax.plot_wireframe(x,y,z, rstride=10, cstride=10)
    plt.show()





def SpeakersPlotting(phi,theta,rho):
    from functions import PlotOverSphere
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    x,y,z = PlotOverSphere(phi,theta,rho)

    ax.scatter(x,y,z)
   # ax.plot_wireframe(x,y,z, rstride=10, cstride=10)
    plt.show()



def Polar(title,angle,*variables):
    from constants import DEC, DEG
    import matplotlib
    import pylab
    from matplotlib.pyplot import figure, show, rc, grid, legend, savefig
    
    pylab.ion()
    # radar black, solid grid lines
    rc('grid', color='k', linewidth=1, linestyle='-')
    rc('xtick', labelsize=15)
    rc('ytick', labelsize=15)
    
    # force square figure and square axes looks better for polar, IMO
    width, height = matplotlib.rcParams['figure.figsize']
    size = min(width, height)
    # make a square figure
    fig = figure(figsize=(size, size))
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], polar=True, axisbg='0.85') # axisbg is the background colour
    # white canvas colour
    rect = fig.patch
    rect.set_facecolor('w')
    
    colorlist = ['r','g','b','c','m','y','k','w']
    i=0
    for var in variables:
		if type(var)!=tuple: ax.plot(angle, var, color=colorlist[i], lw=3)
		if type(var)==tuple: leg=var
		i+=1
        
    ax.set_rmax(1.2)
    grid(True)
    
    legend( leg, loc = 'upper right', bbox_to_anchor = (1.125, 1.13))
    
    ax.set_title(title, fontsize=20)
    show()
    fig.savefig(str(DEG)+"-"+str(DEC)+"-"+title+".eps")


def PlSpherePt(NP):
    # it is different from SpherePt because has some redundancy at 0 and 2*pi   
    # you can increase redundancy to ease the plot task...
    # (plots in vertical plane now are a bit messy)
    import sys
    import numpy as np


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
    import sys
    import numpy as np


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




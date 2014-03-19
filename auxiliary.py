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
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program. If not, see <http://www.gnu.org/licenses/>,
* or write to the Free Software Foundation, Inc.,
* 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
'''


import numpy as np

## defining points over the sphere
def SpherePt(NP):

    thetaPrev = [(np.pi/2.0-(float(i)/NP)*np.pi) for i in range(NP)]
    theta = []
    phi = []
    
    for i in range(len(thetaPrev)):
        n = max(int(2*NP*np.cos(thetaPrev[i])),1)
        phi.append( [(float(jj)/n*2)*np.pi for jj in range(n)] )
        temp = [thetaPrev[i]] *(n)
        theta.append(temp)

    phiok = [item for sublist in phi for item in sublist]
    thetaok = [item for sublist in theta for item in sublist]

    if len(phiok)!=len(thetaok): sys.exit("Died generating points on the sphere")
    return np.asarray(phiok), np.asarray(thetaok)



#############
## weights ##
#############

def Wfront(phiT,thetaT):
    wvect = 1. + np.cos(phiT) * np.cos(thetaT)/2.
    return wvect

def Wplane(thetaT):
    wvect = 1. + np.cos(thetaT)**2.
    return wvect

def Wbinary(thetaT,tThresh):
    return  thetaT > tThresh


def angDist(phi1,e1,phi2,e2):
    import numpy as np
    # Implement haversine formula (see Wikipedia)
    dist = 2.0*np.arcsin(np.sqrt( (np.sin((e1 - e2)/2.0))**2 + np.cos(e1)*np.cos(e2)*(np.sin((phi1 - phi2)/2.0))**2 ) )
    return dist

def spkDist():
    from constants import PHI,THETA
    import numpy as np
    
    dvec = np.array([])
    mins = np.array([])
    for i in range(len(PHI)):
        for j in range(len(PHI)): # you could do probably something like: for j in range(i+1,len(PHI))
            dvec= np.append(dvec,[angDist(PHI[i],THETA[i],PHI[j],THETA[j])])

        mins = np.append(mins,[min(dvec[dvec!=0])])
        dvec = np.array([]) # reset dvec
    mean = np.mean(mins) # calculates the mean only of the smalles values - the closest speakers
    return mean



def autoremoval(PHI,THETA,phiT,thetaT):

    meanSpkDist = spkDist()
    phit = []
    thetat = []
    for i in range(len(phiT)):
        for j in range(len(PHI)):
            if angDist(phiT[i],thetaT[i],PHI[j],THETA[j]) < meanSpkDist*1.5:
                phit.append(phiT[i])
                thetat.append(thetaT[i])
                break
    return phit, thetat



def Wautoremoval(PHI,THETA,phiT,thetaT):

    meanSpkDist = spkDist()
    
    wvect = []
    for i in range(len(phiT)):
        for j in range(len(PHI)):
            if angDist(phiT[i],thetaT[i],PHI[j],THETA[j]) < meanSpkDist*1.5:
                temp = True
                break
            else: temp = False
        wvect.append(temp)

    return np.asarray(wvect)



def autoinit(PHI,THETA,SEED,AUTOREM,thetaThreshold):
    ###############################
    ## automatic initializations ##
    ###############################
    # Threshold for binary masking the lower region
    thetaThreshold *= np.pi/180.0

    phiTest,thetaTest = SpherePt(SEED) ## generating sampling space points

    if AUTOREM: phiTest, thetaTest = autoremoval(PHI,THETA,phiTest,thetaTest) # autoremoving test points without speakers around

    WfrontVec = Wfront(phiTest,thetaTest)
    WplaneVec = Wplane(thetaTest)
    WbinVec = Wbinary(thetaTest,thetaThreshold)
    WremVec = Wautoremoval(PHI,THETA,phiTest,thetaTest) # only necessary if WAUTOREM is true...

    NPOINTS = len(phiTest)

    return phiTest, thetaTest, NPOINTS, WfrontVec, WplaneVec, WbinVec, WremVec



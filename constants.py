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


import auxiliary as aux
import numpy as np

#######################
## Speakers position ##
#######################
# studio
PHI = [27.0415460692, 2.04283140615, 332.958453931, 1.03372198009, 90.5036008348, 269.496399165, 166.80995235 , 193.19004765 , 54.0379157019, 302.679984191, 124.309021931, 235.690978069, 92.4265791571, 267.573420843, 35.1695965277, 324.830403472, 15.3392440494, 344.660755951, 30.5463301933, 329.453669807, 143.681356818, 216.318643182, 301.738222232]          # position of the speakers
THETA = [-1.35939557882, -1.39079054354, -1.35939557882, 25.2255330043 , -2.05829450955, 2.05829450955 , 6.89524790513 , 6.89524790513 , -4.49039914903, -1.34181305608, 2.73330540595 , 2.73330540595 , 27.4212399195 , 27.4212399195 , 21.065228491  , 21.065228491  , -22.9266013259, -22.9266013259, 30.6450023912 , 30.6450023912 , 27.2300782225 , 27.2300782225 , 88.0711642016]        # position of the speakers

# icosahedron
#PHI = [45.0, -45.0,-135.0, 135.0, 45.0, -45.0,-135.0, 135.0, 90.0,  90.0, -90.0, -90.0, 0.0,   0.0, 180.0, 180.0, 69.1, -69.1, 110.9,-110.9]
#THETA = [35.3, 35.3, 35.3, 35.3,-35.3,-35.3,-35.3,-35.3, 69.1,-69.1, 69.1,-69.1, 20.9,-20.9, 20.9,-20.9, 0.0, 0.0,  0.0, 0.0]


################
## Parameters ##
################

DEC = 'maxRe' # available: basic, maxRe, phase
DEG = 1


# parameters
SEED = 17   # number of maximum horizontal points in the sampling function (divided by two)
            # from 14 to 20 it's a reasonable number

# FLAGS
MUTESPKRS = 0 # tries to mute "unnecessary" speakers
AUTOREM = 1 # removes some points from the sphere sampling (the farthest from the speakers)
PREFHOM = 0
# weights
WFRONT = 0
WPLANE = 0
WBIN = 0
WAUTOREM = 0 # the same as AUTOREM but with weights
thetaThreshold = -30 # degrees

print "You choose ", DEC, " at order: ",DEG


## number of speakers
NSPK = len(PHI)

# convert in radiants if necessary
if (np.asarray(PHI)>2.0*np.pi).any(): 
    ## conversion in radiants
    PHI = [PHI[i]*np.pi/180.0 for i in range(len(PHI))]
    THETA = [THETA[i]*np.pi/180.0 for i in range(len(THETA))]

phiTest, thetaTest, NPOINTS, WfrontVec, WplaneVec, WbinVec, WremVec = aux.autoinit(PHI,THETA,SEED,AUTOREM,thetaThreshold)


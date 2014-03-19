#!/usr/bin/python

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
import math as mh
import sys
import nlopt
import time
from constants import *
from functions import *
from plotting import *


#controls(NSPK,DEG) ## checks before running the main program
#SpeakersPlotting(PHI,THETA,1)
#SpeakersPlotting(phiTest,thetaTest,1)
#SpeakersPlotting(phiTest*WremVec,thetaTest*WremVec,1)

start = time.time()


## INITIAL PARAMETERS 
coeffDir = ambisoniC(phiTest,thetaTest,'basic',DEG,0) # calculating ambisonics coefficients in test directions
Guess0 = ambisoniC(PHI,THETA,DEC,DEG,0)     # initial guess
GuessPinv = ambisoniC(PHI,THETA,DEC,DEG,1)     # initial guess


f0 = function(Guess0,coeffDir)
fPV = function(GuessPinv,coeffDir)

if f0 > fPV and np.asarray([abs(GuessPinv)<1.]).all(): Guess0 = GuessPinv     # This way we choose the starting point closest to the minimum


sij = Sij(Guess0,coeffDir,NPOINTS)
pressure0, V0, energyD0, J0, Vlongit0, Jlongit0, Vtrans0, Jtrans0 = physDir(sij,phiTest,thetaTest)

####################################
## Initial PLOTTING
####################################
phiPl,thetaPl = PlSpherePtRotat(SEED)
Npt = len(phiPl)
phi = np.asarray(phiPl)
theta = np.asarray(thetaPl)

coeffPl = ambisoniC(phiPl,thetaPl,'basic',DEG,0)
sij = Sij(Guess0,coeffPl,Npt)
#sij = Sij(GuessPinv,coeffPl,Npt) # if you want to look at the pinv

pressureZ, VZ, energyDZ, JZ, VlongitZ, JlongitZ, VtransZ, JtransZ = physDir(sij,phiPl,thetaPl)

if DEC!='basic':
    Polar("Naive-Horizontal",phi[theta==0.],energyDZ[theta==0.],JlongitZ[theta==0.],JtransZ[theta==0.],('energy','intensity L','intensity T'))
if DEC=='basic':
    Polar("Naive-Horizontal",phi[theta==0.],pressureZ[theta==0.],VlongitZ[theta==0.],VtransZ[theta==0.],('pressure','velocity L','velocity T'))


if DEC!='basic':
    Polar("Naive-Vertical",theta[phi==0],energyDZ[phi==0],JlongitZ[phi==0],JtransZ[phi==0],('energy','intensity L','intensity T'))
if DEC=='basic':
    Polar("Naive-Vertical",theta[phi==0],pressureZ[phi==0],VlongitZ[phi==0],VtransZ[phi==0],('pressure','velocity L','velocity T'))




#####################################
## MINIMIZATION
#####################################

initvect = Guess0.reshape((1,(DEG+1)**2*NSPK))[0]
## Global Optimization
# GN_DIRECT, GN_DIRECT_L, GN_ORIG_DIRECT, GN_ORIG_DIRECT_L, GN_DIRECT_NOSCAL, GN_DIRECT_L_NOSCAL, GN_DIRECT_L_RAND_NOSCAL
# GN_CRS2_LM
# G_MLSL_LDS, G_MLSL
# GD_STOGO, GD_STOGO_RAND
# GN_ISRES
prog = []
run = 0
minstart = time.time()

while True:
    ## Local Optimization
    # LN_COBYLA, LN_BOBYQA, LN_NEWUOA, LN_NEWUOA_BOUND, LN_PRAXIS, LN_NELDERMEAD, LN_SBPLX
    if run > 0: initvect = ResCoeff.reshape((1,(DEG+1)**2*NSPK))[0]
    opt = nlopt.opt(nlopt.LN_SBPLX,len(initvect)) 
    opt.set_min_objective(function)
    tol = np.asarray([0.1]*len(initvect))
    #opt.add_equality_mconstraint(eqconstr, tol)    
    #opt.add_inequality_mconstraint(inconstr, tol)  
    opt.set_initial_step([0.2]*len(initvect))
    if run > 0: opt.set_initial_step([0.9]*len(initvect))
    #opt.set_initial_step(np.random.rand(len(initvect))); 


    upbound = [1.]*len(initvect)
    lowbound = [-1.]*len(initvect)
    
    nodes = ambisoniC(PHI,THETA,'basic',DEG,0)      # calculates the coefficients in the basic or naive scheme. 
    # These values are used to figure out which speakers lay in the nodes of some spherical harmonic
    for i in range((DEG+1)**2*NSPK):
        if np.asarray(abs(nodes)<3.*10e-4).reshape(1,((DEG+1)**2*NSPK))[0,i]: upbound[i] = 0.; lowbound[i] = 0.; # putting to zero the speakers that are in a node of a SH
        if run>0: 
            if np.asarray(abs(res)<3.*10e-4).reshape(1,((DEG+1)**2*NSPK))[0,i]: upbound[i] = 0.; lowbound[i] = 0.; # putting to zero the speakers that are in a node of a SH

    # number of non-zero elements 
    print len(np.nonzero(initvect)[0]), " non-zero elements"

    opt.set_upper_bounds(upbound)
    opt.set_lower_bounds(lowbound)

    opt.set_xtol_abs(10e-6)
    opt.set_ftol_abs(10e-6)

    if run > 0: 
        opt.set_xtol_abs(10e-8)
        opt.set_ftol_abs(10e-8)

    res = opt.optimize(initvect)
    result = opt.last_optimize_result() #alternative way to get the result
    
    ResCoeff = np.reshape(res, ((DEG+1)**2,NSPK))
    #print "Function value for Naive: ", function(Guess0,coeffDir), " for Pinv: ", function(GuessPinv,coeffDir), " for NL-min: ",function(ResCoeff,coeffDir), " Elapsed time: ", time.time()-minstart
    print "Function value for Naive: ", function(Guess0,coeffDir), " for Pinv: ", function(GuessPinv,coeffDir), " for NL-min: ",function(ResCoeff,coeffDir), " Elapsed time: ", time.time()-minstart
    minstart = time.time()
    prog.append(function(ResCoeff,coeffDir))
    
    #####################
    ## exit condition
    # if (run>0 and (str(prog[run-1])[0:6]==str(prog[run])[0:6] or prog[run]>min(prog)+1)): break # for PRAXIS use this
    if (run>0 and ("{:7.4f}".format(prog[run-1])=="{:7.4f}".format(prog[run]) or prog[run]>min(prog)+1) or not MUTESPKRS ): break # for SBPLX use this
    run+=1




#####################################
## RESULTS 
#####################################
print "Minimization Results:"
if result == nlopt.SUCCESS: print "Successful minimization"
if result == nlopt.STOPVAL_REACHED: print "Optimization stopped because stopval was reached"
if result == nlopt.FTOL_REACHED: print "Optimization stopped because ftol_rel or ftol_abs was reached"
if result == nlopt.XTOL_REACHED: print "Optimization stopped because xtol_rel or xtol_abs was reached"
if result == nlopt.MAXEVAL_REACHED: print "Optimization stopped because maxeval was reached"
if result == nlopt.MAXTIME_REACHED: print "Optimization stopped because maxtime was reached"
if result == nlopt.FAILURE: print "Minimization FAILED" 
if result == nlopt.INVALID_ARGS: print "Invalid arguments (e.g. lower bounds are bigger than upper bounds, an unknown algorithm was specified, etcetera)"
if result == nlopt.OUT_OF_MEMORY: print "Ran out of memory."
if result == nlopt.ROUNDOFF_LIMITED: print "Halted because roundoff errors limited progress. (In this case, the optimization still typically returns a useful result.)" 
if result == nlopt.FORCED_STOP: print "Halted because of a forced termination: the user called nlopt_force_stop(opt) on the optimization's nlopt_opt object opt from the user's objective function or constraints." 


ResCoeff = np.reshape(res, ((DEG+1)**2,NSPK))
print "\nCoefficients \n", ResCoeff.T
if ResCoeff.T[abs(ResCoeff.T>1.)].any(): print "WARNING: You reached a bad minimum."



##############################
## PLOTTING
##############################
## results plots
coeffPl = ambisoniC(phiPl,thetaPl,'basic',DEG,0)
sij = Sij(ResCoeff,coeffPl,Npt)
pressure, V, energyD, J, Vlongit, Jlongit, Vtrans, Jtrans = physDir(sij,phiPl,thetaPl)
if DEC!='basic':
    Polar("Horizontal",phi[theta==0.],energyD[theta==0.],Jlongit[theta==0.],Jtrans[theta==0.],('energy','intensity L','intensity T'))
if DEC=='basic':
    Polar("Horizontal",phi[theta==0.],pressure[theta==0.],Vlongit[theta==0.],Vtrans[theta==0.],('pressure','velocity L','velocity T'))


if DEC!='basic':
    Polar("Vertical",theta[phi==0],energyD[phi==0],Jlongit[phi==0],Jtrans[phi==0],('energy','intensity L','intensity T'))
if DEC=='basic':
    Polar("Vertical",theta[phi==0],pressure[phi==0],Vlongit[phi==0],Vtrans[phi==0],('pressure','velocity L','velocity T'))

SpeakersPlotting(phi,theta,Jlongit)



# output some files
fuma = np.asarray([ResCoeff[0]*np.sqrt(2), ResCoeff[1]*np.sqrt(3), ResCoeff[2]*np.sqrt(3), ResCoeff[3]*np.sqrt(3)])
Gfuma = np.asarray([Guess0[0]*np.sqrt(2), Guess0[1]*np.sqrt(3), Guess0[2]*np.sqrt(3), Guess0[3]*np.sqrt(3)])
np.savetxt("Guess0-Fuma.idhoa",np.asarray(Gfuma).T,fmt="%f",delimiter="  ")
np.savetxt("Fuma-1st.idhoa",np.asarray(fuma).T,fmt="%f",delimiter="  ")

filename = str(DEG)+"-"+str(DEC)+"-rem"+str(AUTOREM)
np.savetxt(filename+".idhoa",ResCoeff.T,fmt="%f",delimiter="  ")
np.savetxt(filename+"-Gpinv.idhoa",GuessPinv.T,fmt="%f",delimiter="  ")
np.savetxt(filename+"-G0.idhoa",Guess0.T,fmt="%f",delimiter="  ")
##
###############################

print "Elapsed time ", time.time()-start
print "Function value for Naive: ", function(Guess0,coeffDir), " for Pinv: ", function(GuessPinv,coeffDir), " for NL-min: ",function(ResCoeff,coeffDir)

wait = raw_input()

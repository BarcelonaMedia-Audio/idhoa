#!/usr/bin/python

## checks
def controls(NSPK,DEG):
    import sys
    if NSPK < (DEG+1)**2:
        sys.exit("\n##!!!!!!!!!!!!!!!!!!!!!!!!!!##\nYou don't have enough speakers for this decoding order.")


## defining points over the sphere
def SpherePt(NP):
    import sys
    import numpy as np

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




def gauss_legendre(n):
    import scipy as sp
    import scipy.linalg
    # http://www.scientificpython.net/1/post/2012/04/gausslegendre1.html
    k=sp.arange(1.0,n)       
    a_band = sp.zeros((2,n)) 
    a_band[1,0:n-1] = k/sp.sqrt(4*k*k-1) 
    x,V=sp.linalg.eig_banded(a_band,lower=True) 
    w=2*sp.real(sp.power(V[0,:],2)) 
    return x, w

def maxReG1(order):
    x,w = gauss_legendre(order+1)
    return max(x)

def maxReCoeffs(order):
    from scipy.special import legendre
    coeffs = [legendre(i)(maxReG1(order)) for i in range(order+1)]
    return coeffs


def ambisoniC(phi,theta,DEC,DEG,inversion):
    import numpy as np
    import math as mh
    import sys

   
    NUM = len(phi)         # number of speakers
  
    #things to be initializated
    g0 = [None]*(DEG+1)
    g1 = [None]*(DEG+1)
    W  = [None]*(len(phi))
    if DEG>=1:
        X  = [None]*(len(phi))
        Y  = [None]*(len(phi))
        Z  = [None]*(len(phi))
    if DEG>=2:
        V  = [None]*(len(phi))
        T  = [None]*(len(phi))
        R  = [None]*(len(phi))
        S  = [None]*(len(phi))
        U  = [None]*(len(phi))
    if DEG>=3:
        Q  = [None]*(len(phi))
        O  = [None]*(len(phi))
        M  = [None]*(len(phi))
        K  = [None]*(len(phi))
        L  = [None]*(len(phi))
        N  = [None]*(len(phi))
        P  = [None]*(len(phi))
    if DEG>=4:
        c16  = [None]*(len(phi))
        c17  = [None]*(len(phi))
        c18  = [None]*(len(phi))
        c19  = [None]*(len(phi))
        c20  = [None]*(len(phi))
        c21  = [None]*(len(phi))
        c22  = [None]*(len(phi))
        c23  = [None]*(len(phi))
        c24  = [None]*(len(phi))
    if DEG>=5:
        c25  = [None]*(len(phi))
        c26  = [None]*(len(phi))
        c27  = [None]*(len(phi))
        c28  = [None]*(len(phi))
        c29  = [None]*(len(phi))
        c30  = [None]*(len(phi))
        c31  = [None]*(len(phi))
        c32  = [None]*(len(phi))
        c33  = [None]*(len(phi))
        c34  = [None]*(len(phi))
        c35  = [None]*(len(phi))
    if DEG>=6:
        sys.exit("DEG =",DEG," is not implemented yet\n")


    #####################################################
    ## Calculating the decoding dependent coefficients ##
    #####################################################
    if DEC == 'basic':
        for i in range(0,DEG+1):
            g0[i] = 1.
            g1[i] = 1.

    elif DEC == 'maxRe':
        g1p = maxReCoeffs(DEG) 
        Egm = sum([(2.*i+1.)*g1p[i]**2 for i in range(DEG+1)])
        g0 = [np.sqrt(NUM/Egm)]*4
        g1 = [g1p[i]*g0[0] for i in range(len(g1p))]


    elif DEC == 'phase':
        g0d = np.sqrt(3. *NUM/4.)            # from Dani 1st order Ambisonics
        g1d = g0d*1./3.                    # from Dani 1st order Ambisonics
        g1p = [None]*(DEG+1)

        for i in range(0,DEG+1):
            g0[i] = np.sqrt(NUM*(2*DEG+1)) / (DEG+1)
            g1[i] = g0[i]*(mh.factorial(DEG)*mh.factorial(DEG+1)) / (mh.factorial(DEG+i+1)*mh.factorial(DEG-i))  # note that g1(0) is equal to g0
            g1p[i]= (mh.factorial(DEG)*mh.factorial(DEG+1)) / (mh.factorial(DEG+i+1)*mh.factorial(DEG-i))          # just for debugging purposes
            

    else:
        sys.exit('Decoding scheme unknow: ', DEC, ' Possible ones are: basic, maxRe, phase')



    ##########################################
    ## Calculating the "naive" coefficients ##
    ##########################################
    for i in range(0,len(phi)):
        if DEG>=0:
            W[i] = 1.0/NUM 
            # W has this sqrt(2) factor to take into account that W is not directly the pressure

        if DEG>=1: 
            # 1st order
            X[i] = np.sqrt(3)* np.cos(phi[i])*np.cos(theta[i])/NUM;    # factor sqrt(3) comes from sqrt(2m+1) - see tab.3.3 p.156 Jerome Daniel
            Y[i] = np.sqrt(3)* np.sin(phi[i])*np.cos(theta[i])/NUM;   # The decoding matrix is given for SN3D (aka 'semi-normalized') so to go to N3D
            Z[i] = np.sqrt(3)* np.sin(theta[i])/NUM;                # you have to use the alphas. 

        if DEG>=2:
            # 2nd order
            V[i] = np.sqrt(5.)* np.sqrt(3.)/2.*np.sin(2.*phi[i])*np.cos(theta[i])**2/NUM ;
            T[i] = np.sqrt(5.)* np.sqrt(3.)/2.*np.sin(phi[i])*np.sin(2.*theta[i])/NUM ;
            R[i] = np.sqrt(5.)* (3.*np.sin(theta[i])**2.-1.)/2./NUM ;   
            S[i] = np.sqrt(5.)* np.sqrt(3.)/2.*np.cos(phi[i])*np.sin(2.*theta[i])/NUM ;
            U[i] = np.sqrt(5.)* np.sqrt(3.)/2.*np.cos(2.*phi[i])*np.cos(theta[i])**2/NUM ; 

        if DEG>=3:
            # 3rd order
            Q[i] = np.sqrt(7.)* np.sqrt(5./8.)*np.sin(3.*phi[i])*np.cos(theta[i])**3/NUM;
            O[i] = np.sqrt(7.)* np.sqrt(15.)/2.*np.sin(2.*phi[i])*np.sin(theta[i])*np.cos(theta[i])**2/NUM; 
            M[i] = np.sqrt(7.)* np.sqrt(3./8.)*np.sin(phi[i])*np.cos(theta[i])*(5.*np.sin(theta[i])**2-1.)/NUM; 
            K[i] = np.sqrt(7.)* np.sin(theta[i])*(5.*np.sin(theta[i])**2-3.)/2./NUM;
            L[i] = np.sqrt(7.)* np.sqrt(3./8.)*np.cos(phi[i])*np.cos(theta[i])*(5.*np.sin(theta[i])**2-1.)/NUM; 
            N[i] = np.sqrt(7.)* np.sqrt(15.)/2.*np.cos(2.*phi[i])*np.sin(theta[i])*np.cos(theta[i])**2/NUM;  
            P[i] = np.sqrt(7.)* np.sqrt(5./8.)*np.cos(3.*phi[i])*np.cos(theta[i])**3/NUM;

        if DEG>=4:
            # 4th order
            c16[i] = np.sqrt(9.)* np.sqrt(35./2.)*3./8.* np.sin(4.*phi[i])*np.cos(theta[i])**4  /NUM;  # (4,-4)
            c17[i] = np.sqrt(9.)* np.sqrt(35.)*3./4.* np.sin(3.*phi[i])*np.sin(theta[i])*np.cos(theta[i])**3  /NUM;  # (4,-3)
            c18[i] = np.sqrt(9.)* np.sqrt(5./2.)/4.*( -3.* np.sin(2.*phi[i])*np.sin(theta[i])*np.cos(theta[i])**3 + 21.*np.sin(2.*phi[i])*np.sin(theta[i])**2*np.cos(theta[i])**2 ) /NUM;  # (4,-2)
            c19[i] = np.sqrt(9.)* np.sqrt(5.)/4.* ( 21.*np.sin(phi[i])*np.cos(theta[i])*np.sin(theta[i])**3 - 9.*np.sin(phi[i])*np.cos(theta[i])*np.sin(theta[i]) )  /NUM;  # (4,-1)
            c20[i] = np.sqrt(9.)* 3./64.*( 20.*np.cos(2.*theta[i]) - 35.*np.cos(4.*theta[i]) -9. )  /NUM;  # (4,-0)
            c21[i] = np.sqrt(9.)* np.sqrt(5.)/4.* ( 21.*np.cos(phi[i])*np.cos(theta[i])*np.sin(theta[i])**3 - 9.*np.cos(phi[i])*np.cos(theta[i])*np.sin(theta[i]) )  /NUM;  # (4, 1)
            c22[i] = np.sqrt(9.)* np.sqrt(5./2.)/4.*( -3.* np.cos(2.*phi[i])*np.sin(theta[i])*np.cos(theta[i])**3 + 21.*np.cos(2.*phi[i])*np.sin(theta[i])**2*np.cos(theta[i])**2 )  /NUM;  # (4, 2)
            c23[i] = np.sqrt(9.)* np.sqrt(35.)*3./4.* np.cos(3.*phi[i])*np.sin(theta[i])*np.cos(theta[i])**3  /NUM;  # (4, 3)
            c24[i] = np.sqrt(9.)* np.sqrt(35./2.)*3./8.* np.cos(4.*phi[i])*np.cos(theta[i])**4  /NUM;  # (4, 4)

        if DEG>=5:
            # 5th order
            c25[i] = np.sqrt(11.)* 3./16.*np.sqrt(77.)* np.sin(5.*phi[i])*np.cos(theta[i])**5  /NUM;  # (5,-5)
            c26[i] = np.sqrt(11.)* 3./8.*np.sqrt(385./2.)* np.sin(4.*phi[i])*np.sin(theta[i])*np.cos(theta[i])**4  /NUM;  # (5,-4) 
            c27[i] = np.sqrt(11.)* np.sqrt(385.)/16.* ( 9.*np.sin(3.*phi[i])*np.cos(theta[i])**3*np.sin(theta[i])**2 )  /NUM;  # (5, -3)
            c28[i] = np.sqrt(11.)* np.sqrt(1155./2.)/4.* ( 3.*np.sin(2.*phi[i])*np.cos(theta[i])**2*np.sin(theta[i])**3 - np.sin(2.*phi[i])*np.cos(theta[i])**2*np.sin(theta[i]) )  /NUM;  # (5,-2)
            c29[i] = np.sqrt(11.)* np.sqrt(165./2.)/8.* ( np.sin(phi[i])*np.cos(theta[i]) - 7.*np.sin(phi[i])*np.cos(theta[i])*np.sin(theta[i])**2 + 21.*np.sin(phi[i])*np.cos(theta[i])*np.sin(theta[i])**4 )  /NUM;  # (5,-1)
            c30[i] = np.sqrt(11.)* np.sqrt(11.)/8.* ( 15.*np.sin(theta[i]) - 70.*np.sin(theta[i])**3 + 63.*np.sin(theta[i])**5 ) /NUM;  # (5, 0)
            c31[i] = np.sqrt(11.)* np.sqrt(165./2.)/8.* ( np.cos(phi[i])*np.cos(theta[i]) - 7.*np.cos(phi[i])*np.cos(theta[i])*np.sin(theta[i])**2 + 21.*np.cos(phi[i])*np.cos(theta[i])*np.sin(theta[i])**4 )  /NUM;  # (5, 1)
            c32[i] = np.sqrt(11.)* np.sqrt(1155./2.)/4.* ( 3.*np.cos(2.*phi[i])*np.cos(theta[i])**2*np.sin(theta[i])**3 - np.cos(2.*phi[i])*np.cos(theta[i])**2*np.sin(theta[i]) )  /NUM;  # (5, 2)
            c33[i] = np.sqrt(11.)* np.sqrt(385.)/16.* ( 9.*np.cos(3.*phi[i])*np.cos(theta[i])**3*np.sin(theta[i])**2 )  /NUM;  # (5, 3)
            c34[i] = np.sqrt(11.)* 3./8.*np.sqrt(385./2.)* np.cos(4.*phi[i])*np.sin(theta[i])*np.cos(theta[i])**4   /NUM;  # (5, 4)
            c35[i] = np.sqrt(11.)* 3./16.*np.sqrt(77.)* np.cos(5.*phi[i])*np.cos(theta[i])**5  /NUM;  # (5, 5)

        if DEG>6:
            print "DEG =",DEG," is not implemented yet\n"

 




    if DEG==1:
        coeffs = np.array([W, Y, Z, X])
    elif DEG==2:
        coeffs = np.array([W, Y, Z, X, V, T, R, S, U])
    elif DEG==3:
        coeffs = np.array([W, Y, Z, X, V, T, R, S, U, Q, O, M, K, L, N, P])
    elif DEG==4:
        coeffs = np.array([W, Y, Z, X, V, T, R, S, U, Q, O, M, K, L, N, P, c16, c17, c18, c19, c20, c21, c22, c23, c24])
    elif DEG==5:
        coeffs = np.array([W, Y, Z, X, V, T, R, S, U, Q, O, M, K, L, N, P, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27, c28, c29, c30, c31, c32, c33, c34, c35])
 
    if inversion ==1:
        coeffs =  np.linalg.pinv(coeffs, rcond=1e-8).T /NUM
        coeffs[abs(coeffs)<1e-8] = 0. # because inversion gives some very small values somewhere
	



    ### MULTIPLYING FOR THE SELECTED DECODING SCHEME
    coeff = np.empty(coeffs.shape,dtype=np.float64)
    g1[0] = g0[0]
    for i in range(DEG+1):
        for jj in range(i**2,(i+1)**2):
            coeff[jj]=g1[i]*coeffs[jj]
    

    return coeff





##########################
## supporting functions ##
##########################

def Sij(coeffSpk, coeffDir,NSphPt):
    import numpy as np
    import sys
    from constants import NSPK

    sij = np.dot(coeffSpk.T, coeffDir*NSphPt) # this will have the dimensions of NSPK*NPOINTS
    if sij.shape != (NSPK,NSphPt): sys.exit("Wrong dimensions in Sij\n")
    return sij



def physOmni(Sij):
   
    # LOW FREQUENCIES
    # pressure
    pressure = sum(Sij)
    # velocity
    V = velocity(Sij)


    # HIGH FREQUENCIES
    # energy density
    energyD = sum(Sij*Sij)
    # intensity
    J = velocity(Sij*Sij)

    return pressure, V, energyD, J



def velocity(Sij):
    import numpy as np
    from constants import PHI,THETA
    phi = np.asarray(PHI)
    theta = np.asarray(THETA)
    Sx = np.cos(phi) * np.cos(theta)
    Sy = np.sin(phi) * np.cos(theta)
    Sz = np.sin(theta)
    # rewrite these five lines using a call to ambisoniC (first order)
    # and fixing the coefficient of X,Y,Z

    if Sx.shape[0] == Sij.shape[0] and Sy.shape[0] == Sij.shape[0] and Sz.shape[0] == Sij.shape[0]:

        # since Sij and Sx are numpy arrays, the * is the element-wise multiplication
        # http://wiki.scipy.org/NumPy_for_Matlab_Users#head-e9a492daa18afcd86e84e07cd2824a9b1b651935
        Vx = sum((Sij.T * Sx).T) / sum(Sij)
        Vy = sum((Sij.T * Sy).T) / sum(Sij)
        Vz = sum((Sij.T * Sz).T) / sum(Sij)

    else: 
        sys.exit("Something went wrong calculating velocity\n")

    return Vx, Vy, Vz



def Longit(A,B):
    Sum = 0
    if len(A)!=len(B): sys.exit("I'm dying in Longit function. Arrays with different shapes.")
    for i in range(len(A)):
        Sum += A[i]*B[i]

    return Sum


def Transv(A,B):    # aka Cross
    import numpy as np
    Sum = 0
    if len(A)!=len(B): sys.exit("I'm dying in Transv function. Arrays with different shapes.")
    for i in range(len(A)):
        Sum += (A[i]*B[i-1]-A[i-1]*B[i])**2.0

    return np.sqrt(Sum)


def physDir(Sij,phi,theta):
    import numpy as np
    phi = np.asarray(phi)
    theta = np.asarray(theta)
    Zx = np.cos(phi) * np.cos(theta)
    Zy = np.sin(phi) * np.cos(theta)
    Zz = np.sin(theta)
    # rewrite these five lines using a call to ambisoniC (first order)
    # and fixing the coefficient of X,Y,Z
    
    Z = Zx, Zy, Zz

    press, V, energyD, J = physOmni(Sij)
    
    Vlongitudinal = Longit(V,Z)
    Jlongitudinal = Longit(J,Z)

    Vtransverse = Transv(V,Z)
    Jtransverse = Transv(J,Z)

    return press, V, energyD, J, Vlongitudinal, Jlongitudinal, Vtransverse, Jtransverse



def oppGain(Sij):
    import numpy as np
    from constants import NSPK,NPOINTS
    # FOR IN-PHASE DECODING
    oppGain = np.zeros((NSPK,NPOINTS),dtype=np.float64)
    oppGain[Sij<0] = Sij[Sij<0]
    if oppGain.shape != (NSPK,NPOINTS): sys.exit("Wrong dimensions in oppGain\n")
    oppGain = sum(oppGain*oppGain)

    return oppGain




def sph2cart(phi,theta):
    import numpy as np
    """Acoustics convention!"""
    x = np.cos(phi) * np.cos(theta)
    y = np.sin(phi) * np.cos(theta)
    z = np.sin(theta)
    return x, y, z

def Physph2cart(phi,theta):
    import numpy as np
    """Physics convention!"""
    x = np.cos(phi) * np.sin(theta)
    y = np.sin(phi) * np.sin(theta)
    z = np.cos(theta)
    return x, y, z

def PlotOverSphere(phi,theta,rho):
    import numpy as np
    """Acoustics convention!"""
    x = np.cos(phi) * np.cos(theta)*rho
    y = np.sin(phi) * np.cos(theta)*rho
    z = np.sin(theta)*rho
    return x, y, z


def eval_grad(f, theta):
    import algopy
    theta = algopy.UTPM.init_jacobian(theta)
    return algopy.UTPM.extract_jacobian(f(theta))

def eval_hess(f, theta):
    import algopy
    theta = algopy.UTPM.init_hessian(theta)
    return algopy.UTPM.extract_hessian(len(theta), f(theta))


################################
# The function to be minimized #
################################
def function(VarCoeffSpk,grad):
    import numpy as np
    from constants import NSPK,DEG,DEC,PREFHOM,coeffDir,NPOINTS,phiTest,thetaTest,WFRONT,WPLANE,WBIN,WAUTOREM,WfrontVec,WplaneVec,WbinVec,WremVec
    from numpy.linalg import norm
    
    if grad.size>0:
        print "You have to implement the gradient or use another algorithm"

    """The function to be minimized"""
    CP = 400. # arbitrary coefficients
    CE = 400.
    CL = 100.
    CT = 100.
    CPH= 10000.
    CV = 25.
    Wj = np.ones(coeffDir.shape[1]) # biasing factor dependent on direction j
    Mj = np.ones(coeffDir.shape[1]) # biasing factor dependent on direction j

    if VarCoeffSpk.shape != ((DEG+1)**2,NSPK):
        VarCoeffSpk = np.reshape(VarCoeffSpk,((DEG+1)**2,NSPK)) 
    # this because the minimization library destroys the shape of VarCoeffSpk

    sij = Sij(VarCoeffSpk,coeffDir,NPOINTS)
    pressure, V, energyD, J, Vlongit, Jlongit, Vtrans, Jtrans = physDir(sij,phiTest,thetaTest)

    # weighting functions
    if WFRONT: Wj = Wj*WfrontVec # maybe it's better if you calculate Wfront Wplane and Wbinary outside the function, so that you calculate them only once...
    if WPLANE: Wj = Wj*WplaneVec
    if WBIN:   Wj = Wj*WbinVec
    if WAUTOREM: Mj = Mj*WremVec
    if WFRONT or WPLANE or WBIN: Wj = Wj*float(len(Wj))/sum(Wj)


    if DEC=='basic':
        Tpressure = ((1.-pressure)**2*Wj)/NPOINTS
        TVlon = ((1.-Vlongit)**2*Wj*Mj)/NPOINTS
        TVtrans = ((Vtrans)**2*Wj*Mj)/NPOINTS
        Tvar = 0
        if PREFHOM:
            Tvar = np.var(VarCoeffSpk[0])/(np.mean(VarCoeffSpk[0]))**2
        
        target = np.sum(CP*Tpressure + CL*TVlon + CT* TVtrans) + CV*Tvar 
        # one sum instead of three (in Tpressure, TVlon,...)
        # anyway norm is (much) faster...

    elif DEC=='maxRe':
        TenergyD = ((1.-energyD)**2*Wj)/NPOINTS
        TJlon = ((1.-Jlongit)**2*Wj*Mj)/NPOINTS
        TJtrans = ((Jtrans)**2*Wj*Mj)/NPOINTS
        Tvar = 0
        if PREFHOM:
            Tvar = np.var(VarCoeffSpk[0])/(np.mean(VarCoeffSpk[0]))**2
        
        target = np.sum(CE*TenergyD + CL*TJlon + CT*TJtrans) + CV*Tvar
        
        #TenergyD = norm((1.-energyD)*Wj)**2/NPOINTS
        #TJlon = norm((1.-Jlongit)*Wj*Mj)**2/NPOINTS
        #TJtrans = norm((Jtrans)*Wj*Mj)**2/NPOINTS
        #Tvar = 0
        #if PREFHOM:
        #    Tvar = np.var(VarCoeffSpk[0])/(np.mean(VarCoeffSpk[0]))**2
        #
        #target = CE*TenergyD + CL*TJlon + CT*TJtrans + CV*Tvar


    elif DEC=='phase':
        TenergyD = ((1.-energyD)**2*Wj)/NPOINTS
        TJlon = ((1.-Jlongit)**2*Wj*Mj)/NPOINTS
        TJtrans = ((Jtrans)**2*Wj*Mj)/NPOINTS
        Tvar = 0
        if PREFHOM:
            Tvar = np.var(VarCoeffSpk[0])/(np.mean(VarCoeffSpk[0]))**2
        ToppGain = norm(oppGain(sij))**2/NPOINTS

        target = np.sum(CE*TenergyD + CL*TJlon + CT*TJtrans) + CPH*ToppGain * CV*Tvar  # missing some extra factors (see dani's paper page 5)

    return target


## with nlopt you can also write linear and non linear constraints... have a look at the reference
## http://ab-initio.mit.edu/wiki/index.php/NLopt_Python_Reference
def eqconstr(result,x,grad):
	from constants import NSPK,PHI,THETA
	if grad.size > 0:
		print "gradient to be implemented"
		
		
	nodes = ambisoniC(PHI,THETA,'basic',DEG,0)	# calculates the coefficients in the basic or naive scheme. 
							# These values are used to figure out which speakers lay in the nodes of some spherical harmonic
	for i in range((DEG+1)**2*NSPK):
		if np.asarray(abs(nodes)<10e-8).reshape(1,((DEG+1)**2*NSPK))[0,i]: result[i] = x[i] # putting to zero the speakers that are in a node of a SH
	return


def inconstr(result,x,grad):
	from constants import NSPK,PHI,THETA
	if grad.size > 0:
		print "gradient to be implemented"
		
		
	nodes = ambisoniC(PHI,THETA,'basic',DEG,0)	# calculates the coefficients in the basic or naive scheme. 
							# These values are used to figure out which speakers lay in the nodes of some spherical harmonic
	for i in range((DEG+1)**2*NSPK):
		if np.asarray(nodes>0.).reshape(1,((DEG+1)**2*NSPK))[0,i]: result[i] > x[i]		# keeping the sign
		if np.asarray(nodes<0.).reshape(1,((DEG+1)**2*NSPK))[0,i]: result[i] < x[i]		# keeping the sign
	return



#############
## weights ##
#############

def Wfront():
    from constants import phiTest,thetaTest
    import numpy as np
    wvect = 1. + np.cos(phiTest) * np.cos(thetaTest)/2.
    return wvect

def Wplane():
    from constants import thetaTest
    import numpy as np
    wvect = 1. + np.cos(thetaTest)**2.
    return wvect

def Wbinary():
    from constants import thetaTest,thetaThreshold
    import numpy as np
    return  thetaTest > thetaThreshold


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

def autoremoval():
    from constants import phiTest,thetaTest,PHI,THETA

    meanSpkDist = spkDist()
    phit = []
    thetat = []
    for i in range(len(phiTest)):
        for j in range(len(PHI)):
            if angDist(phiTest[i],thetaTest[i],PHI[j],THETA[j]) < meanSpkDist*1.5:
                phit.append(phiTest[i])
                thetat.append(thetaTest[i])
                break
    return phit, thetat


def Wautoremoval():
    from constants import phiTest,thetaTest,PHI,THETA
    import numpy as np

    meanSpkDist = spkDist()
    
    wvect = []
    for i in range(len(phiTest)):
        for j in range(len(PHI)):
            if angDist(phiTest[i],thetaTest[i],PHI[j],THETA[j]) < meanSpkDist*1.5:
                temp = True
                break
            else: temp = False
        wvect.append(temp)

    return np.asarray(wvect)


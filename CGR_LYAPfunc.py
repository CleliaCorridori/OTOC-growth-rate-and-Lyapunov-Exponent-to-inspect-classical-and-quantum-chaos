from numpy import *
from scipy.optimize import curve_fit
import random as rand

# Here there is a set of functions to compute:
# 1- the OTOC growth rate, called CGR, for the quantum kicked rotator;
# 2- the OTOC growth rate, called CGR, for the classical kicked rotator;
# 3- the Lyapunv exponent computed numerically for the classical kicked rotator.

# 1. Quantum CGR numerically computed
#computation of the last kick to use to obtain lambda fitting the OTOC
def quantumCGR_num(C,kicks):
    '''Function to compute numerically the CGR of the OTOC for each Kick
    input:
    -C: the values of the OTOC
    -kicks: vector of kicks
    output:
    -quantum CGR'''
    quant_CGR=zeros(len(kicks)-1)
    for ii in range(1,len(kicks)):
        quant_CGR[ii-1] = mean(log(C[ii]/C[ii-1]))
    return(quant_CGR)

def CGR_fit_lin(x,a):
    ''' Function y=a to fit the quantum CGR function'''
    return(a)

def quantumCGR_fit(Kick, K, C_val):   
    ''' Function to set the intervals for the fit of the quantum CGR and to compute the fit.
    input:
    -Kick: vector of kics
    -K: kick strength
    -C_val: values of the OTOC
    -''' 
    CGR_fit=zeros(len(K))
    # intervals for the fit
    for kk in range(len(K)):
        if kk < 7:
            ii=25
            jj=-1
        if (kk==7 or kk==8):
            ii=16
            jj=23
        if (kk==9 or kk==10):
            ii=15
            jj=20
        if (kk <13 and kk>=11):
            ii=14
            jj=18
        if (kk <16 and kk>=13):
            ii=11
            jj=16
        if (kk==16 or kk==17):
            ii=8
            jj=13
        if (kk <20 and kk>=18):
            ii=5
            jj=12
        if kk==20:
            ii=5
            jj=9
        if kk==21:
            ii=5
            jj=8
        if kk==22:
            ii=4
            jj=8
        if (kk <26 and kk>=23):
            ii=3
            jj=7
        if (kk <29 and kk>=25):
            ii=3
            jj=5
        if (kk <33 and kk>=29):
            ii=3
            jj=4
        if (kk <38 and kk>=33):
            ii=2
            jj=3
        if (kk <44 and kk>=38):
            ii=1
            jj=2
        if (kk <len(K) and kk>=44):
            ii=0
            jj=2

        CGR_num=quantumCGR_num(C_val[:,kk],Kick) # numerical CGR
        CGR_fit[kk], cov = curve_fit(CGR_fit_lin, Kick[ii:jj], CGR_num[ii:jj]) #gitted CGR*2

    CGR_Q = CGR_fit/2
    return(CGR_Q)

############
# 2. Classical CGR numerically computed

def map_kr_CGR(x_i, p_i, dx_i, dp_i, Nkicks, K):
    ''' Function to compute the evolution of the classical KR, it returns only the square momentum displacement.
    input:
    -x_i: initial position
    -p_i: initial momentum
    -dx_i: initial displacement for the position
    -dp:i: initial displacement for the momentum
    -Nkicks: number of kikcs
    -K: kick strength
    output:
    final square displacement'''
    #initialize useful quantities
    p=zeros(Nkicks)
    x=zeros(Nkicks)
    dp=zeros(Nkicks)
    dx=zeros(Nkicks)
    # set the initial values
    x[0]=x_i
    p[0]=p_i
    dx[0]=dx_i
    dp[0]=dp_i
    # comupte the evolution
    for tt in range(Nkicks-1):
        # map of the KR
        p[tt+1]=(p[tt]+K*sin(x[tt]))%(2*pi)
        x[tt+1]=(x[tt]+p[tt+1])%(2*pi)
        # tangent map
        dp[tt+1] = (K*cos(x[tt])*dx[tt]+dp[tt]) # momentum displacement
        dx[tt+1] = (dx[tt]+dp[tt+1])            # position displacement

    return(dp**2)

############
# 3. Classical LE numerically computed

def map_kr(x_i, p_i, dx_i, dp_i, Nkicks, K):
    '''Function to compute the evolution of N steps of the classical KR
    input:
    -x_i: initial position
    -p_i: initial momentum
    -dx_i: initial displacement for the position
    -dp:i: initial displacement for the momentum
    -Nkicks: number of kikcs
    -K: kick strength
    output:
    -dx_f: final displacement for the position
    -dp_f: final displacement for the momentum
    -n_d: mean square value of the displacement
    -x[-1]: final position after Nkikcs
    -p[-1]: final momentum after Nkikcs'''
    #initialize useful quantities
    p=zeros(Nkicks+1)
    x=zeros(Nkicks+1)
    dp=zeros(Nkicks+1)
    dx=zeros(Nkicks+1)
    # set the initial values
    x[0]=x_i
    p[0]=p_i
    dx[0]=dx_i
    dp[0]=dp_i
    # comupte the evolution
    for tt in range(Nkicks):
        p[tt+1]=(p[tt]+K*sin(x[tt]))%(2*pi)
        x[tt+1]=(x[tt]+p[tt+1])%(2*pi)
        # tangent map
        dp[tt+1] = (K*cos(x[tt])*dx[tt]+dp[tt])
        dx[tt+1] = (dx[tt]+dp[tt+1])
    #norm of the displacement
    n_d = (sqrt(dx[-1]**2+dp[-1]**2))

    #final normalized displacements
    dx_f = dx[-1]/n_d
    dp_f = dp[-1]/n_d
    return(dx_f,dp_f,n_d,x[-1],p[-1])


def evoluz(P, x_i, p_i, dx_i, dp_i, Nkicks, K):
    '''Function to compute the evolution of P steps (eachone made of Nsteps) for the classical KR.
    input:
    -x_i: initial position
    -p_i: initial momentum
    -dx_i: initial displacement for the position
    -dp:i: initial displacement for the momentum
    -Nkicks: number of kikcs
    -K: kick strength
    output:
    n_d: mean square value of the displacement
    '''
    p= zeros(P)
    x= zeros(P)
    dp=zeros(P)
    dx=zeros(P)
    n_d=zeros(P)
    x[0]=x_i
    p[0]=p_i
    n_d[0]=sqrt(dx_i**2+dp_i**2)
    dx[0]=dx_i/n_d[0]
    dp[0]=dp_i/n_d[0]
    
    for pp in range(P-1):
        dx[pp+1], dp[pp+1], n_d[pp+1], x[pp+1], p[pp+1] = map_kr(x[pp], p[pp], dx[pp], dp[pp], Nkicks, K) 
        
    return(n_d)
    

def lyapFunc(P,Nkicks,nd):
    '''Function to compute the Lyapunov exponents for the classical KR'''
    lyap=(1/(Nkicks*P))*sum(log(nd))
    return(lyap)

    

############
# 4. Classical LE analitically computed
# -> is omitted here, it is simply log(K/2)
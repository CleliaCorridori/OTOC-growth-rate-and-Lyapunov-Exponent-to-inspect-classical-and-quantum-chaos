# Functions to compute OTOC v and Lyapunov exponent
from numpy import *
import math
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter
from scipy.signal import medfilt
import random as rand
############
# 1. Quantum CGR numerically computed
#computation of the last kick to use to obtain lambda fitting the OTOC
def lyapunov_num(C,kicks):
    '''Function to compute numerically the Lyapunov exponent for each Kick from OTOC'''
    lyap=zeros(len(kicks)-1)
    for ii in range(1,len(kicks)):
        lyap[ii-1] = log(C[ii]/C[ii-1])
    return(lyap)

def lyap_fit_lin(x,a):
    ''' Function y=a to fit'''
    return(a)

def lyapQ_comp(Kick, K,mean_C):   
    ''' Function to set the intervals for the fit and to compute the fit''' 
    lyap=zeros(len(K))
    for kk in range(len(K)):
        if kk < 9:
            ii=25
            jj=-1
        if (kk <16 and kk>=9):
            ii=10
            jj=15
        if (kk <20 and kk>=16):
            ii=7
            jj=10
        if (kk <23 and kk>=20):
            ii=5
            jj=8
        if (kk <26 and kk>=23):
            ii=3
            jj=6
        if (kk <29 and kk>=26):
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

        lyap_num=lyapunov_num(mean_C[:,kk],Kick)
        lyap[kk], cov = curve_fit(lyap_fit_lin, Kick[ii:jj], lyap_num[ii:jj])

    lyap_Q = lyap/2
    return(lyap_Q)

############
# 2. Classical CGR numerically computed

def map_kr_CGR(x_i, p_i, dx_i, dp_i, Nkicks, K):
    ''' Function to compute the evolution of the classical KR, this return only dp'''
    p=zeros(Nkicks)
    x=zeros(Nkicks)
    dp=zeros(Nkicks)
    dx=zeros(Nkicks)
    
    x[0]=x_i
    p[0]=p_i
    dx[0]=dx_i
    dp[0]=dp_i
    
    for tt in range(Nkicks-1):
        p[tt+1]=(p[tt]+K*sin(x[tt]))%(2*pi)
        x[tt+1]=(x[tt]+p[tt+1])%(2*pi)

        # tangent map
        dp[tt+1] = (K*cos(x[tt])*dx[tt]+dp[tt])
        dx[tt+1] = (dx[tt]+dp[tt+1])

    return(dp**2)

############
# 3. Classical LE numerically computed
#map function and tangent map
def map_kr(x_i, p_i, dx_i, dp_i, Nkicks, K):
    '''Function to compute the evolution of N steps of the classical KR'''
    p=zeros(Nkicks)
    x=zeros(Nkicks)
    dp=zeros(Nkicks)
    dx=zeros(Nkicks)
    
    x[0]=x_i
    p[0]=p_i
    dx[0]=dx_i
    dp[0]=dp_i
    
    for tt in range(Nkicks-1):
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
    '''Function to compute the evolution of P steps (made of Nsteps) of the classical KR'''
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
# 3. Classical LE analitically computed
# -> is omitted here, it is simply loj(K/2)
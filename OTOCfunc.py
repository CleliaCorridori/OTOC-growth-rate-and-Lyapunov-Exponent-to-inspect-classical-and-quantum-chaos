from numpy import *
import math
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter
from scipy.signal import medfilt
import random as rand

def nextpow2(x):
    #computes the exponent of the first power of 2 >=x
    return 1 if x == 0 else math.ceil(math.log2(x))

def save(mean_C,name):
    with open(name, 'w') as out:
        for line in mean_C:
            savetxt(out, line)

#####################################################################################################################
#####################################################################################################################
#####################################################################################################################

def OTOC(p_state,t,dp,N,K,T):
    #Computes the OTOC on the state p_state at the time t with P(0)=dp

    #U*=Udagger -> exp sign + -> sign=1
    #First term: <i| P U* P^2 U P |i>
    A=sum(conj(p_state)*dp*Ut(dp**2*Ut(dp*p_state,2,t,dp,N,K,T),1,t,dp,N,K,T))

    #Second term: <i| U* P U P^2 U* P U |i>
    B=sum(conj(p_state)*Ut(dp*Ut(dp**2*Ut(dp*Ut(p_state,2,t,dp,N,K,T),1,t,dp,N,K,T),2,t,dp,N,K,T),1,t,dp,N,K,T))

    #Third term: -<i| U* P U P U* P U P |i>
    C=-sum(conj(p_state)*Ut(dp*Ut(dp*Ut(dp*Ut(dp*p_state,2,t,dp,N,K,T),1,t,dp,N,K,T),2,t,dp,N,K,T),1,t,dp,N,K,T))

    #Fourth term: -<i| P U* P U P U* P U |i>
    D=-sum(conj(p_state)*dp*Ut(dp*Ut(dp*Ut(dp*Ut(p_state,2,t,dp,N,K,T),1,t,dp,N,K,T),2,t,dp,N,K,T),1,t,dp,N,K,T))

    otoc=(A+B+C+D)
    return(sqrt((otoc.real)**2+(otoc.imag)**2))#otoc.real)

#####################################################################################################################

def Ut(p_state,sign,Nkicks,dp,N,K,T=2**-7):
    #xvec definito dentro dp definito fuori 
    k=K/T
  
    # discretization x
    xmax=2*pi
    dx=xmax/N
    xvec=arange(0,xmax,dx) #0=2pi
    x_state=zeros([N],dtype = 'complex_')

    # hamiltonian and temporal evolution
    H0=(dp)**2/2 # mass=1
    V=k*cos(xvec) 
    if sign==2:
        U0=exp(-1j*H0*T)
        Uv=exp(-1j*V) 
    else:
        U0=exp(+1j*H0*T)
        Uv=exp(+1j*V)   

    temp_p=zeros([N],dtype = 'complex_') # Temporary states
    temp_x=zeros([N],dtype = 'complex_')

    #Alterno calci a propagazione libera

    if sign==1:
        x_state=(sqrt(N)*fft.ifft(p_state))
        

        #così se voglio un kick solo lo fa
        for j in range(1,Nkicks+1):
            #Kick
            temp_x=Uv*x_state
            #Fourier
            temp_p=(1/sqrt(N)*fft.fft(temp_x)) #La normalizzazione ad N è necessaria perchè altrimenti la trasformata non converge
            #Free evolution
            p_state=U0*temp_p
            #AntiFourier
            x_state=(sqrt(N)*fft.ifft(p_state)) 
    if sign==2:
        for j in range(1,Nkicks+1):
            #Free evolution
            temp_p=U0*p_state
            #AntiFourier
            temp_x=sqrt(N)*fft.ifft(temp_p)
            
            #Kick
            x_state=Uv*temp_x
            #Fourier
            p_state=(1/sqrt(N)*fft.fft(x_state))
           

    return(p_state)

#######################################################################################

def mean_otoc(N,T,K,Nkicks,trials,meanOTOC=True):
    '''function to set the useful discretizations and compute OTOC and its mean (set meanOTOC to True if you want to compute it)'''
    trials=int(trials)
    # discretization p
    dp=zeros([N])
    for i in range(0,N):
        dp[i]=i-(N/2)

    # initial conditions
    p_state=zeros(N,dtype = 'complex_')
    sigma=4
    # OTOC computation
    C=zeros([Nkicks,len(K),trials])
    #computes slightly different inital states in order to compute the mean OTOC 
    for ii in range(trials):
        #random initial states
        p0 = random.uniform(-pi,pi,N)
        p0=0
        p_state=(1/sqrt(2*pi*sigma**2*N**2))*exp(-((p0-dp)**2)/(2*sigma**2)) # gaussian
        
        for jj in range(len(K)):
            C[:,jj,ii]=Parallel(n_jobs=1)(delayed(OTOC)(p_state,aa,dp,N,K[jj],T) for aa in range(1,Nkicks+1))
        print('Trial #',ii)

    mean_C = abs(mean(C,2))
    if meanOTOC==True:
        out=mean_C
    else:
        out=C
    return(out)

################################################################
################################################################
################################################################


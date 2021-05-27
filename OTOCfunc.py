from numpy import *
import math
from joblib import Parallel, delayed

def save(input_info,name):
    '''Functio to save the input in a file .txt'''
    with open(name, 'w') as out:
        for line in input_info:
            savetxt(out, line)

#####################################################################################################################
# Here there is a set of function to compute the OTOC for the quantum kiched rotator

def OTOC(p_state,t,dp,N,K,T):
    '''Function to compute the OTOC on the state p_state at the time t with P(0)=dp.
    input: 
    p_state: State in momentum representation
    t: discretized time as number of kicks
    dp: initial momentum at time t=0
    N: number of eigenstate
    K: kick strength
    T: period
    output: computed OTOC'''


    #First term: <i| P U* P^2 U P |i>
    A=sum(conj(p_state)*dp*Ut(dp**2*Ut(dp*p_state,2,t,dp,N,K,T),1,t,dp,N,K,T))

    #Second term: <i| U* P U P^2 U* P U |i>
    B=sum(conj(p_state)*Ut(dp*Ut(dp**2*Ut(dp*Ut(p_state,2,t,dp,N,K,T),1,t,dp,N,K,T),2,t,dp,N,K,T),1,t,dp,N,K,T))

    #Third term: -<i| U* P U P U* P U P |i>
    C=-sum(conj(p_state)*Ut(dp*Ut(dp*Ut(dp*Ut(dp*p_state,2,t,dp,N,K,T),1,t,dp,N,K,T),2,t,dp,N,K,T),1,t,dp,N,K,T))

    #Fourth term: -<i| P U* P U P U* P U |i>
    D=-sum(conj(p_state)*dp*Ut(dp*Ut(dp*Ut(dp*Ut(p_state,2,t,dp,N,K,T),1,t,dp,N,K,T),2,t,dp,N,K,T),1,t,dp,N,K,T))

    otoc=(A+B+C+D)
    return(sqrt((otoc.real)**2+(otoc.imag)**2))

#####################################################################################################################

def Ut(p_state,sign,Nkicks,dp,N,K,T=2**-7):
    ''' Function to compute the Floquet operator implementing the time evolution operator, U, and its 
    Hermitian coniugate, U dagger, and using the Fourier transform and antitransform.
    The computation are implemented considering the mass m=1
    input: 
    p_state: State in momentum representation
    sign: if it is equal to 1 we compute U, if it is equal to 2 we compute U dagger
    Nkicks: number of kicks
    dp: initial momentum at time t=0
    N: number of eigenstate
    K: kick strength
    T: period, fixed to 2**-7 by default
    '''
    #useful quantities:
    k=K/T
    # position x discretization 
    xmax=2*pi
    dx=xmax/N # step size for the position space
    xvec=arange(0,xmax,dx) #0=2pi
    x_state=zeros([N],dtype = 'complex_') # state in position representation

    # hamiltonian and temporal evolution
    H0=(dp)**2/2
    V=k*cos(xvec) 
    if sign==2: # U
        U0=exp(-1j*H0*T)
        Uv=exp(-1j*V) 
    else:       # U dagger
        U0=exp(+1j*H0*T)
        Uv=exp(+1j*V)   

    temp_p=zeros([N],dtype = 'complex_') # Temporary states
    temp_x=zeros([N],dtype = 'complex_')

    if sign==1: # U dagger
        x_state=(sqrt(N)*fft.ifft(p_state))
        
        for j in range(1,Nkicks+1):
            #Kick
            temp_x=Uv*x_state
            #Fourier
            temp_p=(1/sqrt(N)*fft.fft(temp_x)) 
            #Free evolution
            p_state=U0*temp_p
            #AntiFourier
            x_state=(sqrt(N)*fft.ifft(p_state)) 
    if sign==2: #U
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
    '''Function to set the useful discretizations and to compute OTOC and its mean for different values of K.
    Set meanOTOC=True if you want to compute the mean of the OTOC
    input: 
    N: number of eigenstate
    T: period, fixed to 2**-7 by default
    K: kick strength
    Nkicks: number of kicks
    trials: number of trials to average
    meanOTOC: if it is True the function return the mean value of the OTOC, if it is False 
    function return the computed OTOC for each trial
    output:
    out: it is the the OTOC for each trial or the mean value of the OTOC, depending on the variable "meanOTOC"   
    '''
    # discretization of the momentum p
    N=2*N
    dp=zeros([N])
    for i in range(0,N):
        dp[i]=i-(N/2)

    # set initial conditions
    p_state=zeros(N, dtype = 'complex_')
    sigma=4 # as in the paper, it is fixed
    # OTOC computation
    C=zeros([Nkicks,len(K),trials])
    # Compute different inital state for each trial: different point in the phase space
    for ii in range(trials):
        #random initial states
        p0 = random.uniform(-pi,pi,N)
        p_state=(1/sqrt(2*pi*sigma**2*N**2))*exp(-((p0-dp)**2)/(2*sigma**2)) # gaussian as in the paper
        
        # Comupute the OTOC for different values of K, optimization of the code parallelizing over different number of kicks
        for jj in range(len(K)):
            C[:,jj,ii]=Parallel(n_jobs=1)(delayed(OTOC)(p_state,aa,dp,N,K[jj],T) for aa in range(1,Nkicks+1))
        print('Trial #',ii) # To observe the computed trials when the program is running

    mean_C = mean(C,2) # mean of the OTOC
    if meanOTOC==True:
        out=mean_C
    else:
        out=squeeze(C)
    return(out)


def mean_otoc_heff(N,T,K,Nkicks,trials,meanOTOC=True):
    '''Function to set the useful discretizations and to compute OTOC and its mean for 
    different values of h_eff(=T because here we consider h_bar=1).
    Set meanOTOC=True if you want to compute the mean of the OTOC
    input: 
    N: number of eigenstate
    T: period, fixed to 2**-7 by default
    K: kick strength
    Nkicks: number of kicks
    trials: number of trials to average
    meanOTOC: if it is True the function return the mean value of the OTOC, if it is False 
    function return the computed OTOC for each trial
    output:
    out: it is the the OTOC for each trial or the mean value of the OTOC, depending on the variable "meanOTOC"   
    '''
    # discretization p
    dp=zeros([N])
    for i in range(0,N):
        dp[i]=i-(N/2)

    # initial conditions
    p_state=zeros(N,dtype = 'complex_')
    sigma=4
    # OTOC computation
    C=zeros([Nkicks,len(T),trials])
    #computes slightly different inital states in order to compute the mean OTOC 
    for ii in range(trials):
        #random initial states
        p0 = random.uniform(-pi,pi,N)
        p_state=(1/sqrt(2*pi*sigma**2*N**2))*exp(-((p0-dp)**2)/(2*sigma**2)) # gaussian
        
        for jj in range(len(T)):
            C[:,jj,ii]=Parallel(n_jobs=1)(delayed(OTOC)(p_state,aa,dp,N,K,T[jj]) for aa in range(1,Nkicks+1))
        print('Trial #',ii)

    mean_C = abs(mean(C,2))
    if meanOTOC==True:
        out=mean_C
    else:
        out=C
    return(out)



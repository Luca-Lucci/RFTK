#!/usr/bin/env python
# coding: utf-8

# # uW

# ** DO NOT EDIT THIS .py FILE. EDIT THE JUPYTER NOTEBOOK **

# library with routines ported from T. Reveyrand free codes found in `http://www.microwave.fr/uW.html` \
# local copies of the original codes are stored in `C:\home\Zibaldone\TRL_Deembedding_Googling\TRL_Deembedding_Googling\macros`
# 
# To actually create the module, export as an executable python script

# Changelog

# - v 0.2: fix S2Y Y22 bug and also minor typos in comments.
# - v 0.1: porting of first macros

# In[ ]:


# this line is useless, will be added by export filter
# -*- coding: utf-8 -*-

# import re   as re 
import numpy as np  # so far is the only dependable module to be included
# import pickle

################################################################################


# ## Utility 

# Utility functions are found for the moment at the header of each library, and are copied over with inefficient copy-paste. All routines sholud be found in common.py lib.

# In[ ]:


def isfloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False

def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n


# Basic notes on porting scilab math into PYthon. 
# - `./` and `.*` are element-wise division and multiplication of two vectors with same size/length.
# - some constant are predefined: `%i`, `%pi` and `%e` are imaginary unit, pi-number and neper

# The key function for matrix conversion is a sort of wrapper that is abstract for a list of checks that would have been copied over to each and every conversion function. A simple function will be created to convert a 2x2 matrix say from impedance `Z` to admirttances `Y`. The NextConv function will handle the cases of a nested matrix (a frequency sweep) a doble nested matrix (a bias sweep and the a frequency sweep). We avoid using objects for the moment and we ignore the frequency sweeep. the simple function that is actually performing the conversion is passed as an argument to the NestConv function.

# In[ ]:


def NestConv(M, conv):
    # conversion of M to a target format specified by chosing rigth conv function
    # supposed to be internal function to the module
    shape=np.shape(M)
    if (shape[-1] != 2) or (shape[-2] != 2):
        raise ValueError("Function is defined for 2x2 matrix exclusively")
    if len(shape)>4:
        raise ValueError("Matrix is nested more than 2 times (bias and frequency. aborting")
    if len( np.shape(M) )==2: # intercept a matrix conversion without f sweep, i,e Z[i,j]           
        N=conv(M)
        return N
    N=np.empty(np.shape(M),dtype=np.complex128) # simple frequency sweep
    if len( np.shape(M) )==3: # intercept a matrix conversion with f sweep, i,e Z[f,i,j]   
        for f in range(len(M)): # each frequency sweep
            N[f]=conv(M[f])
        return N      
    if len( np.shape(M) )==4: # bias and frequency sweep 
        bn,fn,x,y=np.shape(M)  # intercept a matrix conversion with bias and f sweep, i,e Z[b,f,i,j] 
        for b in range(bn): # bias sweep
            for f in range(fn): # frequency sweep
                N[b,f]=conv(M[b,f])
        return N
    if len( np.shape(M) )==5: # die or temperature, bias and frequency sweep 
        dn, bn,fn,x,y=np.shape(M)  # intercept a matrix conversion with bias and f sweep, i,e Z[t,b,f,i,j]
        for d in range(dn):
            for b in range(bn): # bias sweep
                for f in range(fn): # frequency sweep
                    N[d,b,f]=conv(M[b,f])
        return N    


# Elementary matrix transformations. These are fast routines to transform matrices in 50-Ohms system
# - ZtoY and YtoZ are the same function, basically one is the inverse of the other
# - ZtoS
# - TtoS
# - StoZ
# - StoY
# - StoT
# - StpH

# In[ ]:


def matrixZ2Y(Z):
    # no error checks. suppose to have a 2x2 matrix to convert and that's it
    Y=np.empty((2,2), dtype=np.complex128)
    z11=Z[0,0]
    z12=Z[0,1]
    z21=Z[1,0]
    z22=Z[1,1]
    Y[0,0]=z22/(z11*z22-z12*z21)
    Y[0,1]=-z12/(z11*z22-z12*z21)
    Y[1,0]=-z21/(z11*z22-z12*z21)
    Y[1,1]=z11/(z11*z22-z12*z21)
    return Y


# In[ ]:


def matrixZ2S(Z):
    # no error checks. suppose to have a 2x2 matrix to convert and that's it
    S=np.empty((2,2), dtype=np.complex128)
    z11=Z[0,0]
    z12=Z[0,1]
    z21=Z[1,0]
    z22=Z[1,1]
    z0=50
    c=(z11+z0)*(z22+z0)-(z12*z21)
    S[1,1]=((z11+z0)*(z22-z0)-(z12*z21))/c
    S[1,0]=(2*z0*z21)/c
    S[0,1]=(2*z0*z12)/c
    S[0,0]=((z11-z0)*(z22+z0)-(z12*z21))/c
    return S


# In[ ]:


def matrixT2S(T):
    # no error checks. suppose to have a 2x2 matrix to convert and that's it
    S=np.empty((2,2), dtype=np.complex128)
    t11=T[0,0]
    t12=T[0,1]
    t21=T[1,0]
    t22=T[1,1]
    S[0,0]=t21/t11
    S[0,1]=(t11*t22-t12*t21)/t11
    S[1,0]=1/t11
    S[1,1]=-t12/t11
    return S


# In[ ]:


def matrixS2Z(S):
    # no error checks. suppose to have a 2x2 matrix to convert and that's it
    Z=np.empty((2,2), dtype=np.complex128)
    s11=S[0,0]
    s12=S[0,1]
    s21=S[1,0]
    s22=S[1,1]
    z0=50
    c=z0/((1-s11)*(1-s22)-s12*s21)
    Z[0,0]=c*((1+s11)*(1-s22)+s12*s21)
    Z[0,1]=c*(2*s12)
    Z[1,0]=c*(2*s21)
    Z[1,1]=c*((1-s11)*(1+s22)+s12*s21)
    return Z


# In[ ]:


def matrixS2Y(S):
    # no error checks. suppose to have a 2x2 matrix to convert and that's it
    Y=np.empty((2,2), dtype=np.complex128)
    s11=S[0,0]
    s12=S[0,1]
    s21=S[1,0]
    s22=S[1,1]
    z0=50
    c=1/(z0*((1+s11)*(1+s22)-s12*s21))
    Y[0,0]=c*((1-s11)*(1+s22)+s12*s21)
    Y[0,1]=c*(-2*s12)
    Y[1,0]=c*(-2*s21)
    Y[1,1]=c*((1+s11)*(1-s22)+s12*s21)   # bug in V0.1 corrected
    return Y


# In[ ]:


def matrixS2T(S):
    # no error checks. suppose to have a 2x2 matrix to convert and that's it
    T=np.empty((2,2), dtype=np.complex128)
    s11=S[0,0]
    s12=S[0,1]
    s21=S[1,0]
    s22=S[1,1]
    T[0,0]=1/s21
    T[0,1]=-s22/s21
    T[1,0]=s11/s21
    T[1,1]=(s12*s21-s11*s22)/s21
    return T


# In[ ]:


def matrixS2H(S):
    # no error checks. suppose to have a 2x2 matrix to convert and that's it
    H=np.empty((2,2), dtype=np.complex128)
    s11=S[0,0]
    s12=S[0,1]
    s21=S[1,0]
    s22=S[1,1]  
    d=(1-s11)*(1+s22)+s12*s21
    y0=1/50
    H[0,0]=(1/y0)*((1+s11)*(1+s22)-s12*s21)/d
    H[0,1]=(2*s12)/d
    H[1,0]=(-2*s21)/d
    H[1,1]=y0*((1-s11)*(1-s22)-s12*s21)/d
    return H


# In[ ]:


def matrixZ2ABCD(Z):
    # addded by inspecting ael in ADS code. not in original TR code!
    A=np.empty((2,2), dtype=complex)
    A[0,0] = 1.0*Z[0,0]/Z[1,0]
    A[0,1] = 1.0*(Z[0,0]*Z[1,1]-Z[1,0]*Z[0,1])/Z[1,0]
    A[1,0] = 1.0/Z[1,0]
    A[1,1] = 1.0*Z[1,1]/Z[1,0]
    return A


# In[15]:


def matrixABCD2Y(A):
    # addded by inspecting ael in ADS code. not in original TR code!
    Y=np.empty((2,2), dtype=complex)
    Y[0,0] = A[1,1]/A[0,1]
    Y[0,1] = -1.0*(A[0,0]*A[1,1]-A[1,0]*A[0,1])/A[0,1]
    Y[1,0] = -1.0/A[0,1]
    Y[1,1] =A[0,0]/A[0,1]
    return Y


# In[ ]:


def matrixS2ABCD(S):
    # addded by inspecting ael in ADS code. not in original TR code!
    # convert S to Z then Z to ABCD
    return matrixZ2ABCD(matrixS2Z(S))  


# In[1]:


def matrixY2ABCD(Y):
    # addded by inspecting ael in ADS code. not in original TR code!
    A=np.empty((2,2), dtype=complex)   
    A[0,0] = -1.0*Y[1,1]/Y[1,0]
    A[0,1] = -1.0/Y[1,0]
    A[1,0] = -1.0*(Y[0,0]*Y[1,1]-Y[1,0]*Y[0,1])/Y[1,0]
    A[1,1] = -1.0*Y[0,0]/Y[1,0]
    return A


# In[1]:


def matrixY2S(Y):
    # addded after, in ADS is different, need to specify the impedance
    # it is an inversion so matrixY2Z == matrixZ2Y
    return matrixZ2S(matrixZ2Y(Y))


# In[16]:


def Z2Y(Z):
    # conversion of Z to Y matrix of 50 Ohm data
    # valid for a a matrix, a matrix with a frequency sweep and a matrix with a bias and frequency sweep
    return NestConv(Z, matrixZ2Y)

def Y2Z(Y):
    # conversion of Y to Z matrix of 50 Ohm data
    # valid for a a matrix, a matrix with a frequency sweep and a matrix with a bias and frequency sweep
    return NestConv(Y, matrixZ2Y)

def Z2S(Z):
    # conversion of Z to S matrix of 50 Ohm data
    # valid for a a matrix, a matrix with a frequency sweep and a matrix with a bias and frequency sweep
    return NestConv(Z, matrixZ2S)

def T2S(T):
    # conversion of S to Z matrix of 50 Ohm data
    # valid for a a matrix, a matrix with a frequency sweep and a matrix with a bias and frequency sweep
    return NestConv(T, matrixT2S)

def S2Z(S):
    # conversion of S to Z matrix of 50 Ohm data
    # valid for a a matrix, a matrix with a frequency sweep and a matrix with a bias and frequency sweep
    return NestConv(S, matrixS2Z)

def S2Y(S):
    # conversion of S to Y matrix of 50 Ohm data
    # valid for a a matrix, a matrix with a frequency sweep and a matrix with a bias and frequency sweep
    return NestConv(S, matrixS2Y)

def S2T(S):
    # conversion of S to T matrix of 50 Ohm data
    # valid for a a matrix, a matrix with a frequency sweep and a matrix with a bias and frequency sweep
    return NestConv(S, matrixS2T)

def S2H(S):
    # conversion of S to H matrix of 50 Ohm data
    # valid for a a matrix, a matrix with a frequency sweep and a matrix with a bias and frequency sweep
    return NestConv(S, matrixS2H)

def Z2ABCD(Z):
    # conversion of S to H matrix of 50 Ohm data
    # valid for a a matrix, a matrix with a frequency sweep and a matrix with a bias and frequency sweep
    return NestConv(Z, matrixZ2ABCD)

def S2ABCD(S):
    # conversion of S to H matrix of 50 Ohm data
    # valid for a a matrix, a matrix with a frequency sweep and a matrix with a bias and frequency sweep
    return NestConv(S, matrixS2ABCD)

def ABCD2Y(A):
    # conversion of S to H matrix of 50 Ohm data
    # valid for a a matrix, a matrix with a frequency sweep and a matrix with a bias and frequency sweep
    return NestConv(A, matrixABCD2Y)

def Y2ABCD(Y):
    # conversion of S to H matrix of 50 Ohm data
    # valid for a a matrix, a matrix with a frequency sweep and a matrix with a bias and frequency sweep
    return NestConv(Y, matrixY2ABCD)


def Y2S(Y):
    # conversion of S to H matrix of 50 Ohm data
    # valid for a a matrix, a matrix with a frequency sweep and a matrix with a bias and frequency sweep
    return NestConv(Y, matrixY2S)


# simply convert offset in mm to delay in sec

# In[ ]:


def offset2delay(length_mm, approx=False):
    """
    From TReveyrand code uW_cal_offset2delay(l)
    length of offset is in [mm],
    the routine from Tibault uses rather crude approximantion, so I added a switch
    """
    if approx:
        c=3.0e8 # speed of light in m/s in vacuum
    else:
        c=299,792,458.0
    length=length_mm/1000.0 # in meters
    delay=length/c
    return delay


# Generates the $S_{11}$ or $\Gamma$ for an ideal SOL calibration elements from its parameters

# In[ ]:


# // #######################################
# // # Génération du S11 d'un standard
# // #######################################

# function S=uW_cal_standard(f,std)
#     Zo=50;
#     select typeof(std),
#     case "Open standard",
#         Ceff=std.C0+std.C1*f+std.C2*f.^2+std.C3*f^3;
#         S11=(1-%i*100*%pi*f.*Ceff)./(1+%i*100*%pi*f.*Ceff);
#     case "Short standard",
#         Leff=std.L0+std.L1*f+std.L2*f.^2+std.L3*f^3;
#         S11=(%i*2*%pi*f.*Leff-50)./(%i*2*%pi*f.*Leff+50);
#     case "Match standard",
#         Z=std.R0+%i*2*%pi*std.L0;
#         S11=(Z-50)./(Z+50);
#     end;
    
#     Loss=exp( (-std.delay/Zo)*std.loss*sqrt(f./(10^9)) );
#     Delay=exp(-%i*4*%pi*std.delay*f);
#     S11=Loss.*Delay.*S11;
    
#     S=tlist(['S parameters';'frequency';'S11'],f,S11);
# endfunction


# In[ ]:


def GammaShort(L0):
    # TBD
    return L0


# My take of the scilab code for a simple TRL

# There are two version of this code. The code saved in the unwrap routines

# In[ ]:


def unwrap(vec, threshold=200,delta=360):
    # a=[1,2,3,4]; a[1:]->[2, 3, 4]; a[:-1]->[1, 2, 3]
    c=np.zeros(len(vec), dtype=np.float64)
    unwrap=np.zeros(len(vec), dtype=np.float64)
    for i in range(len(c)):
        if i==0: continue
        c[i]=deg(vec[i])-deg(vec[i-1])
        if c[i]<-threshold: a[i]=delta
        elif c[i]>threshold: a[i]=-delta
        else: c[i]=0
        c[i]+=c[i-1]  # cumulative sum
        unwrap[i]=deg[i]+c[i]
    return unwrap


# In[ ]:


def S2P_deembedding(Stotal,Sin,Sout):
    # I skip this part I think it just manages the case where the freequency sweep 
    # is different in Stotoal an Sin, Sout
#     for k=1:4,
#         Sij=Sin(2+k);
#         M_Sij=abs(Sij);
#         P_Sij=uW_unwarp(Sij);
        
#         M2_Sij=interp(frequence, f_Sij,M_Sij, splin(f_Sij,M_Sij,"monotone")); 
#         P2_Sij=interp(frequence, f_Sij,P_Sij, splin(f_Sij,P_Sij,"monotone")); 
            
#         Sin(2+k)=M2_Sij.*exp(%i*%pi*P2_Sij/180);
#     end;
    
#     // Sout
#         f_Sij=Sout.frequency;
#     for k=1:4,
#         Sij=Sout(2+k);
#         M_Sij=abs(Sij);
#         P_Sij=uW_unwarp(Sij);
        
#         M2_Sij=interp(frequence, f_Sij,M_Sij, splin(f_Sij,M_Sij,"monotone")); 
#         P2_Sij=interp(frequence, f_Sij,P_Sij, splin(f_Sij,P_Sij,"monotone")); 
            
#         Sout(2+k)=M2_Sij.*exp(%i*%pi*P2_Sij/180);
#     end;
    Sm=np.empty(np.shape(Stotal), dtype=np.complex128)
    for f in len(Stotal):
        #// Calcul de TA^(-1) et TB^(-1) => Sdut
        S11A=Sin[f,0,0]
        S12A=Sin[f,0,1]
        S21A=Sin[f,1,0]
        S22A=Sin[f,1,1]

        S11B=Sout[f,0,0]
        S12B=Sout[f,0,1]
        S21B=Sout[f,1,0]
        S22B=Sout[f,1,1]
        
        TAI=(1/S12A)*np.array([[(-S11A*S22A+S12A*S21A),S22A],[-S11A,1]])
        TBI=(1/S12B)*np.array([[(-S11B*S22B+S12B*S21B),S22B],[-S11B,1]])
        Ttotal=(1/Stotal[f,1,0])*np.array([[1,-Stotal[f,1,1]],[Stotal[f,0,0],(Stotal[f,0,1]*Stotal[f,1,0]-Stotal[f,0,0]*Stotal[f,1,1])]])
        Tdut=np.dot(TAI,np.dot(Ttotal,TBI))
        Sdut=(1/Tdut[0,0])*np.array([[Tdut[1,0],Tdut[0,0]*Tdut[1,1]-Tdut[1,0]*Tdut[0,1]],[1,-Tdut[0,1]]])
        
        Sm[f]=Sdut
    return Sdut


# In[ ]:


###### simple TRL calibration
def TibaultR_TRL(Sthru, Sline, Sreflect, Reflect_type='Open'):
    Tthru=S2T(Sthru)
    Tline=S2T(Tline)
    T1=Tthru
    T2=Tline
    K   =np.zeros(len(Sthru), dtype=np.complex128)
    Tin =np.zeros(np.shape(Sthru), dtype=np.complex128)
    Tout=np.zeros(np.shape(Sthru), dtype=np.complex128)
    # works at a single frequency f 
    for f in range(len(Thru)):
        T1=Tthru[f]
        T2=Tthru[f]
        M=np.matmul(np.linalg.inv(T1),T2)
        N=np.matmul(T2,np.linalg.inv(T1))
        # equation TRL M
        delta=(M[0,0]-M[1,1])^2-4*M[1,0]*(-M[0,1])
        X1=((M[1,1]-M[0,0])-sqrt(delta))/(2*M[1,0])
        X2=((M[1,1]-M[0,0])+sqrt(delta))/(2*M[1,0])
        Sol=np.array([X1,X2])
        # these two lines are blowing my mind ...
        # C1=Sol(find(abs(Sol)==max(abs(Sol))))  #;  // OUTPUT - T22/T21
        # C2=Sol(find(abs(Sol)==min(abs(Sol))))  # ;  // OUTPUT - T12/T11    #
        C1=max(Sol)
        C2=min(Sol)
        
        # equation TRL N
        delta=(N[0,0]-N[1,1])^2-4*N[1,0]*(-N[0,1])
        X1=((N[1,1]-N[0,0])-sqrt(delta))/(2*N[1,0])
        X2=((N[1,1]-N[0,0])+sqrt(delta))/(2*N[1,0])
        Sol=np.array([X1,X2])
        C3=max(Sol)
        C4=max(Sol)
        # Equation THRU FORWARD
        C5=(1+C4*Sthru[f,0,0])/Sthru[f,1,0] # T11/A
        # Equation THRU REVERSE    
        C6=Sthru[f,0,1]/(Sthru[f,1,1]+C1) # T21/D
        # Equation Reflect
        _X=(Sreflect[f,1,1]+C2)/(Sreflect[f,0,0]+C1)
        _Y=(1+C3*Sreflect[f,0,0])/(1+C4*Sreflect[f,0,0])
        sol=sqrt((C5*_X)/(_Y*C6*C3)) #  C

        A=1
        B=C4
        C=sol
        D=C3*C
        GAMMA_STD=(C+D*Sreflect[f,0,0])/(A+B*Sreflect[f,0,0])
        if ((np.real(GAMMA_STD)>0) and  (REFLECT_STD=="SHORT")) or ( (np.real(GAMMA_STD)<0) and (REFLECT_STD=="OPEN") ):  ##    // If OPEN, chose the other solution
             sol=-sol 
        K[f]=np.sqrt(1/(A*D-B*C))
        Tin[f] =np.linalg.inv(np.array[[A,B],[C,D]])
        Tout[f]=np.array([[C5,C2*C5],[C6*D,C1*C6*D]])
        uwp=unwrap(K, threshold=90,delta=180)
####    mind blowing code:    
#        cc=pinv([S_thru.frequency/10^9,uwp*0+1])*uwp;
#       ### // if  (cc(2)>180) then, uwp=uwp-360;cc(2)=cc(2)-360;end;
#       ### // if  (cc(2)<-180) then, uwp=uwp+360;cc(2)=cc(2)+360;end;
#        uwp=uwp-cc(2)+modulo(cc(2),180);
        K=np.abs(K)*exp(1j*np.pi*uwp/180)
    Pin =np.zeros(np.shape(Sthru), dtype=np.complex128)
    Pout=np.zeros(np.shape(Sthru), dtype=np.complex128)
    # works at a single frequency f 
    for f in range(len(Thru)):
        Pin[f,0,0]=T[f,0,0]/K[f]
        Pin[f,0,1]=T[f,0,1]/K[f]
        Pin[f,1,0]=T[f,1,0]/K[f]
        Pin[f,1,1]=T[f,1,1]/K[f]
        
        Pout[f,0,0]=T[f,0,0]*K[f]
        Pout[f,0,1]=T[f,0,1]*K[f]
        Pout[f,1,0]=T[f,1,0]*K[f]
        Pout[f,1,1]=T[f,1,1]*K[f]
        
    Sin=T2S(Pin)
    Sout=T2S(Pout)
    
    Sr=S2P_deembedding(Sreflect,Sin,Sout)
    return Sin, Sout, Sr


# In[ ]:





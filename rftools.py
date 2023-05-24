# -*- coding: utf-8 -*-
# v1.2 28/3/2023 fix np.complex into np.complex128 since deprecation
# v1.1 28/3/2023 improve touchstone export
# 31/3/2021 - Update of findHz. Backward compatible

#import skrf as rf
import re   as re 
import numpy as np
import _pickle

################################################################################
### UTILITY
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

################################################################################
### Spectre interface
def spectre2S(freq, s11, s12, s21, s22):
    S=np.empty([len(freq), 2, 2],dtype=np.complex128)
    f=np.empty(len(freq),dtype=np.float_)
    for i,ff in enumerate(freq):
        f[i]=freq[i]
        S[i,0,0]=s11[i]
        S[i,0,1]=s12[i]
        S[i,1,0]=s21[i]
        S[i,1,1]=s22[i]
    return f, S

################################################################################
## simple cases : 2 ports, 50Ohms

def StoT_2P(S):
    T=np.empty(np.shape(S),dtype=np.complex128)
    for i in range(len(S)):
        s11=S[i,0,0]
        s12=S[i,0,1]
        s21=S[i,1,0]
        s22=S[i,1,1]
        T[i,0,0]=1/s21
        T[i,0,1]=-s22/s21
        T[i,1,0]=s11/s21
        T[i,1,1]=(s12*s21-s11*s22)/s21
    return T

def TtoS_2P(T):
    S=np.empty(np.shape(T),dtype=np.complex128)
    for i in range(len(T)):
        t11=T[i,0,0]
        t12=T[i,0,1]
        t21=T[i,1,0]
        t22=T[i,1,1]
        S[i,0,0]=t21/t11
        S[i,0,1]=(t11*t22-t12*t21)/t11
        S[i,1,0]=1/t11
        S[i,1,1]=-t12/t11
    return S

def YtoA_2P(Y):
    A=np.empty(np.shape(Y),dtype=np.complex128)
    for i in range(len(Y)):
        y11=Y[i,0,0]
        y12=Y[i,0,1]
        y21=Y[i,1,0]
        y22=Y[i,1,1]
        A[i,0,0]=-y22/y21
        A[i,0,1]=-1/y21
        A[i,1,0]=(y12*y21-y11*y22)/y21
        A[i,1,1]=-y11/y21
    return A

def HtoA_2P(H):
    A=np.empty(np.shape(H),dtype=np.complex128)
    for i in range(len(H)):
        h11=H[i,0,0]
        h12=H[i,0,1]
        h21=H[i,1,0]
        h22=H[i,1,1]
        A[i,0,0]=(h12*h21-h11*h22)/h21
        A[i,0,1]=-h11/h21
        A[i,1,0]=-h22/h21
        A[i,1,1]=-1/h21
    return A

def ZtoA_2P(Z):
    A=np.empty(np.shape(Z),dtype=np.complex128)
    for i in range(len(Z)):
        z11=Z[i,0,0]
        z12=Z[i,0,1]
        z21=Z[i,1,0]
        z22=Z[i,1,1]
        A[i,0,0]=z11/z21
        A[i,0,1]=(z11*z22-z12*z21)/z21
        A[i,1,0]=1/z21
        A[i,1,1]=z22/z21
    return A

def AtoH_2P(A):
    H=np.empty(np.shape(A),dtype=np.complex128)
    for i in range(len(A)):
        a=A[i,0,0]
        b=A[i,0,1]
        c=A[i,1,0]
        d=A[i,1,1]
        H[i,0,0]=b/d
        H[i,0,1]=(a*d-b*c)/d
        H[i,1,0]=-1/d
        H[i,1,1]=c/d
    return H

def AtoY_2P(A):
    Y=np.empty(np.shape(A),dtype=np.complex128)
    for i in range(len(A)):
        a=A[i,0,0]
        b=A[i,0,1]
        c=A[i,1,0]
        d=A[i,1,1]
        Y[i,0,0]=d/b
        Y[i,0,1]=(b*c-a*d)/b
        Y[i,1,0]=-1/b
        Y[i,1,1]=a/b
    return Y

def YtoZ_2P(Y):
    Z=np.empty(np.shape(Y),dtype=np.complex128)
    for i in range(len(Y)):
        y11=Y[i,0,0]
        y12=Y[i,0,1]
        y21=Y[i,1,0]
        y22=Y[i,1,1]
        det=(y11*y22-y12*y21)
#        if det==0:
#               det=1e-15
        Z[i,0,0]=y22/det
        Z[i,0,1]=-y12/det
        Z[i,1,0]=-y21/det
        Z[i,1,1]=y11/det
    return Z

def AtoZ_2P(A):
    Z=np.empty(np.shape(A),dtype=np.complex128)
    for i in range(len(A)):
        a=A[i,0,0]
        b=A[i,0,1]
        c=A[i,1,0]
        d=A[i,1,1]
        Z[i,0,0]=a/c
        Z[i,0,1]=(a*d-b*c)/c
        Z[i,1,0]=1/c
        Z[i,1,1]=d/c
    return Z

def HtoZ_2P(H):
    Z=np.empty(np.shape(H),dtype=np.complex128)
    for i in range(len(H)):
        h11=H[i,0,0]
        h12=H[i,0,1]
        h21=H[i,1,0]
        h22=H[i,1,1]
        Z[i,0,0]=(h11*h22-h12*h21)/h22
        Z[i,0,1]=h12/h22
        Z[i,1,0]=-h21/h22
        Z[i,1,1]=1/h22
    return Z

def ZtoY_2P(Z):
    Y=np.empty(np.shape(Z),dtype=np.complex128)
    for i in range(len(Z)):
        if len(Z)==2:
            
            z11=Z[0,0]
            z12=Z[0,1]
            z21=Z[1,0]
            z22=Z[1,1]
            Y[0,0]=z22/(z11*z22-z12*z21)
            Y[0,1]=-z12/(z11*z22-z12*z21)
            Y[1,0]=-z21/(z11*z22-z12*z21)
            Y[1,1]=z11/(z11*z22-z12*z21)
            
        else : 
            
            z11=Z[i,0,0]
            z12=Z[i,0,1]
            z21=Z[i,1,0]
            z22=Z[i,1,1]
            Y[i,0,0]=z22/(z11*z22-z12*z21)
            Y[i,0,1]=-z12/(z11*z22-z12*z21)
            Y[i,1,0]=-z21/(z11*z22-z12*z21)
            Y[i,1,1]=z11/(z11*z22-z12*z21)
            
    return Y

def ZtoH_2P(Z):
    H=np.empty(np.shape(Z),dtype=np.complex128)
    for i in range(len(Z)):
        z11=Z[i,0,0]
        z12=Z[i,0,1]
        z21=Z[i,1,0]
        z22=Z[i,1,1]
        H[i,0,0]=(z11*z22-z12*z21)/z22
        H[i,0,1]=z12/z22
        H[i,1,0]=-z21/z22
        H[i,1,1]=1/z22
    return H

def YtoH_2P(Y):
    H=np.empty(np.shape(Y),dtype=np.complex128)
    for i in range(len(Y)):
        y11=Y[i,0,0]
        y12=Y[i,0,1]
        y21=Y[i,1,0]
        y22=Y[i,1,1]
        H[i,0,0]=1/y11
        H[i,0,1]=-y12/y11
        H[i,1,0]=y21/y11
        H[i,1,1]=(y11*y22-y12*y21)/y11
    return H


def HtoY_2P(H):
    Y=np.empty(np.shape(H),dtype=np.complex128)
    for i in range(len(H)):
        h11=H[i,0,0]
        h12=H[i,0,1]
        h21=H[i,1,0]
        h22=H[i,1,1]
        Y[i,0,0]=1/h11
        Y[i,0,1]=-h12/h11
        Y[i,1,0]=h21/h11
        Y[i,1,1]=(h11*h22-h12*h21)/h11
    return Y


################################################################################
### not so simple cases, 2 ports, generic Z01 Z02

def YtoS(Y, z0=[50.0+0.0j, 50.0+0.0j]):
   # dovrebbe essere coerente con 
   # la funzione s2y del rf-scikit
    S=np.empty(np.shape(Y),dtype=np.complex128)
    
    z1=z0[0]
    z1c=np.conj(z1)
    r1=np.real(z1)
    z2=z0[1]
    z2c=np.conj(z2)
    r2=np.real(z2)

    for i in range(len(Y)):
        y11=Y[i,0,0]
        y12=Y[i,0,1]
        y21=Y[i,1,0]
        y22=Y[i,1,1]
        den=(1+y11*z1)*(1+y22*z2)-y12*y21*z1*z2
        s11=(1-y11*z1c)*(1+y22*z2)+y12*y21*z1c*z2
        s12=-2.0*y12*np.sqrt(r1*r2)
        s21=-2.0*y21*np.sqrt(r1*r2)
        s22=(1+y11*z1)*(1-y22*z2c)+y12*y21*z1*z2c

        s11=s11/den
        s12=s12/den
        s21=s21/den
        s22=s22/den

        S[i,0,0]=s11
        S[i,0,1]=s12
        S[i,1,0]=s21
        S[i,1,1]=s22
    return S

def StoY(S, z0=[50.0+0.0j, 50.0+0.0j]):
   # dovrebbe essere coerente con 
   # la funzione s2y del rf-scikit
    Y=np.empty(np.shape(S),dtype=np.complex128)
    
    z1=z0[0]
    z1c=np.conj(z1)
    r1=np.real(z1)
    z2=z0[1]
    z2c=np.conj(z2)
    r2=np.real(z2)

    for i in range(len(S)):
        s11=S[i,0,0]
        s12=S[i,0,1]
        s21=S[i,1,0]
        s22=S[i,1,1]
        den=(z1c+s11*z1)*(z2c+s22*z2)-s12*s21*z1*z2
        y11=(1-s11)*(z2c+s22*z2)+s12*s21*z2
        y12=-2*s12*np.sqrt(r1*r2)
        y21=-2*s21*np.sqrt(r1*r2)
        y22=(z1c+s11*z1)*(1-s22)+s12*s21*z1

        y11=y11/den
        y12=y12/den
        y21=y21/den
        y22=y22/den

        Y[i,0,0]=y11
        Y[i,0,1]=y12
        Y[i,1,0]=y21
        Y[i,1,1]=y22
    return Y

def StoH(S):
    H=np.empty(np.shape(S),dtype=np.complex128)
    Z0=50
    for i in range(len(S)):
        s11=S[i,0,0]
        s12=S[i,0,1]
        s21=S[i,1,0]
        s22=S[i,1,1]
        den=(1-s11)*(Z0+s22*Z0)+s12*s21*Z0
        h11=((Z0+s11*Z0)*(Z0+s22*Z0)-s12*s21*Z0*Z0)/den
        h12=2*s12*Z0/den
        h21=-2*s21*Z0/den
        h22=((1-s11)*(1-s22)-s12*s21)/den
        H[i,0,0]=h11
        H[i,0,1]=h12
        H[i,1,0]=h21
        H[i,1,1]=h22
    return H

def StoABCD(S, Z0=[complex(50,0), complex(50,0)]):
    A=np.empty(np.shape(S),dtype=np.complex128)
    R01=np.real(Z0[0])
    Z01=Z0[0]
    Z01c=np.conj(Z01)
    R02=np.real(Z0[1])
    Z02=Z0[1]
    Z02c=np.conj(Z02)
    for i in range(len(S)):
        s11=S[i,0,0]
        s12=S[i,0,1]
        s21=S[i,1,0]
        s22=S[i,1,1]
        den=2*s21*np.sqrt(R01*R02)
        a11=( (Z01c+s11*Z01)*(1.0-s22)+s12*s21*Z01 )/den
        a12=( (Z01c+s11*Z01)*(Z02c+s22*Z02)-s21*s12*Z01*Z02 )/den
        a21=( (1.0-s11)*(1.0-s22)-s21*s12 )/den
        a22=( (1.0-s11)*(Z02c+s22*Z02)+s12*s21*Z02 )/den 
        A[i,0,0]=a11
        A[i,0,1]=a12
        A[i,1,0]=a21
        A[i,1,1]=a22
    return A

################################################################################
#### 2 port extractions

def H2Ptoft(H, fHz, ft_start=0, ft_stop=-1):
    #ftv =dict() # da tutto H21
    #ft_start=18 
    #ft_stop=35
    #ft=dict() # solo tra ft_start ft_stop
    h21=H[:,1,0]
    ftv=(fHz*np.abs(h21))
    ft =sum(ftv[ft_start:ft_stop])/1e9/len(ftv[ft_start:ft_stop]) # /1e9 per GHz
    return ft, ftv


################################################################################
## 4 ports
# dovrebbe essere coerente con 
# la funzione s2y del rf-scikit

def StoZ4P(S, z0=[complex(50.0,0.0), complex(50.0,0.0), complex(50.0,0.0), complex(50.0,0.0)]):
    Z=np.empty(np.shape(S),dtype=np.complex128)
    GG =np.diag(z0)
    FF =np.diag(1.0/(2.0*np.sqrt(np.real(np.array(z0)))))
    II =np.diag([1.0, 1.0, 1.0, 1.0])
    for i in range(len(S)):
            SS=S[i]
            Z[i]=np.dot( np.linalg.inv(FF), 
                   np.dot( np.linalg.inv(II-SS),
                     np.dot(II+SS, 
                       np.dot(GG, FF) ) ) )
    return Z

def StoY4P(S, z0=[complex(50.0,0.0), complex(50.0,0.0), complex(50.0,0.0), complex(50.0,0.0)]):
    Y=np.empty(np.shape(S),dtype=np.complex128)
    GG =np.diag(z0)
    FF =np.diag(1.0/(2.0*np.sqrt(np.real(np.array(z0)))))
    II =np.diag([1.0, 1.0, 1.0, 1.0])
    for i in range(len(S)):
            SS=S[i]
            Y[i]=np.dot( np.linalg.inv(FF), 
                   np.dot( np.linalg.inv(GG),
                     np.dot( np.linalg.inv(II+SS), 
                       np.dot(II-SS, FF) ) ) )
    return Y

def YtoS4P(Y, z0=[complex(50.0,0.0), complex(50.0,0.0), complex(50.0,0.0), complex(50.0,0.0)]):
    S=np.empty(np.shape(Y),dtype=np.complex128)
    GG =np.diag(z0)
    FF =np.diag(1.0/(2.0*np.sqrt(np.real(np.array(z0)))))
    II =np.diag([1.0, 1.0, 1.0, 1.0])
    for i in range(len(Y)):
            YY=Y[i]
            S[i]=np.dot( FF, 
                   np.dot( II - np.dot(GG,YY),
                     np.dot( np.linalg.inv( II +  np.dot(GG+YY)) , np.linalg.inv(FF) ) ) ) 
    return S

def ZtoS4P(Z, z0=[complex(50.0,0.0), complex(50.0,0.0), complex(50.0,0.0), complex(50.0,0.0)]):
    S=np.empty(np.shape(Z),dtype=np.complex128)
    GG =np.diag(z0)
    FF =np.diag(1.0/(2.0*np.sqrt(np.real(np.array(z0)))))
    II =np.diag([1.0, 1.0, 1.0, 1.0])
    for i in range(len(Z)):
            ZZ=Z[i]
            S[i]=np.dot( FF, 
                   np.dot( ZZ - GG,
                     np.dot( np.linalg.inv( ZZ +  GG) , np.linalg.inv(FF) ) ) ) 
    return S

def ZtoY4P(Z):
    Y=np.empty(np.shape(Z),dtype=np.complex128)
    for i in range(len(Z)):
            ZZ=Z[i]
            Y[i]=np.linalg.inv(ZZ) 
    return Y

def YtoZ4P(Y):
    Z=np.empty(np.shape(Y),dtype=np.complex128)
    for i in range(len(Y)):
            YY=Y[i]
            Z[i]=np.linalg.inv(YY) 
#            try:
#                Z[i]=np.linalg.inv(YY)
#            except np.linalg.linalg.LinAlgError as err:
#                if 'Singular matrix' in err.message:
#                     print "Singular in ",
#                     print i
#                    II=np.diag([1.0, 1.0, 1.0, 1.0])*1e-18
#                    YYY=YY+II
#                    Z[i]=np.linalg.inv(YYY)
#                else:
#                    raise
    return Z

def Z4PtoH4(Z):
    #
    #  |V1|    | H11 H12 H13 H14 | |I1| !!G
    #  |I2|  - | H21 H22 H23 H24 | |V2| !!D 
    #  |I3|  - | H31 H32 H33 H34 | |V3| !!S
    #  |V4|    | H41 H42 H43 H44 | |I4| !!BG
    #
    H=np.empty(np.shape(Z),dtype=np.complex128)
    for i in range(len(Z)):
            ZZ=Z[i]
            z11=ZZ[0,0]
            z12=ZZ[0,1]
            z13=ZZ[0,2]
            z14=ZZ[0,3]
            z21=ZZ[0,0]
            z22=ZZ[1,1]
            z23=ZZ[1,2]
            z24=ZZ[1,3]
            z31=ZZ[2,0]
            z32=ZZ[2,1]
            z33=ZZ[2,2]
            z34=ZZ[2,3]
            z41=ZZ[3,0]
            z42=ZZ[3,1]
            z43=ZZ[3,2]
            z44=ZZ[3,3]
            H[i,0,0]=(z13*z22*z31 - z12*z23*z31 - z13*z21*z32 + z11*z23*z32 + z12*z21*z33 - z11*z22*z33)/(z23*z32 - z22*z33)
            H[i,0,1]=(z13*z32 - z12*z33)/(z23*z32 - z22*z33)
            H[i,0,2]=(z13*z22 - z12*z23)/(-(z23*z32) + z22*z33)
            H[i,0,3]=(z14*z23*z32 - z13*z24*z32 - z14*z22*z33 + z12*z24*z33 + z13*z22*z34 - z12*z23*z34)/(z23*z32 - z22*z33)
            H[i,1,0]=(z23*z31 - z21*z33)/(-(z23*z32) + z22*z33)
            H[i,1,1]=z33/(-(z23*z32) + z22*z33)
            H[i,1,2]=z23/(z23*z32 - z22*z33)
            H[i,1,3]=(z24*z33 - z23*z34)/(z23*z32 - z22*z33)
            H[i,2,0]=(z22*z31 - z21*z32)/(z23*z32 - z22*z33)
            H[i,2,1]=z32/(z23*z32 - z22*z33)
            H[i,2,2]=z22/(-(z23*z32) + z22*z33)
            H[i,2,3]=(z24*z32 - z22*z34)/(-(z23*z32) + z22*z33)
            H[i,3,0]=(z23*z32*z41 - z22*z33*z41 - z23*z31*z42 + z21*z33*z42 + z22*z31*z43 - z21*z32*z43)/(z23*z32 - z22*z33)
            H[i,3,1]=(z33*z42 - z32*z43)/(-(z23*z32) + z22*z33) 
            H[i,3,2]=(z23*z42 - z22*z43)/(z23*z32 - z22*z33) 
            H[i,3,3]=(z24*z33*z42 - z23*z34*z42 - z24*z32*z43 + z22*z34*z43 + z23*z32*z44 - z22*z33*z44)/(z23*z32 - z22*z33) 
    return H

def Y4PtoH4(Y):
    #
    #  |V1|    | H11 H12 H13 H14 | |I1| G
    #  |I2|  - | H21 H22 H23 H24 | |V2| D
    #  |I3|  - | H31 H32 H33 H34 | |V3| S
    #  |V4|    | H41 H42 H43 H44 | |I4| BG
    #
    H=np.empty(np.shape(Y),dtype=np.complex128)
    for i in range(len(Y)):
            YY=Y[i]
            y11=YY[0,0]
            y12=YY[0,1]
            y13=YY[0,2]
            y14=YY[0,3]
            y21=YY[0,0]
            y22=YY[1,1]
            y23=YY[1,2]
            y24=YY[1,3]
            y31=YY[2,0]
            y32=YY[2,1]
            y33=YY[2,2]
            y34=YY[2,3]
            y41=YY[3,0]
            y42=YY[3,1]
            y43=YY[3,2]
            y44=YY[3,3]
            H[i,0,0]=y44/(-(y14*y41) + y11*y44)
            H[i,0,1]=(y14*y42 - y12*y44)/(-(y14*y41) + y11*y44)
            H[i,0,2]=(y14*y43 - y13*y44)/(-(y14*y41) + y11*y44)
            H[i,0,3]=y14/(y14*y41 - y11*y44)
            H[i,1,0]=(y24*y41 - y21*y44)/(y14*y41 - y11*y44)
            H[i,1,1]=(y14*y22*y41 - y12*y24*y41 - y14*y21*y42 + y11*y24*y42 + y12*y21*y44 - y11*y22*y44)/(y14*y41 - y11*y44)
            H[i,1,2]=(y14*y23*y41 - y13*y24*y41 - y14*y21*y43 + y11*y24*y43 + y13*y21*y44 - y11*y23*y44)/(y14*y41 - y11*y44)
            H[i,1,3]=(y14*y21 - y11*y24)/(y14*y41 - y11*y44)
            H[i,2,0]=(y34*y41 - y31*y44)/(y14*y41 - y11*y44)
            H[i,2,1]=(y14*y32*y41 - y12*y34*y41 - y14*y31*y42 + y11*y34*y42 + y12*y31*y44 - y11*y32*y44)/(y14*y41 - y11*y44)
            H[i,2,2]=(y14*y33*y41 - y13*y34*y41 - y14*y31*y43 + y11*y34*y43 + y13*y31*y44 - y11*y33*y44)/(y14*y41 - y11*y44)
            H[i,2,3]=(y14*y31 - y11*y34)/(y14*y41 - y11*y44)
            H[i,3,0]=y41/(y14*y41 - y11*y44)
            H[i,3,1]=(y12*y41 - y11*y42)/(-(y14*y41) + y11*y44)
            H[i,3,2]=(y13*y41 - y11*y43)/(-(y14*y41) + y11*y44)
            H[i,3,3]=y11/(-(y14*y41) + y11*y44)
    return H

def Y4PtoH3(Y):
    #
    #  |V1|    | H11 H12 H13 H14 | |I1| G
    #  |I2|  - | H21 H22 H23 H24 | |V2| D
    #  |V3|  - | H31 H32 H33 H34 | |I3| BG
    #  |I4|    | H41 H42 H43 H44 | |V4| S
    #
    H=np.empty(np.shape(Y),dtype=np.complex128)
    for i in range(len(Y)):
            YY=Y[i]
            y11=YY[0,0]
            y12=YY[0,1]
            y13=YY[0,2]
            y14=YY[0,3]
            y21=YY[0,0]
            y22=YY[1,1]
            y23=YY[1,2]
            y24=YY[1,3]
            y31=YY[2,0]
            y32=YY[2,1]
            y33=YY[2,2]
            y34=YY[2,3]
            y41=YY[3,0]
            y42=YY[3,1]
            y43=YY[3,2]
            y44=YY[3,3]
            H[i,0,0]=y44/(-(y14*y41) + y11*y44)
            H[i,0,1]=(y14*y42 - y12*y44)/(-(y14*y41) + y11*y44)
            H[i,0,2]=(y14*y43 - y13*y44)/(-(y14*y41) + y11*y44)
            H[i,0,3]=y14/(y14*y41 - y11*y44)
            H[i,1,0]=(y24*y41 - y21*y44)/(y14*y41 - y11*y44)
            H[i,1,1]=(y14*y22*y41 - y12*y24*y41 - y14*y21*y42 + y11*y24*y42 + y12*y21*y44 - y11*y22*y44)/(y14*y41 - y11*y44)
            H[i,1,2]=(y14*y23*y41 - y13*y24*y41 - y14*y21*y43 + y11*y24*y43 + y13*y21*y44 - y11*y23*y44)/(y14*y41 - y11*y44)
            H[i,1,3]=(y14*y21 - y11*y24)/(y14*y41 - y11*y44)
            H[i,2,0]=(y34*y41 - y31*y44)/(y14*y41 - y11*y44)
            H[i,2,1]=(y14*y32*y41 - y12*y34*y41 - y14*y31*y42 + y11*y34*y42 + y12*y31*y44 - y11*y32*y44)/(y14*y41 - y11*y44)
            H[i,2,2]=(y14*y33*y41 - y13*y34*y41 - y14*y31*y43 + y11*y34*y43 + y13*y31*y44 - y11*y33*y44)/(y14*y41 - y11*y44)
            H[i,2,3]=(y14*y31 - y11*y34)/(y14*y41 - y11*y44)
            H[i,3,0]=y41/(y14*y41 - y11*y44)
            H[i,3,1]=(y12*y41 - y11*y42)/(-(y14*y41) + y11*y44)
            H[i,3,2]=(y13*y41 - y11*y43)/(-(y14*y41) + y11*y44)
            H[i,3,3]=y11/(-(y14*y41) + y11*y44)
    return H

################################################################################

def findHz(fHz, fstart, fstop=0.0):
    if (fstart == fstop) or (fstop==0.0): # we ask for a precise value
       if (fHz[0]== fstart): # Azz it is the lower bound
            return 0, 1
       elif (fHz[-1]== fstart): # Azz it is the upper bound
            return len(fHz)-1, len(fHz)
       elif (fHz[0]<= fstart and fHz[-1]>= fstart): # if it is included
           op=-1
           for i  in range(len(fHz)-1):
              if fHz[i] <= fstart and fHz[i+1]> fstart:
                 op=i
           return op, op+1
       else:  # the asked value is outside the range. return all the range
           return 0, len(fHz) 
    elif (fstop > fstart):
       ar=0
       op=len(fHz)
       for i in range(len(fHz)-1):
          low=fHz[i]
          hi =fHz[i+1]
          if (fstart >= low and fstart<hi):
             ar=i
          if (fstop > low and fstop<=hi):
             op=i+2
       if op>=len(fHz):
          op=len(fHz)
       return ar, op 
    else:
       return 0,len(fHz)


################################################################################
##### 2 Port extractions #######################################################

def Y2PtoGds(Y, fHz, fstart=10, fstop=1000e9):
    start, stop = findHz(fHz, fstart, fstop)
    Xv=(np.real(Y[:,1,1]))
    X=sum(Xv[start:stop])/len(Xv[start:stop])
    return X, Xv

def Y2PtoLgAvg(Y, fHz, fstart=10, fstop=1000e9):
    start, stop = findHz(fHz, fstart, fstop)
    ZZ=YtoZ_2P(Y)
    Xv=(np.imag(ZZ[:,0,0]-ZZ[:,0,1]))/(fHz*2*np.pi)
    X=sum(Xv[start:stop])/len(Xv[start:stop])
    return X, Xv

def Y2PtoLdAvg(Y, fHz, fstart=10, fstop=1000e9):
    start, stop = findHz(fHz, fstart, fstop)
    ZZ=YtoZ_2P(Y)
    Xv=(np.imag(ZZ[:,1,1]-ZZ[:,0,1]))/(fHz*2*np.pi)
    X=sum(Xv[start:stop])/len(Xv[start:stop])
    return X, Xv

def Y2PtoLsAvg(Y, fHz, fstart=10, fstop=1000e9):
    start, stop = findHz(fHz, fstart, fstop)
    ZZ=YtoZ_2P(Y)
    Xv=(np.imag(ZZ[:,0,1]))/(fHz*2*np.pi)
    X=sum(Xv[start:stop])/len(Xv[start:stop])
    return X, Xv

def Y2PtoLg(Y, fHz, fstart=10, fstop=1000e9):
    start, stop = findHz(fHz, fstart, fstop)
    ZZ=YtoZ_2P(Y)
    Xv=(np.imag(ZZ[:,0,0]-ZZ[:,0,1]))/(fHz*2*np.pi)
    X=Xv[start:stop]
    denom=1/((fHz*2*np.pi)**2)
    return X, Xv, denom

def Y2PtoLd(Y, fHz, fstart=10, fstop=1000e9):
    start, stop = findHz(fHz, fstart, fstop)
    ZZ=YtoZ_2P(Y)
    Xv=(np.imag(ZZ[:,1,1]-ZZ[:,0,1]))/(fHz*2*np.pi)
    X=Xv[start:stop]
    return X, Xv

def Y2PtoLs(Y, fHz, fstart=10, fstop=1000e9):
    start, stop = findHz(fHz, fstart, fstop)
    ZZ=YtoZ_2P(Y)
    Xv=(np.imag(ZZ[:,0,1]))/(fHz*2*np.pi)
    X=Xv[start:stop]
    return X, Xv

def Y2PtoRgBracale(Y, fHz, fstart=10, fstop=1000e9):
    start, stop = findHz(fHz, fstart, fstop)
    ZZ=YtoZ_2P(Y)
    Xv=(np.real(ZZ[:,0,0]-ZZ[:,0,1]))
    X=sum(Xv[start:stop])/len(Xv[start:stop])
    return X, Xv

def Y2PtoRgEnz(Y, fHz, fstart=10, fstop=1000e9):
    start, stop = findHz(fHz, fstart, fstop)
    Xv=(np.real(Y[:,0,0])/((np.imag(Y[:,0,0]))**2))
    X=Xv[start:stop]
    return X, Xv

def Y2PtoRgJen(Y, fHz, fstart=10, fstop=1000e9):
    start, stop = findHz(fHz, fstart, fstop) 
    Xv=np.real(Y[:,0,1])/((np.imag(Y[:,0,0]))*(np.imag(Y[:,0,1])))
    X=Xv[start:stop]
    return X, Xv

def Y2PtoRgDormieu(Y, fHz, fstart=10, fstop=1000e9):
    start, stop = findHz(fHz, fstart, fstop) 
    Xv=np.real(1/Y[:,0,0])
    X=Xv[start:stop]
    return X, Xv

def Y2PtoRgEnzAvg(Y, fHz, fstart=10, fstop=1000e9):
    start, stop = findHz(fHz, fstart, fstop)
    Xv=(np.real(Y[:,0,0])/((np.imag(Y[:,0,0]))**2))
    X=sum(Xv[start:stop])/len(Xv[start:stop])
    return X, Xv

def Y2PtoRgJenAvg(Y, fHz, fstart=10, fstop=1000e9):
    start, stop = findHz(fHz, fstart, fstop) 
    Xv=np.real(Y[:,0,1])/((np.imag(Y[:,0,0]))*(np.imag(Y[:,0,1])))
    X=sum(Xv[start:stop])/len(Xv[start:stop])
    return X, Xv

def Y2PtoRgDormieuAvg(Y, fHz, fstart=10, fstop=1000e9):
    start, stop = findHz(fHz, fstart, fstop) 
    Xv=np.real(1/Y[:,0,0])
    X=sum(Xv[start:stop])/len(Xv[start:stop])
    return X, Xv


def Y2PtoRd(Y, fHz, fstart=10, fstop=1000e9):
    start, stop = findHz(fHz, fstart, fstop)
    ZZ=YtoZ_2P(Y)
    Xv=(np.real(ZZ[:,1,1]-ZZ[:,0,1]))
    X=sum(Xv[start:stop])/len(Xv[start:stop])
    return X, Xv

def Y2PtoRs(Y, fHz, fstart=10, fstop=1000e9):
    start, stop = findHz(fHz, fstart, fstop)
    ZZ=YtoZ_2P(Y)
    Xv=(np.real(ZZ[:,0,1]))
    X=sum(Xv[start:stop])/len(Xv[start:stop])
    return X, Xv

def Y2PtoRs1(Y, fHz, fstart=10, fstop=1000e9):
    start, stop = findHz(fHz, fstart, fstop)
    ZZ=YtoZ_2P(Y)
    Xv=(np.real(ZZ[:,0,1]))
    X=Xv[start:stop]
    return X, Xv

def Y2PtoRs2(Y, fHz, fstart=10, fstop=1000e9):
    start, stop = findHz(fHz, fstart, fstop)
    ZZ=YtoZ_2P(Y)
    Xv=(np.real(ZZ[:,1,0]))
    X=Xv[start:stop]
    return X, Xv

def Y2PtoCopen1(Y, fHz, fstart=10, fstop=1000e9):
    start, stop = findHz(fHz, fstart, fstop)
    Copenv=(np.imag(Y[:,0,0])/(fHz*2*np.pi))
    Copen=Copenv[start:stop]
    return Copen, Copenv

def Y2PtoCopen2(Y, fHz, fstart=10, fstop=1000e9):
    start, stop = findHz(fHz, fstart, fstop)
    Copenv=(np.imag(Y[:,1,1])/(fHz*2*np.pi))
    Copen=Copenv[start:stop]
    return Copen, Copenv

def Y2PtoC12(Y, fHz, fstart=10, fstop=100e9):
    start, stop = findHz(fHz, fstart, fstop)
    Cgsv=(np.imag(Y[:,0,1]))/(fHz*2*np.pi)
    Cgs=Cgsv[start:stop]
    return Cgs, Cgsv

def Y2PtoC21(Y, fHz, fstart=10, fstop=100e9):
    start, stop = findHz(fHz, fstart, fstop)
    Cgsv=(np.imag(Y[:,1,0]))/(fHz*2*np.pi)
    Cgs=Cgsv[start:stop]
    return Cgs, Cgsv

def Y2PtoCgs(Y, fHz, fstart=10, fstop=100e9):
    start, stop = findHz(fHz, fstart, fstop)
    Cgsv=(np.imag(Y[:,0,0]+Y[:,0,1]))/(fHz*2*np.pi)
    Cgs=sum(Cgsv[start:stop])/len(Cgsv[start:stop])
    return Cgs, Cgsv

def Y2PtoCgd(Y, fHz, fstart=10, fstop=100e9):
    start, stop = findHz(fHz, fstart, fstop)
    Cgdv=-(np.imag(Y[:,0,1]))/(fHz*2*np.pi)
    Cgd=sum(Cgdv[start:stop])/len(Cgdv[start:stop])
    return Cgd, Cgdv

def Y2PtoGm(Y, fHz, fstart=10, fstop=1000e9):
    start, stop = findHz(fHz, fstart, fstop)
    Xv=(np.real(Y[:,1,0]-Y[:,0,1]))
    X=sum(Xv[start:stop])/len(Xv[start:stop])
    return X, Xv

def Y2PtoGmAlt(Y, fHz, fstart=10, fstop=1000e9):
    # da usare in caso ci sia resistenza di Gate
    start, stop = findHz(fHz, fstart, fstop)
    Xnum=np.abs(Y[:,1,0]-Y[:,0,1])
    Xden=np.abs(Y[:,0,0]+Y[:,0,1])
    Xima=np.imag(-1/(Y[:,0,0]+Y[:,0,1]))
    Xv=np.abs(Xnum/Xden)/Xima
    Xv=(np.real(Y[:,1,0]-Y[:,0,1]))
    X=sum(Xv[start:stop])/len(Xv[start:stop])
    return X, Xv

#####3 GAIN and STABILITY for a two port

def S2PUMasondB(S):
    S11=S[:,0,0]
    S22=S[:,1,1]
    S12=S[:,0,1]
    S21=S[:,1,0]
    Ratio=S21/S12
    Det=S11*S22-S12*S21
    kstb=1-np.abs(S11)**2-np.abs(S22)**2+np.abs(Det)**2
    kstb=kstb/(2*np.abs(S12*S21))
    mason = (np.abs(Ratio-1)**2)/(2*kstb*np.abs(Ratio)-2*np.real(Ratio))
    return 10*np.log10(mason)

def S2PMAGdB(S):
    S11=S[:,0,0]
    S22=S[:,1,1]
    S12=S[:,0,1]
    S21=S[:,1,0]
    Ratio=S21/S12
    Det=S11*S22-S12*S21
    kstb=1-np.abs(S11)**2-np.abs(S22)**2+np.abs(Det)**2
    kstb=kstb/(2*np.abs(S12*S21))
    mag = np.abs(Ratio)*(kstb-np.sqrt(kstb**2-1)) # l'articolo Gupta mette
		    # k +- sqrt(k^2-1)  ma la maggior parte dei paper prende
		    # solo la version - radice non quella + radice
		    # inoltre e' normale che sia NaN se k<1
    return 10*np.log10(mag)

def S2PMEGdB(S):
    S11=S[:,0,0]
    S22=S[:,1,1]
    S12=S[:,0,1]
    S21=S[:,1,0]
    Ratio=S21/S12
    Det=S11*S22-S12*S21
    kstb=1-np.abs(S11)**2-np.abs(S22)**2+np.abs(Det)**2
    kstb=kstb/(2*np.abs(S12*S21))
    meg = (np.abs(Ratio)**2-1)/(2*kstb*np.abs(Ratio)-1)
    return 10*np.log10(meg)

def S2PMSGdB(S):
    S11=S[:,0,0]
    S22=S[:,1,1]
    S12=S[:,0,1]
    S21=S[:,1,0]
    Ratio=S21/S12
    return 10*np.log10(np.abs(Ratio))

def S2Pstability(S):
    S11=S[:,0,0]
    S22=S[:,1,1]
    S12=S[:,0,1]
    S21=S[:,1,0]
    Det=S11*S22-S12*S21
    kstb=1-np.abs(S11)**2-np.abs(S22)**2+np.abs(Det)**2
    kstb=kstb/(2*np.abs(S12*S21))
    return Det, kstb

def Y2PRolletk(Y, fHz):
    # aggiunto solo per testare pedissecuamente
    # che estrarlo da S o Y da lo stesso risultato
    Y11=Y[:,0,0]
    Y22=Y[:,1,1]
    Y12=Y[:,0,1]
    Y21=Y[:,1,0]
    num=2.0*np.real(Y11)*np.real(Y22)
    num=num-np.real(Y12*Y21)
    Xv=num/np.abs(Y12*Y21)
    X=min(Xv)
    return Xv

def Y2PMAGdB(Y, fHz):
    num1=2.0*np.real(Y[:,0,0])*np.real(Y[:,1,1])
    num2=num1-np.real(Y[:,1,0]*Y[:,0,1])
    ksf=num2/np.abs(Y[:,1,0]*Y[:,0,1])
    num3=np.abs(Y[:,1,0]/Y[:,0,1])
    num=num3
    den=ksf+np.sqrt(ksf*ksf-1)
    return 10*np.log10(num/den)

def Y2PMSGdB(Y,fHz):
    return  10*np.log10(np.abs(Y[:,1,0]/Y[:,0,1]))

def Y2PMEGdB(Y,fHz):
    Y11=Y[:,0,0]
    Y22=Y[:,1,1]
    Y12=Y[:,0,1]
    Y21=Y[:,1,0]
    Ya=Y11*Y22
    Yb=Y12*Y21
    num=np.abs(Y21)**2-np.abs(Y12)**2
    den=4*np.real(Ya)-2*np.real(Yb)-2*np.abs(Y12)**2
    return 10*np.log10(num/den)

def Y2PMUGdB(Y,fHz):
    Y11=Y[:,0,0]
    Y22=Y[:,1,1]
    Y12=Y[:,0,1]
    Y21=Y[:,1,0]
    num=Y12-Y21
    num=num*np.conjugate(num)
    den=4*np.real(Y11)*np.real(Y22)
    den=den-4*np.real(Y12)*np.real(Y21)
    return 10*np.log10(num/den)

def ADS_stab_fact(S):
   # Returns the Rollett stability factor.
   # RETURNED_VALUE: Real
   # CATEGORY: S-Parameter
   # EXAMPLE: k = stab_fact(S)
   # SYNTAX: k = stab_fact(S)
   # ARGUMENT
   #   ARG_NAME: S
   #   ARG_DESCRIPTION: scattering matrix of 2-port network
   #   ARG_DEFAULT: None
   #   ARG_RANGE: (-inf:inf)
   #   ARG_TYPE: Complex
   #  ARG_REQUIRED: YES
   # NOTES: Given a 2 x 2 scattering matrix between the input and measurement ports,
   # this function calculates the stability factor.
   # The Rollett stability factor is given by
   # k = {1- |S11|**2 - |S22|**2 + |S11*S22 - S12*S21| **2} / {2*|S12*S21|}
   # The necessary and sufficient conditions for unconditional stability are that
   # the stability factor is greater than unity and the stability measure is positive.
   # Reference
   # [1] Guillermo Gonzales, Microwave Transistor Amplifiers, second edition, Prentice-Hall, 1997.
   S11=S[:,0,0]
   S22=S[:,1,1]
   S12=S[:,0,1]
   S21=S[:,1,0]
   a1 = np.abs(S11)*np.abs(S11)+np.abs(S22)*np.abs(S22)
   a2 = np.abs(S11*S22-S12*S21)*np.abs(S11*S22-S12*S21)
   numer = 1.0-a1+a2
   denom = 2*np.abs(S12*S21)
   return numer/denom

def ADS_stab_meas(S):
   # Returns the stability measure
   # RETURNED_VALUE: Real
   # CATEGORY: S-Parameter
   # EXAMPLE: k = stab_meas(S)
   # SYNTAX: k = stab_meas(S)
   # NOTES: Given a 2 x 2 scattering matrix between the input and measurement ports,
   # this function calculates the stability measure.
   # The stability measure is given by
   # b = 1+ |S11|**2 - |S22|**2 - |S11*S22 - S12*S21| **2
   # The necessary and sufficient conditions for unconditional stability are that
   # the stability factor is greater than unity and the stability measure is positive.
   # Reference
   # Guillermo Gonzales, Microwave Transistor Amplifiers, second edition, Prentice-Hall, 1997.
   S11=S[:,0,0]
   S22=S[:,1,1]
   S12=S[:,0,1]
   S21=S[:,1,0]    
   a1 = np.abs(S11)*np.abs(S11) - np.abs(S22)*np.abs(S22)
   a2 = np.abs(S11*S22-S12*S21) * np.abs(S11*S22-S12*S21)
   return (1.0+a1-a2)


def ADS_max_gain(S):
   ##  Given a 2 x 2 scattering matrix, this measurement returns the maximum available and stable gain (in dB) between the input and the measurement ports.
   ## RETURNED_VALUE: Real
   ## CATEGORY: S-Parameter
   ## EXAMPLE: y = max_gain(S)
   # SYNTAX:  y = max_gain(S)
   # ARGUMENT
   #  ARG_NAME: S
   #  ARG_DESCRIPTION: scattering matrix of 2-port network
   #  ARG_DEFAULT: None
   #  ARG_RANGE: (-inf:inf)
   #  ARG_TYPE: Complex
   #  ARG_REQUIRED: YES
   S11=S[:,0,0]
   S22=S[:,1,1]
   S12=S[:,0,1]
   S21=S[:,1,0]
   magS21=np.abs(S21)
   magS12=np.abs(S12)
   magSqS11=np.abs(S11)*np.abs(S11)
   magSqS22=np.abs(S22)*np.abs(S22)
   Bstab = ADS_stab_meas(S)
   k = ADS_stab_fact(S)
   k[k<=1]=1.0 # cancella i valori di k<1 e li sostituisce con 1
   gain0  = (magS21*magS21)/((1.0 - magSqS11)*(1.0 - magSqS22)) # 
   gain1  = (magS21/magS12)*(k - np.sqrt(k*k-1.0)) # MAG + 
   gain2  = (magS21/magS12)*(k + np.sqrt(k*k-1.0)) # MAG - 
   gainMSG =(magS21/magS12)
   maxGain = np.empty(np.shape(gain0), dtype=np.float_)
   for f in range(len(S11)):
      if magS12[f]==0: 
         tmp = gain0[f]
      else: # magS12[f] > 0
        if Bstab[f] > 0 :
            tmp = gain1[f]
        else:
            tmp = gain2[f]
      if tmp < 1e-304:
        tmp=1e-304
      maxGain[f]=tmp
   return 10.0*np.log10(maxGain)

def Unilateral_FOM(S):
# returns the unilateral figure of merit as a real(f)
# Gonzalez pg. 239
# http://edadocs.software.keysight.com/pages/viewpage.action?pageId=5920607
# Se Gt is the transducer power gain
# Se Gtu is the trasducer power gain of the unilateralized device
# then
# U is the unilateral FOM
# and 
#
#                  |   1    |    Gt    |   1    |
#   U(plus)=10 log |--------| < ---- < |--------|= U(minus)
#                  | (1-U)^2|   Gtu    | (1+U)^2|
#
# Used in Small-signal S-parameter simulations.
# This FOM determines whether the simplification can be 
# made in neglecting the effect of S12 
# (unilateral behavior of device).
    S11a=np.abs(S[:,0,0])
    S22a=np.abs(S[:,1,1])
    S12a=np.abs(S[:,0,1])
    S21a=np.abs(S[:,1,0])
    U=S11a*S12a*S21a*S22a
    U=U/(1-S11a*S11a)
    U=U/(1-S22a*S22a)
    Uminus=10*np.log10(1/((1-U)*(1-U)))
    Uplus =10*np.log10(1/((1+U)*(1+U)))
    return U, Uminus, Uplus
# to test:
# SS=np.zeros((1,2,2),dtype=np.complex_)
# d2r = lambda x: np.deg2rad(x)
# r2d = lambda x: np.rad2deg(x)
# polar = lambda x,y: x*np.exp(1j*d2r(y))
# SS[0]=((polar(0.55,-50), polar(0.02,10)), (polar(3.82,80),polar(0.15,-20)))
# U, Up, Um=rftools.Unilateral_FOM(SS)
# print "U", U[0]
# print "U- (dB)" , Um[0]
# print "U+ (dB)" , Up[0]
# deve dare: U returns 0.009, Uplus returns 0.081 dB Uminus returns -0.08 dB
# so The error in using Gt is less than 0.1 dB 
# and the device is approximatively unilateral

def Y4PtoY2P(Y, Port=3):
    # for the BG in some case simply reduce the Y matrix fronm 4 to 2 ports
    # riduce le matrice Y a 4 porte ad una identica a 2 Porte.
    # chiaramente conservando drain e Back-Gate
    # per poter applicare le trasformate dei due porte al BGate
    
    PG=Port-1 # Port=1 -> PG=0 Port=4 -> PG=3 and Port=4 -> PG=2
    Y2P=np.empty((len(Y),2,2),dtype=np.complex128)
    for i in range(len(Y)):
            YY=Y[i]
            Y2P[i,0,0]=YY[PG,PG]
            Y2P[i,0,1]=YY[PG,1]
            Y2P[i,1,0]=YY[1,PG]
            Y2P[i,1,1]=YY[1,1]
    return Y2P

###################################################################
#### Pavageau - bracale

def Y2PtoCpg(Y, fHz, fstart=10, fstop=1000e9):
    # function used to extract residual Capacitance to ground
    # from Port1 and Port2. Use on a zero bias (Vds=0V, Vg=0V)
    # device. Extraction as in Pavageau PhD Thesis, pg 88.
    start, stop = findHz(fHz, fstart, fstop)
    Copenv=(np.imag(Y[:,0,0]+2*Y[:,0,1])/(fHz*2*np.pi))
    Copen=sum(Copenv[start:stop])/len(Copenv[start:stop])
    return Copen, Copenv

def Y2PtoCpd(Y, fHz, fstart=10, fstop=1000e9):
    start, stop = findHz(fHz, fstart, fstop)
    Copenv=(np.imag(Y[:,1,1]+Y[:,0,1])/(fHz*2*np.pi))
    Copen=sum(Copenv[start:stop])/len(Copenv[start:stop])
    return Copen, Copenv

def Y2PtoCb(Y, fHz, fstart=10, fstop=1000e9):
    start, stop = findHz(fHz, fstart, fstop)
    Copenv=(np.imag(-Y[:,0,1])/(fHz*2*np.pi))
    Copen=sum(Copenv[start:stop])/len(Copenv[start:stop])
    return Copen, Copenv


################################################################################
#### 4-port extractions


def Y4Ptoft(Y, fHz, fstart=10, fstop=1000e9):
    start, stop = findHz(fHz, fstart, fstop)
    # suppone che FG sia Port 1, drain port 2
    num = Y[:,1,0]
    den = Y[:,0,0]
    h21=num/den
    ftv=(fHz*np.abs(h21))
    ft=sum(ftv[start:stop])/1e9/len(ftv[start:stop])
    return ft, ftv


def Y4Ptoftbg(Y, fHz, Port=3, fstart=10, fstop=1000e9, method=1, clean=False):
    ft_start, ft_stop = findHz(fHz, fstart, fstop)
    if Port==4:
       num = Y[:,1,3]
       den = Y[:,3,3]
    else:
       num = Y[:,1,2]
       den = Y[:,2,2]
    hdg = np.abs(num/den)
    ftbgv=(fHz*hdg)
    nn=len(hdg)-1 ## nn index of last value

    zerocrossing=False
    for i in range(nn,1,-1): # zc first value <=1
       if hdg[i-1] > 1 and hdg[i]<=1:
          zc=i
          zerocrossing=True
          break

    if not zerocrossing and hdg[nn]<=1: # the whole curve is <1, no gain
       if clean:
          return 0, ftbgv[0:2], hdg
       else:
          return 0, ftbgv, hdg

    if not zerocrossing and hdg[nn]>=1: # the whole curve is > 1
       zc=nn # we need it instantiated

    if clean and zerocrossing and zc<nn: # there is crossing
       # hdg=hdg[:(zc+1)]
       fHz=fHz[:(zc+1)]    
       ftbgv=ftbgv[:(zc+1)]

    if ft_stop>=0:
        ft_stop=min(ft_stop,zc)
    else:
        ft_stop=min(nn+1+ft_stop,zc)

    if ft_stop<=ft_start : # bad interval!
       ftbg=0
    else: 
       ftbg=sum(ftbgv[ft_start:ft_stop])/1e9/len(ftbgv[ft_start:ft_stop])

    if method==2 and zerocrossing and zc<nn: # ricalcola ftbg con l'intercetta, se l'intercetta c'e'
       lhdg=np.log10(hdg)
       lf=np.log10(fHz)
       ftbg=lf[zc-1]+(lf[zc]-lf[zc-1])*(-lhdg[zc-1]/(lhdg[zc]-lhdg[zc-1]))
       ftbg=(10**ftbg)/1e9
    return ftbg, ftbgv, hdg

def Y4Ptoftbg4_old_wrong(Y, fHz, ft_start=0, ft_stop=-1):
    # suppone che BG sia Port 4, drain port 2
    num = Y[:,1,3]*Y[:,0,0]-Y[:,1,0]*Y[:,0,3] # 24 11 - 21 14
    den = Y[:,3,3]*Y[:,0,0]-Y[:,0,3]*Y[:,3,0] # 44 11 - 14 41
    h24 = num/den
    ftbgv1=(fHz*np.abs(h24))
    ftbg1=sum(ftbgv1[ft_start:ft_stop])/1e9/len(ftbgv1[ft_start:ft_stop])
    return ftbg1, ftbgv1

def Y4Ptofmaxbg(Y, fHz, Port=3, fstart=10, fstop=1000e9, method=1, clean=False):
    ft_start, ft_stop = findHz(fHz, fstart, fstop)
    if Port==4:
       # suppone che BG sia Port 4, drain port 2
       num=np.abs(Y[:,1,3]-Y[:,3,1])**2
       den=4*(Y[:,3,3].real*Y[:,1,1].real-Y[:,3,1].real*Y[:,1,3].real)
    else:
       # suppone che BG sia Port 3, drain port 2
       num=np.abs(Y[:,1,2]-Y[:,2,1])**2
       den=4*(Y[:,2,2].real*Y[:,1,1].real-Y[:,2,1].real*Y[:,1,2].real)
   
    sqU=np.sqrt(np.abs(num/den)) 
    fmaxbgv=(fHz*sqU)    
    nn=len(sqU)-1 ## nn index of last value

    zerocrossing=False
    for i in range(nn,1,-1): # zc first value <=1
       if sqU[i-1] > 1 and sqU[i]<=1:
          zc=i
          zerocrossing=True
          break

    if not zerocrossing and sqU[nn]<=1: # the whole curve is <1, no gain
       if clean:
          return 0, fmaxbgv[0:2], sqU
       else:
          return 0, fmaxbgv, sqU

    if not zerocrossing and sqU[nn]>=1: # the whole curve is > 1
       zc=nn # we need it instantiated

    if clean and zerocrossing and zc<nn: # there is crossing
       # sqU=sqU[:(zc+1)]
       fHz=fHz[:(zc+1)]    
       fmaxbgv=fmaxbgv[:(zc+1)]

    if ft_stop>=0:
        ft_stop=min(ft_stop,zc)
    else:
        ft_stop=min(nn+1+ft_stop,zc)

    if ft_stop<=ft_start : # bad interval!
       fmaxbg=0
    else: 
       fmaxbg=sum(fmaxbgv[ft_start:ft_stop])/1e9/len(fmaxbgv[ft_start:ft_stop])

    if method==2 and zerocrossing and zc<nn: # ricalcola ftbg con l'intercetta, se l'intercetta c'e'
       lsqU=np.log10(sqU)
       lf=np.log10(fHz)
       fmaxbg=lf[zc-1]+(lf[zc]-lf[zc-1])*(-lsqU[zc-1]/(lsqU[zc]-lsqU[zc-1]))
       fmaxbg=(10**fmaxbg)/1e9
    return fmaxbg, fmaxbgv, sqU
    
def Y4Ptofmax(Y, fHz, fstart=10, fstop=1000e9):
    start, stop = findHz(fHz, fstart, fstop)
    num=np.abs(Y[:,1,0]-Y[:,0,1])**2
    den=4*(Y[:,0,0].real*Y[:,1,1].real-Y[:,0,1].real*Y[:,1,0].real)
    fmaxv=(fHz*np.sqrt(np.abs(num/den)))    
    fmax=sum(fmaxv[start:stop])/1e9/len(fmaxv[start:stop])
    return fmax, fmaxv

def Y4Ptofmaxbg_old_but_correct(Y, fHz, Port=3, fstart=10, fstop=1000e9):
    if Port==4:
       (s, v)=Y4Ptofmaxbg4(Y, fHz, fstart=fstart, fstop=fstop)
       return s, v
    else:
       (s, v)=Y4Ptofmaxbg3(Y, fHz, fstart=fstart, fstop=fstop)
       return s, v

def Y4Ptofmaxbg3(Y, fHz, fstart=10, fstop=100e9):
    # suppone che BG sia Port 3, drain port 2
    start, stop = findHz(fHz, fstart, fstop)
    num=np.abs(Y[:,1,2]-Y[:,2,1])**2
    den=4*(Y[:,2,2].real*Y[:,1,1].real-Y[:,2,1].real*Y[:,1,2].real)
    fmaxv=(fHz*np.sqrt(np.abs(num/den)))    
    fmax=sum(fmaxv[start:stop])/1e9/len(fmaxv[start:stop])
    return fmax, fmaxv

def Y4Ptofmaxbg4(Y, fHz, fstart=10, fstop=100e9):
    start, stop = findHz(fHz, fstart, fstop)
    # suppone che BG sia Port 4, drain port 2
    # recentemente cambiato i limiti start/stop per passare in Hz
    num=np.abs(Y[:,1,3]-Y[:,3,1])**2
    den=4*(Y[:,3,3].real*Y[:,1,1].real-Y[:,3,1].real*Y[:,1,3].real)
    fmaxv=(fHz*np.sqrt(np.abs(num/den)))    
    fmax=sum(fmaxv[start:stop])/1e9/len(fmaxv[start:stop])
    return fmax, fmaxv

def Y4PtoRgg(Y, fHz, fstart=10, fstop=100e9):
    start, stop = findHz(fHz, fstart, fstop)
    Xv=(np.real(1/Y[:,0,0]))
    X=sum(Xv[start:stop])/len(Xv[start:stop])
    return X, Xv

def Y4PtoRw(Y, fHz, fstart=10, fstop=1000e9, Port=3):
    start, stop = findHz(fHz, fstart, fstop)
    ZZ=YtoZ4P(Y)
    if Port==4:
       Xv=(np.real(ZZ[:,3,3]-ZZ[:,2,3]))
    else: # BG sulla porta 3
       Xv=(np.real(ZZ[:,2,2]-ZZ[:,2,3]))
    X=sum(Xv[start:stop])/len(Xv[start:stop])
    return X, Xv

def Y4PurgeRgRw(Y, Rgv, Rwv, Port=3):
    ZZ=YtoZ4P(Y)
    ZZ[:,0,0]=ZZ[:,0,0]-Rgv
    if Port==4:
        ZZ[:,3,3]=ZZ[:,3,3]-Rwv
    else:
        ZZ[:,2,2]=ZZ[:,2,2]-Rwv
    YY=ZtoY4P(ZZ)
    return YY

def Y4PtoCgg(Y, fHz, fstart=10, fstop=100e9):
    start, stop = findHz(fHz, fstart, fstop)
    Cggv=-(1/np.imag(1/Y[:,0,0]))/(fHz*2*np.pi)
    Cgg=sum(Cggv[start:stop])/len(Cggv[start:stop])
    return Cgg, Cggv

def Y4PtoCgd(Y, fHz, fstart=10, fstop=100e9):
    start, stop = findHz(fHz, fstart, fstop)
    Xv=(1/np.imag(1/Y[:,0,1]))/(fHz*2*np.pi)
    X=sum(Xv[start:stop])/len(Xv[start:stop])
    return X, Xv

def Y4PtoCgs(Y, fHz, fstart=10, fstop=100e9,Port=3):
    start, stop = findHz(fHz, fstart, fstop)
    if Port==4:
        Xv=(1/np.imag(1/Y[:,0,2]))/(fHz*2*np.pi)
    else:
        Xv=(1/np.imag(1/Y[:,0,3]))/(fHz*2*np.pi)
    X=sum(Xv[start:stop])/len(Xv[start:stop])
    return X, Xv

def Y4PtoCgb(Y, fHz, fstart=10, fstop=100e9,Port=3):
    start, stop = findHz(fHz, fstart, fstop)
    if Port==4:
        Xv=(1/np.imag(1/Y[:,0,3]))/(fHz*2*np.pi)
    else:
        Xv=(1/np.imag(1/Y[:,0,2]))/(fHz*2*np.pi)
    X=sum(Xv[start:stop])/len(Xv[start:stop])
    return X, Xv

def Y4PtoCbg(Y, fHz, fstart=10, fstop=100e9,Port=3):
    start, stop = findHz(fHz, fstart, fstop)
    if Port==4:
        Xv=(1/np.imag(1/Y[:,3,0]))/(fHz*2*np.pi)
    else:
        Xv=(1/np.imag(1/Y[:,2,0]))/(fHz*2*np.pi)
    X=sum(Xv[start:stop])/len(Xv[start:stop])
    return X, Xv

def Y4PtoCbs(Y, fHz, fstart=10, fstop=100e9,Port=3):
    start, stop = findHz(fHz, fstart, fstop)
    if Port==4:
        Xv=(1/np.imag(1/Y[:,3,2]))/(fHz*2*np.pi)
    else:
        Xv=(1/np.imag(1/Y[:,2,3]))/(fHz*2*np.pi)
    X=sum(Xv[start:stop])/len(Xv[start:stop])
    return X, Xv

def Y4PtoCbd(Y, fHz, fstart=10, fstop=100e9,Port=3):
    start, stop = findHz(fHz, fstart, fstop)
    if Port==4:
        Xv=(1/np.imag(1/Y[:,3,1]))/(fHz*2*np.pi)
    else:
        Xv=(1/np.imag(1/Y[:,2,1]))/(fHz*2*np.pi)
    X=sum(Xv[start:stop])/len(Xv[start:stop])
    return X, Xv

def Y4PtoCbbiso(Y, fHz, fstart=10, fstop=100e9,Port=3):
    start, stop = findHz(fHz, fstart, fstop)
    if Port==4:
        Xv=-(1/np.imag(1/Y[:,3,3]))/(fHz*2*np.pi)
    else:
        Xv=-(1/np.imag(1/Y[:,2,2]))/(fHz*2*np.pi)
    X=sum(Xv[start:stop])/len(Xv[start:stop])
    return X, Xv

def Y4PtoGm(Y, fHz, fstart=10, fstop=1000e9):
    start, stop = findHz(fHz, fstart, fstop)
    Xv=(np.real(Y[:,1,0]-Y[:,0,1]))
    X=sum(Xv[start:stop])/len(Xv[start:stop])
    return X, Xv

def Y4PtoGmb(Y, fHz, fstart=10, fstop=1000e9, Port=3):
    start, stop = findHz(fHz, fstart, fstop)
    if Port==4:
        Xv=(np.real(Y[:,1,3]-Y[:,3,1]))
    else:
        Xv=(np.real(Y[:,1,2]-Y[:,2,1]))
    X=sum(Xv[start:stop])/len(Xv[start:stop])
    return X, Xv

def Y4PtoGds(Y, fHz, fstart=10, fstop=100e9, Port=3):
    start, stop = findHz(fHz, fstart, fstop)
    if Port==4:
        Xv=(np.real(-Y[:,2,1]))
    else: # BG su porta 3, lo standard
        Xv=(np.real(-Y[:,3,1]))
    X=sum(Xv[start:stop])/len(Xv[start:stop])
    return X, Xv

def Y4PtoRjd(Y, fHz, fstart=10, fstop=100e9, Port=3):
    start, stop = findHz(fHz, fstart, fstop)
    if Port==4:
        Xv=-(np.real(1/Y[:,3,1]))
    else: # BG su porta 3, lo standard
        Xv=-(np.real(1/Y[:,2,1]))
    X=sum(Xv[start:stop])/len(Xv[start:stop])
    return X, Xv

def Y4PtoRjs(Y, fHz, fstart=10, fstop=100e9, Port=3):
    start, stop = findHz(fHz, fstart, fstop)
    if Port==4:
        Xv=-(np.real(1/Y[:,3,2]))
    else: # BG su porta 3, lo standard
        Xv=-(np.real(1/Y[:,2,3]))
    X=sum(Xv[start:stop])/len(Xv[start:stop])
    return X, Xv

def Y4PtoCjd(Y, fHz, fstart=10, fstop=100e9,Port=3):
    start, stop = findHz(fHz, fstart, fstop)
    if Port==4:
        Xv=(1/np.imag(1/Y[:,3,1]))/(fHz*2*np.pi)
    else:
        Xv=(1/np.imag(1/Y[:,2,1]))/(fHz*2*np.pi)
    X=sum(Xv[start:stop])/len(Xv[start:stop])
    return X, Xv

def Y4PtoCjs(Y, fHz, fstart=10, fstop=100e9,Port=3):
    start, stop = findHz(fHz, fstart, fstop)
    if Port==4:
        Xv=(1/np.imag(1/Y[:,3,2]))/(fHz*2*np.pi)
    else:
        Xv=(1/np.imag(1/Y[:,2,3]))/(fHz*2*np.pi)
    X=sum(Xv[start:stop])/len(Xv[start:stop])
    return X, Xv

def Y4PtoRiso(Y, fHz, fstart=10, fstop=100e9,Port=3):
    start, stop = findHz(fHz, fstart, fstop)
    if Port==4:
        Yiso=Y[:,3,3]+Y[:,3,0]+Y[:,3,1]+Y[:,3,2]
    else:
        Yiso=Y[:,2,2]+Y[:,2,0]+Y[:,2,1]+Y[:,2,3]
    Xv=np.real(1/Yiso)/(fHz*2*np.pi)
    X=sum(Xv[start:stop])/len(Xv[start:stop])
    return X, Xv

def Y4PtoCiso(Y, fHz, fstart=10, fstop=100e9,Port=3):
    start, stop = findHz(fHz, fstart, fstop)
    if Port==4:
        Yiso=Y[:,3,3]+Y[:,3,0]+Y[:,3,1]+Y[:,3,2]
    else:
        Yiso=Y[:,2,2]+Y[:,2,0]+Y[:,2,1]+Y[:,2,3]
    Xv=-(1/np.imag(1/Yiso))/(fHz*2*np.pi)
    X=sum(Xv[start:stop])/len(Xv[start:stop])
    return X, Xv

def Y4PtoRdx(Y, fHz, fstart=10, fstop=100e9):
    start, stop = findHz(fHz, fstart, fstop)
    Ydd=Y[:,1,1]+Y[:,1,0]+Y[:,1,2]+Y[:,1,3]
    Xv=np.real(1/Ydd)/(fHz*2*np.pi)
    X=sum(Xv[start:stop])/len(Xv[start:stop])
    return X, Xv

def Y4PtoRgx(Y, fHz, fstart=10, fstop=100e9):
    start, stop = findHz(fHz, fstart, fstop)
    Ydd=Y[:,0,0]+Y[:,0,1]+Y[:,0,2]+Y[:,0,3]
    Xv=np.real(1/Ydd)/(fHz*2*np.pi)
    X=sum(Xv[start:stop])/len(Xv[start:stop])
    return X, Xv

def Y4PtoRsx(Y, fHz, fstart=10, fstop=100e9, Port=3):
    start, stop = findHz(fHz, fstart, fstop)
    if Port==4:
        Yss=Y[:,2,2]+Y[:,2,0]+Y[:,2,1]+Y[:,2,3]
    else:
        Yss=Y[:,3,3]+Y[:,3,0]+Y[:,3,1]+Y[:,3,2]
    Xv=np.real(1/Yss)/(fHz*2*np.pi)
    X=sum(Xv[start:stop])/len(Xv[start:stop])
    return X, Xv

def Y4PtoRbx(Y, fHz, fstart=10, fstop=100e9, Port=3):
    start, stop = findHz(fHz, fstart, fstop)
    if Port==4:
        Yss=Y[:,3,3]+Y[:,3,0]+Y[:,3,1]+Y[:,3,2]
    else:
        Yss=Y[:,2,2]+Y[:,2,0]+Y[:,2,1]+Y[:,2,3]
    Xv=np.real(1/Yss)/(fHz*2*np.pi)
    X=sum(Xv[start:stop])/len(Xv[start:stop])
    return X, Xv


#def Y4Ptoftbg3_old(Y, fHz, ft_start=0, ft_stop=-1):
#    # suppone che BG sia Port 3, drain port 2 
#    ## THIS IS OLD AND WRONG!
#    num = Y[:,1,2]*Y[:,0,0]-Y[:,1,0]*Y[:,0,2]
#    den = Y[:,2,2]*Y[:,0,0]-Y[:,0,3]*Y[:,0,3]
#    h23=num/den
#    ftbgv=(fHz*np.abs(h23))
#    ftbg=sum(ftbgv[ft_start:ft_stop])/1e9/len(ftbgv[ft_start:ft_stop])
#    return ftbg, ftbgv

#def Y4Ptoftbg3(Y, fHz, ft_start=0, ft_stop=-1):
#    # suppone che BG sia Port 3, drain port 2 
#    ## THIS IS CORRECT AS FAR AS I UNDERSTAND
#    num = Y[:,1,2]*Y[:,0,0]-Y[:,1,0]*Y[:,0,2]
#    den = Y[:,2,2]*Y[:,0,0]-Y[:,0,3]*Y[:,0,3]
#    h23=num/den
#    ftbgv=(fHz*np.abs(h23))
#    ftbg=sum(ftbgv[ft_start:ft_stop])/1e9/len(ftbgv[ft_start:ft_stop])
#    return ftbg, ftbgv

#def Y4Ptoft_BD(Y, fHz, ft_start=0, ft_stop=-1):
#    # suppone che BG sia Port 3, drain port 2
#    ## THIS IS CORRECT AS FAR AS I UNDERSTAND
#    h21=Y[:,1,0]/Y[:,0,0]
#    ftv=(fHz*np.abs(h21))
#    ft=sum(ftv[ft_start:ft_stop])/1e9/len(ftv[ft_start:ft_stop])
#    return ft, ftv

#def Y4Ptoft(Y, fHz, ft_start=0, ft_stop=-1):
#    # suppone che BG sia Port 3, drain port 2
#    ## THIS IS CORRECT AS FAR AS I UNDERSTAND 
#    num = Y[:,1,3]*Y[:,3,0]-Y[:,1,0]*Y[:,3,3]
#    den = Y[:,0,3]*Y[:,3,0]-Y[:,0,0]*Y[:,3,3]
#    h21=num/den
#    ftv=(fHz*np.abs(h21))
#    ft=sum(ftv[ft_start:ft_stop])/1e9/len(ftv[ft_start:ft_stop])
#    return ft, ftv

#def Y4Ptoftbg4(Y, fHz, ft_start=0, ft_stop=-1):
#    # suppone che fg sia P 1, drain p 2
#    num = Y[:,1,2]*Y[:,0,0]-Y[:,1,0]*Y[:,0,2]
#    den = Y[:,2,2]*Y[:,0,0]-Y[:,0,3]*Y[:,0,3]
#    h23=num/den
##    ftbgv=(fHz*np.abs(h23))
#    ftbg=sum(ftbgv[ft_start:ft_stop])/1e9/len(ftbgv[ft_start:ft_stop])
#    return ftbg, ftbgv



def S4PtoH2P(S):
    (l,x,y)=np.shape(S)
    H=np.empty((l,2,2),dtype=complex128)
    Z0=50
    for i in range(len(S)):
        s11=S[i,0,0]
        s12=S[i,0,1]
        s21=S[i,1,0]
        s22=S[i,1,1]
        den=(1-s11)*(Z0+s22*Z0)+s12*s21*Z0
        h11=((Z0+s11*Z0)*(Z0+s22*Z0)-s12*s21*Z0*Z0)/den
        h12=2*s12*Z0/den
        h21=-2*s21*Z0/den
        h22=(1-s11)*(1-s22)-s12*s21/den
        H[i,0,0]=h11
        H[i,0,1]=h12
        H[i,1,0]=h21
        H[i,1,1]=h22
    return H


################################################################################
### de-embedding raw S-parameters

def SdeembOS(Sdut, Sopen, Sshort):
    try:
       (b, f, x, y)=np.shape(Sdut)
    except ValueError as e:
       b=0
       (f, x, y)=np.shape(Sdut)
    if b>0:
      Sdev=np.empty((b, f, x, y), dtype=np.complex_)
    else:    
      Sdev=np.empty((f, x, y), dtype=np.complex_)

    Yopen2   =StoY(Sopen)  
    Yshort2  =StoY(Sshort) 
    Yshort2  = Yshort2 - Yopen2
    Zshort2  = YtoZ_2P(Yshort2)
    if b>0:
        for i in range(b):
            Ydev     =StoY(Sdut[i])
            Ydev     = Ydev - Yopen2
            Zdev     = YtoZ_2P(Ydev)
            Zdev     = Zdev - Zshort2
            Sdev[i]  = YtoS(ZtoY_2P(Zdev))
    else:
        Ydev  =StoY(Sdut)
        Ydev  = Ydev - Yopen2
        Zdev  = YtoZ_2P(Ydev)
        Zdev  = Zdev - Zshort2
        Sdev  = YtoS(ZtoY_2P(Zdev))

    return Sdev

################################################################################
### de-embedding network objects

def SSdeemb(op,sh): 
# demebedding dell'OPEN dallo short
    Ysh = sh.y.copy()
    Yop = op.y.copy()
    Ysh1=Ysh-Yop
    Ssh1=rf.y2s(Ysh1, z0=50)
    deembedded=rf.Network(f=op.frequency.f*1e-9,s=Ssh1,z0=50)
    return deembedded

def OSdeemb(dut,loc_open,loc_short):
    Yopen2   =loc_open.y.copy() 
    Yshort2  =loc_short.y.copy()
    Ydev     =dut.y.copy()
   
    # OPEN2 CORRECTION
    Yshort2  = Yshort2 - Yopen2 
    Ydev     = Ydev    - Yopen2
   
    # SHORT2 CORRECTION
    Zshort2  = rf.y2z(Yshort2)
    Zdev     = rf.y2z(Ydev)

    Zdev     = Zdev - Zshort2

    Sdev     = rf.z2s(Zdev, z0=[50.0+0.0j, 50.0+0.0j])
    deembedded=rf.Network(f=dut.frequency.f*1e-9,s=Sdev,z0=dut.z0)
    return deembedded

def OSdeemb4P_S(dut,pad_open,gen_open,gen_short,loc_open,loc_short):
    Ypad     =StoY4P(pad_open)
    Yshort   =StoY4P(gen_short) 
    Yopen    =StoY4P(gen_open)
    Yopen2   =StoY4P(loc_open) 
    Yshort2  =StoY4P(loc_short)
    Ydev     =StoY4P(dut)

    # PAD CORRECTION
    Yshort   = Yshort   - Ypad
    Yopen    = Yopen    - Ypad
    Yopen2   = Yopen2   - Ypad
    Yshort2  = Yshort2  - Ypad
    Ydev     = Ydev     - Ypad

    # OPEN CORRECTION
    Yshort   = Yshort   - Yopen
    Yopen2   = Yopen2   - Yopen
    Yshort2  = Yshort2  - Yopen
    Ydev     = Ydev     - Yopen

    # SHORT CORRECTION
    Zshort    = YtoZ4P(Yshort)
    Zopen2    = YtoZ4P(Yopen2)
    Zshort2   = YtoZ4P(Yshort2)
    Zdev      = YtoZ4P(Ydev)
    # Zshort = np.linalg.inv(Yshort) # alterative formulation
    # Zopen2 = np.linalg.inv(Yopen2)
    # Zshort2= np.linalg.inv(Yshort2)
    # Zdut   = np.linalg.inv(Ydut)
    
    Zopen2  = Zopen2  - Zshort
    Zshort2 = Zshort2 - Zshort
    Zdev    = Zdev    - Zshort

    # OPEN2 CORRECTION
    Yopen2   = ZtoY4P(Zopen2)
    Yshort2  = ZtoY4P(Zshort2)
    Ydev     = ZtoY4P(Zdev)

    Yshort2  = Yshort2 - Yopen2 
    Ydev     = Ydev    - Yopen2


    # SHORT2 CORRECTION
    Zshort2  = YtoZ4P(Yshort2)
    Zdev     = YtoZ4P(Ydev)

    Zdev     = Zdev - Zshort2

    Sdev     = ZtoS4P(Zdev)

#    deembedded=rf.Network(f=dut_es.frequency.f,s=Sdev,z0=dut_es.z0)
#    return deembedded
    return Sdev

def OSdeemb2P(dut,pad_open,gen_open,gen_short,loc_open,loc_short):
    Ypad     =pad_open.y.copy()
    Yshort   =gen_short.y.copy() 
    Yopen    =gen_open.y.copy()
    Yopen2   =loc_open.y.copy() 
    Yshort2  =loc_short.y.copy()
    Ydev     =dut.y.copy()

    # PAD CORRECTION
    Yshort   = Yshort   - Ypad
    Yopen    = Yopen    - Ypad
    Yopen2   = Yopen2   - Ypad
    Yshort2  = Yshort2  - Ypad
    Ydev     = Ydev     - Ypad

    # OPEN CORRECTION
    Yshort   = Yshort   - Yopen
    Yopen2   = Yopen2   - Yopen
    Yshort2  = Yshort2  - Yopen
    Ydev     = Ydev     - Yopen

    # SHORT CORRECTION
    Zshort    = rf.y2z(Yshort)
    Zopen2    = rf.y2z(Yopen2)
    Zshort2   = rf.y2z(Yshort2)
    Zdev      = rf.y2z(Ydev)
    # Zshort = np.linalg.inv(Yshort) # alterative formulation
    # Zopen2 = np.linalg.inv(Yopen2)
    # Zshort2= np.linalg.inv(Yshort2)
    # Zdut   = np.linalg.inv(Ydut)
    
    Zopen2  = Zopen2  - Zshort
    Zshort2 = Zshort2 - Zshort
    Zdev    = Zdev    - Zshort

    # OPEN2 CORRECTION
    Yopen2   = rf.z2y(Zopen2)
    Yshort2  = rf.z2y(Zshort2)
    Ydev     = rf.z2y(Zdev)

    Yshort2  = Yshort2 - Yopen2 
    Ydev     = Ydev    - Yopen2


    # SHORT2 CORRECTION
    Zshort2  = rf.y2z(Yshort2)
    Zdev     = rf.y2z(Ydev)

    Zdev     = Zdev - Zshort2

    Sdev     = rf.z2s(Zdev, z0=[50.0+0.0j, 50.0+0.0j])

    deembedded=rf.Network(f=dut.frequency.f,s=Sdev,z0=dut.z0)
    return deembedded

def OSdeemb2P_1(ndut,nopen,nshort):
    Yopen   =nopen.y.copy() 
    Yshort  =nshort.y.copy()
    Ydev     =ndut.y.copy()
   
    # OPEN CORRECTION
    Yshort  = Yshort - Yopen 
    Ydev    = Ydev   - Yopen
   
    # SHORT CORRECTION
    Zshort  = rf.y2z(Yshort) # prova a sostituire con l'inversa della matrice
    Zdev    = rf.y2z(Ydev)   # prova a sostituire con l'inversa della matrice

    Zdev    = Zdev - Zshort

    Sdev    = rf.z2s(Zdev, z0=[50.0+0.0j, 50.0+0.0j])
    deembedded=rf.Network(f=ndut.frequency.f*1e-9,s=Sdev,z0=ndut.z0)
    return deembedded

def experimental_SOdeemb2P_1(ndut,nopen,nshort):
    Yopen   =nopen.y.copy() 
    Yshort  =nshort.y.copy()
    Ydev     =ndut.y.copy()
   
    # Y1 = Yopen[:,0,0]+Yopen[:,0,1]
    # Y2 = Yopen[:,1,1]+Yopen[:,0,1]

    #         |  Yo11+Yo12    0        |
    # Yo2=    |     0       Yo22+Yo12  |
    Yopen2 = np.zeros(Yopen.shape,dtype=np.complex128)
    Yopen2[:,0,0]=Yopen[:,0,0]+Yopen[:,0,1]
    Yopen2[:,1,1]=Yopen[:,1,1]+Yopen[:,0,1]

    Ydev   = Ydev   - Yopen2
    Yshort = Yshort - Yopen2
   
    Zshort  = rf.y2z(Yshort) 
    Zdev    = rf.y2z(Ydev)  
  
    Zb = Zdev - Zshort
    
    Z1=Zshort[:,0,0]-Zshort[:,0,1]
    Z2=Zshort[:,1,1]-Zshort[:,0,1]
    Y3 = (- Yopen[:,0,1])/(1+ Yopen[:,0,1]*(Z1+Z2))
    
    Ydummy=np.zeros(Yopen.shape , dtype=np.complex128)
    Ydummy[:,0,0]=Y3
    Ydummy[:,0,1]=-Y3
    Ydummy[:,1,0]=-Y3
    Ydummy[:,1,1]=Y3

    Yb=rf.z2y(Zb)
    
    Ydeemb=Yb-Ydummy

    Sdev    = rf.y2s(Ydeemb, z0=[50.0+0.0j, 50.0+0.0j])
    deembedded=rf.Network(f=ndut.frequency.f*1e-9,s=Sdev,z0=ndut.z0)
    return deembedded



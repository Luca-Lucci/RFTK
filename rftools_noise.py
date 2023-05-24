################################################################################
##### NOISE
import matplotlib.pyplot as plt
import numpy as np
import copy
import sys
import gzip
from rftools import *
#from Signatone import *

def CrAtoCrY(CrA, Y):
    CrY=np.empty(np.shape(CrA),dtype=np.complex128)
    T=np.zeros((2, 2), dtype=np.complex128)
    Tadj=np.zeros((2, 2), dtype=np.complex128)
    for i in range(len(CrA)):
        T[0,0] =-Y[i,0,0] 
        T[0,1] = complex(1.0,0.0)
        T[1,0] =-Y[i,1,0]
        T[1,1] = complex(0.0,0.0)
   
        Tadj[0,0] = np.conj(T[0,0])
        Tadj[0,1] = np.conj(T[1,0])
        Tadj[1,0] = np.conj(T[0,1])
        Tadj[1,1] = np.conj(T[1,1])
        CrY[i,:,:]=np.dot(np.dot(T,CrA[i,:,:]),Tadj)
    return CrY

def CrYtoCrA(CrY, A):
    CrA=np.empty(np.shape(CrY),dtype=np.complex128)
    T=np.zeros((2, 2), dtype=np.complex128)
    Tadj=np.zeros((2, 2), dtype=np.complex128)
    for i in range(len(CrA)):
        T[0,0] = complex(0.0,0.0)
        T[0,1] = A[i,0,1]
        T[1,0] = complex(1.0,0.0)
        T[1,1] = A[i,1,1]
   
        Tadj[0,0] = np.conj(T[0,0])
        Tadj[0,1] = np.conj(T[1,0])
        Tadj[1,0] = np.conj(T[0,1])
        Tadj[1,1] = np.conj(T[1,1])
        CrA[i,:,:]=np.dot(np.dot(T,CrY[i,:,:]),Tadj)
    return CrA

def CrYtoCrA_fromY(CrY, Y):
    CrA=np.empty(np.shape(CrY),dtype=np.complex128)
    T=np.zeros((2, 2), dtype=np.complex128)
    Tinv=np.zeros((2, 2), dtype=np.complex128)
    Tadj=np.zeros((2, 2), dtype=np.complex128)
    for i in range(len(CrA)): 
        T[0,0] =-Y[i,0,0] 
        T[0,1] = complex(1.0,0.0)
        T[1,0] =-Y[i,1,0]
        T[1,1] = complex(0.0,0.0)
        
        Tinv=np.linalg.inv(T)
        Tadj=np.transpose(np.conj(Tinv))
        
        CrA[i,:,:]=np.dot(np.dot(Tinv,CrY[i,:,:]),Tadj)
        
    return CrA

def CrYtoCrZ(CrY, Y):
    Z=YtoZ_2P(Y)
    CrZ=np.empty(np.shape(CrY),dtype=np.complex128)
    T=np.zeros((2, 2), dtype=np.complex128)
    Tinv=np.zeros((2, 2), dtype=np.complex128)
    Tadj=np.zeros((2, 2), dtype=np.complex128)
    for i in range(len(CrY)): 
        T[0,0] = Z[i,0,0] 
        T[0,1] = Z[i,0,1]
        T[1,0] = Z[i,1,0]
        T[1,1] = Z[i,1,1]
        #Tinv=np.linalg.inv(T)
        #Tadj=np.transpose(np.conj(T))
        Tadj[0,0] = np.conj(T[0,0])
        Tadj[0,1] = np.conj(T[1,0])
        Tadj[1,0] = np.conj(T[0,1])
        Tadj[1,1] = np.conj(T[1,1])
        
        CrZ[i,:,:]=np.dot(np.dot(T,CrY[i,:,:]),Tadj)
        
    return CrZ

def CrAtoCrZ(CrA, Y):
    Z=YtoZ_2P(Y)
    CrZ=np.empty(np.shape(CrA),dtype=np.complex128)
    T=np.zeros((2, 2), dtype=np.complex128)
    Tadj=np.zeros((2, 2), dtype=np.complex128)
    for i in range(len(CrA)):
        T[0,0] = complex(1.0,0.0)
        T[0,1] = -Z[i,0,0]
        T[1,0] = complex(0.0,0.0)
        T[1,1] = -Z[i,1,0]
   
        Tadj[0,0] = np.conj(T[0,0])
        Tadj[0,1] = np.conj(T[1,0])
        Tadj[1,0] = np.conj(T[0,1])
        Tadj[1,1] = np.conj(T[1,1])
        CrZ[i,:,:]=np.dot(np.dot(T,CrA[i,:,:]),Tadj)
    return CrZ

def CrZtoCrA(CrZ, A):
    CrA=np.empty(np.shape(CrZ),dtype=np.complex128)
    T=np.zeros((2, 2), dtype=np.complex128)
    Tadj=np.zeros((2, 2), dtype=np.complex128)
    for i in range(len(CrZ)):
        T[0,0] = complex(1.0,0.0)
        T[0,1] = -A[i,0,0]
        T[1,0] = complex(0.0,0.0)
        T[1,1] = -A[i,1,0]
   
        Tadj[0,0] = np.conj(T[0,0])
        Tadj[0,1] = np.conj(T[1,0])
        Tadj[1,0] = np.conj(T[0,1])
        Tadj[1,1] = np.conj(T[1,1])
        CrA[i,:,:]=np.dot(np.dot(T,CrZ[i,:,:]),Tadj)
    return CrA

def CrZtoCrY(CrZ, Y):
    CrY=np.empty(np.shape(CrZ),dtype=np.complex128)
    T=np.zeros((2, 2), dtype=np.complex128)
    Tadj=np.zeros((2, 2), dtype=np.complex128)
    for i in range(len(CrZ)):
        T[0,0] = Y[i,0,0]
        T[0,1] = Y[i,0,1]
        T[1,0] = Y[i,1,0]
        T[1,1] = Y[i,1,1]
   
        Tadj[0,0] = np.conj(T[0,0])
        Tadj[0,1] = np.conj(T[1,0])
        Tadj[1,0] = np.conj(T[0,1])
        Tadj[1,1] = np.conj(T[1,1])
        CrY[i,:,:]=np.dot(np.dot(T,CrZ[i,:,:]),Tadj)
    return CrY

def CrHtoCrA(CrH, A):
    CrA=np.empty(np.shape(CrH),dtype=np.complex128)
    T=np.zeros((2, 2), dtype=np.complex128)
    Tadj=np.zeros((2, 2), dtype=np.complex128)
    for i in range(len(CrH)):
        T[0,0] = complex(1.0,0.0)
        T[0,1] = A[i,0,1]
        T[1,0] = complex(0.0,0.0)
        T[1,1] = A[i,1,1]
   
        Tadj[0,0] = np.conj(T[0,0])
        Tadj[0,1] = np.conj(T[1,0])
        Tadj[1,0] = np.conj(T[0,1])
        Tadj[1,1] = np.conj(T[1,1])
        CrA[i,:,:]=np.dot(np.dot(T,CrH[i,:,:]),Tadj)
    return CrA


def HtoCrH(H, NFlin, Ys=complex(1.0/50.0,0)):
    Cr=np.empty(np.shape(H),dtype=np.complex128)
    
    for i in range(len(H)):
        h11=H[i,0,0]
        h12=H[i,0,1]
        h21=H[i,1,0]
        h22=H[i,1,1]
        nfl=NFlin[i]
        a22=(np.absolute(h21/(1 + Ys*h11))**2)*(np.real(Ys)*(nfl-1.0)-np.absolute(Ys)**2*np.real(h11))
        a11=np.real(h11)
        a12=complex(0.0)
        a21=complex(0.0)
        Cr[i,0,0]=a11
        Cr[i,0,1]=a12
        Cr[i,1,0]=a21
        Cr[i,1,1]=a22
    return Cr

def Fmin_from_CrA(CrA):
    Fmin=np.zeros(len(CrA),dtype=np.float)
    for i in range(len(CrA)):
        c11=CrA[i,0,0]
        c12=CrA[i,0,1]
        c21=CrA[i,1,0]
        c22=CrA[i,1,1]
        ret=1.0+2.0*np.real(np.real(c12)+np.sqrt(c11*c22-np.imag(c12)**2))
        Fmin[i]=ret
    return Fmin

def Yopt_from_CrA(CrA):
    Yopt=np.zeros(len(CrA), dtype=np.complex128)
    for i in range(len(CrA)):
       CrA11=CrA[i,0,0]
       CrA22=CrA[i,1,1]
       CrA12=CrA[i,0,1]
       Yopt[i]=(np.sqrt(CrA11*CrA22-np.imag(CrA12)**2)+complex(0.0,np.imag(CrA12)))/(CrA11)
    return Yopt

def Rn_from_CrA(CrA):
    return np.real(CrA[:,0,0])

def noise_from_CrA(CrA):
    Fmin=Fmin_from_CrA(CrA);
    NFmin=10*np.log10(Fmin);
    Yopt=Yopt_from_CrA(CrA);
    Rn  =Rn_from_CrA(CrA);
    return (NFmin, Fmin, Yopt, Rn)

def NF50_Hoverride(Sparam, Hparam, NFmeas):
#   The idea is to fix re(H11) to a fixed value (or the a values at Vds=0V)
#        Hdut[bias]=rftools.StoH(Sdut[bias])
#        Hpdut[bias]=rftools.StoH(Sdut[bias])       # Hp is Hprime
#        (fp, nx, ny)=np.shape(Hdut[bias])
#        for jj in range(fp):
#                 Hpdut[bias][jj,0,0]=np.complex( new_reH11, np.imag(Hdut[NFdev][bias][jj,0,0]) )
#        NFmin[bias], Fmin[bias], Yopt[bias], Rn[bias])= \
#            rftools.NF50_Hoverride(Sdut[bias], Hpdut[bias], NFdut[bias])
#   I think both H and S must be already de-embedded
    Adut=StoABCD(Sparam)
    #Hdut=StoH(Sparam)
    NFlin=10**(NFmeas/10.0)
    CrHdut=HtoCrH(Hparam,NFlin)
    CrAdut=CrHtoCrA(CrHdut,Adut)
    (NFmin, Fmin, Yopt, Rn) = noise_from_CrA(CrAdut)
    return (NFmin, Fmin, Yopt, Rn) 
 
def NF50(Sparam, NFmeas):
    # la matrice S deve esere gia stata de-embeddata
    Adut=StoABCD(Sparam)
    Hdut=StoH(Sparam)
    NFlin=10**(NFmeas/10.0)
    CrHdut=HtoCrH(Hdut,NFlin)
    CrAdut=CrHtoCrA(CrHdut,Adut)
    (NFmin, Fmin, Yopt, Rn) = noise_from_CrA(CrAdut)
    return (CrAdut, NFmin, Fmin, Yopt, Rn) 

def PucelNF50(Sdut, Sopen, Sshort, NFmeas, Hparam):
    # analogamente anche qui, la matrice S in Sdut deve gia essere stata de-embeddata!
    # qui si fa il de-embedding solo della contribuzione dei parassiti da NF
    Adut=StoABCD(Sdut)
    # Hdut=StoH(Sdut) # we can override this step!
    Ydut=StoY(Sdut)
    NFlin=10**(NFmeas/10.0)
    #CrHdut=HtoCrH(Hdut,NFlin)
    CrHdut=HtoCrH(Hparam,NFlin)
    CrAdut=CrHtoCrA(CrHdut,Adut)
    CrYdut=CrAtoCrY(CrAdut,Ydut)
    Yopen =StoY(Sopen)
    Yshort=StoY(Sshort)-Yopen
    Zshort=np.linalg.inv(Yshort)
    tmpYdut=StoY(Sdut)-Yopen
    Ydutd=np.linalg.inv(np.linalg.inv(tmpYdut)-Zshort)
    Adutd=YtoA_2P(Ydutd)

    Zc=(Zshort[:,0,1]+Zshort[:,1,0])/2.0
    Za=Zshort[:,0,0]-Zc
    Zb=Zshort[:,1,1]-Zc

    Yc=-(Yopen[:,0,1]+Yopen[:,1,0])/2.0
    Ya=Yopen[:,0,0]-Yc
    Yb=Yopen[:,1,1]-Yc

    Yp=np.zeros((len(Sopen),4,4),dtype=np.complex128)
    dZ=Za*Zb+Zb*Zc+Zc*Za
    Yp[:,0,0]=(Zb+Zc)/dZ+Ya+Yc
    Yp[:,0,1]=(-Zc)/dZ-Yc
    Yp[:,0,2]=-(Zb+Zc)/dZ
    Yp[:,0,3]=(Zc)/dZ

    Yp[:,1,0]=(-Zc)/dZ-Yc
    Yp[:,1,1]=(Za+Zc)/dZ+Yb+Yc
    Yp[:,1,2]=(Zc)/dZ
    Yp[:,1,3]=-(Za+Zc)/dZ

    Yp[:,2,0]=-(Zb+Zc)/dZ
    Yp[:,2,1]=(Zc)/dZ
    Yp[:,2,2]=(Zb+Zc)/dZ
    Yp[:,2,3]=-(Zc)/dZ

    Yp[:,3,0]=(Zc)/dZ
    Yp[:,3,1]=-(Za+Zc)/dZ
    Yp[:,3,2]=-(Zc)/dZ
    Yp[:,3,3]=(Za+Zc)/dZ

    Yee=copy.deepcopy(Yp[:,:2,:2])
    Yii=copy.deepcopy(Yp[:,2:,2:])
    Yei=copy.deepcopy(Yp[:,:2,2:])
    Yie=copy.deepcopy(Yp[:,2:,:2])

    D   =np.zeros(np.shape(Yei),dtype=np.complex128)
    Dinv=np.zeros(np.shape(Yei),dtype=np.complex128)
    Dadj=np.zeros(np.shape(Yei),dtype=np.complex128)
    for i in range(len(Yei)):
        D[i]=-np.dot(Yei[i],np.linalg.inv(Ydutd[i]+Yii[i]))
        Dinv[i]=np.linalg.inv(D[i])
        Dadj[i]=np.conj(Dinv[i])
        tmp=Dadj[i,0,1]
        Dadj[i,0,1]=Dadj[i,1,0]
        Dadj[i,1,0]=tmp

    Cee = np.real(Yee)
    Cei = np.real(Yei)
    Cii = np.real(Yii)
    Cie = np.real(Yie)

    CrY=np.zeros(np.shape(Yei),dtype=np.complex128)
    for i in range(len(CrY)):
        CrY[i]=np.dot(np.dot(Dinv[i],CrYdut[i]-Cee[i]),Dadj[i])-\
               np.dot(Cie[i],Dadj[i])-\
               np.dot(Dinv[i],Cei[i])-\
               Cii[i];

    CrA =CrYtoCrA(CrY,Adutd)
    (NFmin, Fmin, Yopt, Rn) = noise_from_CrA(CrA)

    # calcolo di Zs, non sono molto convinto
    # passo la matrice Yp ai morsetti 1 e 3
    Yp13=np.zeros((len(Yp),2,2), dtype=np.complex128)
    Yp13[:,0,0]=copy.deepcopy(Yp[:,0,0])
    Yp13[:,0,1]=copy.deepcopy(Yp[:,0,2])
    Yp13[:,1,0]=copy.deepcopy(Yp[:,2,0])
    Yp13[:,1,1]=copy.deepcopy(Yp[:,2,2])
 
    Sp13=YtoS(Yp13)
    Gamma=Sp13[:,1,1]
    Zs=50.0*(1+Gamma)/(1-Gamma)
    #return ( NFmin, Fmin, Yopt, Rn, Zs ) # original version before 27-08-20    
    return (CrY, CrA, NFmin, Fmin, Yopt, Rn, Zs)

def getNoisefromS2P (file):
    f=file.readlines()
    Sstop=np.size(f)
    print (Sstop)
    test=False
    for counter, value in enumerate(f): # This loop is to get the frequency range
        vect=value.strip().split(' ')
        tupl=vect
        #print (tupl)
        while '' in tupl:
            tupl.remove('') 
        tupl2=[]
        for i in tupl:
            tupl2.append(i)
        #print (tupl2)
        mode=str.isnumeric((tupl2[0][0]))
        if mode==True and test==False:
            Sstart=counter
            print (Sstart)
            test=True
    numOfFreqPoints=Sstop-Sstart
    NFmin=np.zeros((numOfFreqPoints),dtype=np.float64)
    Rn=np.zeros((numOfFreqPoints),dtype=np.float64)
    Sopt=np.zeros((numOfFreqPoints),dtype=np.complex128)
    NF=np.zeros((numOfFreqPoints),dtype=np.float64)
    Te=np.zeros((numOfFreqPoints),dtype=np.float64)

    for counter, value in enumerate(f):
        if counter>=Sstart  and counter<=Sstop:
            vect=value.strip().split(' ')
            tupl=vect
            while '' in tupl:
                tupl.remove('')  
            tupl2=[]
            for i in tupl:
                tupl2.append(float(i))
            NFmin_tmp=tupl2[1]
            Rn_tmp=tupl2[2]
            Sopt_tmp=tupl2[3]+tupl2[4]*1j
            nf_tmp=tupl2[6]
            te_tmp=tupl2[8]
            freqpoint=counter-Sstart
            NFmin[freqpoint]=NFmin_tmp
            Rn[freqpoint]=Rn_tmp
            Sopt[freqpoint]=Sopt_tmp
            NF[freqpoint]=nf_tmp
            Te[freqpoint]=te_tmp

    return NFmin, Rn, Sopt, NF, Te


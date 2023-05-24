# -*- coding: utf-8 -*-
# v0.2 : - modify and rename parseWincalSOLTraw to handle both 
#          separation and on substrate OPENi
#        - remove ports parameter in parseWincal16ETto12
import os    as os
import re    as re
import numpy as np
import scipy as sp
import pickle 
import gzip
import subprocess
import copy
import matplotlib.pyplot as plt

###### dedicated
import rftools as rft
import spicetools as spt
import smithplot
from smithplot import SmithAxes
from common import *


whatlist=[ 'ED', 'ES', 'ERT', 'EDr', 'ESr', 'ERTr','EL', 'ETT','ELr', 'ETTr', 'EX', 'EXr','Gamma12','Gamma21']   
whatname={ 'ED':'Fwd directivity', 
           'ES':'Fwd source match',
           'ERT':'Fwd reflection tracking',
           'EDr':'Rev directivity',
           'ESr':'Rev source match',
           'ERTr':'Rev reflection tracking',
           'EL':'Fwd load match', 
           'ETT':'Fwd transmission tracking',
           'ELr':'Rev load match', 
           'ETTr':'Rev transmission tracking',
           'EX':'Fwd isolation',
           'EXr':'Rev isolation',
           'Gamma12':'switch term 12',
           'Gamma21':'switch term 21'}


def parseWincal12ET(loc="./",ports=(1,2),Verbose=False, ReadFrequency=False):
    """ 
    v0.1   parse Wincal's error terms if saved in the 12-error term formalism
    ex: 
         ErrTrms=parseWincal12ET(f"{datadir}/102_ErrorTerms/",ports=(3,4), Verbose=False, ReadFrequency=True)
    """
    err =dict()
    
    pin =f"{ports[0]}{ports[0]}"
    pout=f"{ports[1]}{ports[1]}"
    rev =f"{ports[0]}{ports[1]}"    
    fw  =f"{ports[1]}{ports[0]}"
    if Verbose:
        print(f"parsing 12 term errors in {loc} using pin={pin}, pout={pout}, rev={rev} and fw={fw}")
    if ReadFrequency:
       err['F'],err['ED']  =spt.SnPparse(f"{loc}/ErrorTerm_{pin}_EDir.s1p",  Verbose=Verbose, ReadFrequency=ReadFrequency) # forward directivity
    else:
       err['ED']  =spt.SnPparse(f"{loc}/ErrorTerm_{pin}_EDir.s1p",  Verbose=Verbose) # forward directivity
    err['ES']  =spt.SnPparse(f"{loc}/ErrorTerm_{pin}_ESrm.s1p",  Verbose=Verbose) #  source match
    err['ERT'] =spt.SnPparse(f"{loc}/ErrorTerm_{pin}_ERft.s1p",  Verbose=Verbose) # Reflection tracking
    err['EL']  =spt.SnPparse(f"{loc}/ErrorTerm_{fw}_ELdm.s1p",   Verbose=Verbose) # Load match 

    err['EDr'] =spt.SnPparse(f"{loc}/ErrorTerm_{pout}_EDir.s1p", Verbose=Verbose) # reverse directivity
    err['ESr'] =spt.SnPparse(f"{loc}/ErrorTerm_{pout}_ESrm.s1p", Verbose=Verbose) # reverse source match
    err['ERTr']=spt.SnPparse(f"{loc}/ErrorTerm_{pout}_ERft.s1p", Verbose=Verbose) # reverse Reflection tracking
    err['ELr'] =spt.SnPparse(f"{loc}/ErrorTerm_{rev}_ELdm.s1p",  Verbose=Verbose) # reverse Load match 
    
    err['ETT'] =spt.SnPparse(f"{loc}/ErrorTerm_{fw}_ETrt.s1p",   Verbose=Verbose) # reverse Transmission tracking
    err['ETTr']=spt.SnPparse(f"{loc}/ErrorTerm_{rev}_ETrt.s1p",  Verbose=Verbose) # forward Transmission tracking

    err['EX']  =spt.SnPparse(f"{loc}/ErrorTerm_{fw}_EXtlk.s1p",  Verbose=Verbose) # Crosstalk 
    err['EXr'] =spt.SnPparse(f"{loc}/ErrorTerm_{rev}_EXtlk.s1p", Verbose=Verbose) # Crosstalk
    return err

def parseWincal16ETto12(loc='./', Verbose=False, ReadFrequency=False):
    """
    v0.2  no need to have ports=(1,2) parameters. deleted
    v0.1  parse 16 ET and CONVERTS to 12 error terms
    """
    err=dict()
    if ReadFrequency:
        f,ErrCo = spt.SnPparse(f"{loc}/WC4Et16Term.S4P", ReadFrequency=True)
    else:
        ErrCo = spt.SnPparse(f"{loc}/WC4Et16Term.S4P", ReadFrequency=False)
    SwT = spt.SnPparse(f"{loc}/WC4Et16TermSwitch.S2P")

    Gamma12=SwT[:,0,1]
    Gamma21=SwT[:,1,0]
    print("Debug: how to handle this Switch terms!?!?")
    if ReadFrequency:
        err['F']=f
    err['Gamma12']=Gamma12
    err['Gamma21']=Gamma21

    err['ED']=ErrCo[:,0,0]
    err['EDr']=ErrCo[:,1,1]
    err['ES']=ErrCo[:,2,2]
    err['ESr']=ErrCo[:,3,3]
    err['ERT']=ErrCo[:,2,0]*ErrCo[:,0,2]
    err['ERTr']=ErrCo[:,3,1]*ErrCo[:,1,3]

    err['EX']=ErrCo[:,1,0]/(1-Gamma21*ErrCo[:,1,1])
    err['EXr']=ErrCo[:,0,1]/(1-Gamma12*ErrCo[:,0,0])

    err['ETT']=ErrCo[:,2,0]*ErrCo[:,1,3]/(1-Gamma21*ErrCo[:,1,1])
    err['ETTr']=ErrCo[:,3,1]*ErrCo[:,0,2]/(1-Gamma12*ErrCo[:,0,0])

    err['EL'] =ErrCo[:,3,3]+ErrCo[:,1,3]*ErrCo[:,3,1]*Gamma21/(1-Gamma21*ErrCo[:,1,1])
    err['ELr']=ErrCo[:,2,2]+ErrCo[:,0,2]*ErrCo[:,2,0]*Gamma12/(1-Gamma12*ErrCo[:,0,0])

    return err

ErrorTermsNames=[ 'ED', 'ES', 'ERT', 'EL', 'EDr', 'ESr', 'ERTr', 'ELr', 'ETT', 'ETTr', 'EX', 'EXr']

def parseWincalSOLTraw(loc="./",ports=(1,2), SubCal='101-190C', Verbose=False, ReadFrequency=True):
    """
    OBSOLETE - pls use parse_wincal_raw instead
    v 0.1
        parse the RAW SOLT with or without Switch terms data as saved from Wincal
    ex:
        iss15=parseWincalSOLTraw(datadir+"/15_SOLT_PNAX_whileCapella/", ports=(3,4), ReadFrequency=True)
    """
    raw = parse_winca_raw(loc=loc,ports=ports, SubCal=SubCal, Verbose=Verbose, ReadFrequency=False)
    print(" OBSOLETE - pls use parse_wincal_raw instead ")
    return raw

def parse_wincal_raw(loc="./",ports=(1,2), SubCal='101-190C', Verbose=False, ReadFrequency=True):
    """
    v 0.2 fix ReadFrequency bug. Add features to handle separation and on_subs opens
    v 0.1 parse the RAW SOLT with or without Switch terms data as saved from Wincal
    ex:
        iss=parseWincalSOLTraw(datadir+"/15_SOLT_PNAX_whileCapella/", ports=(3,4), ReadFrequency=True)
    """
    from os.path import exists
    raw =dict()
    
    pin =f"{ports[0]}"
    pout=f"{ports[1]}"
    rev =f"{ports[0]}{ports[1]}"
    if Verbose:
        print(f"parsing raw ISS SOLT data (no switcht terms) in {loc} using pin={pin}, pout={pout}, and rev={rev}")
    # separation OPENS
    open1=False
    open2=False
    fopenP1=f"{loc}/Separate^S_PARA^{pin}.s1p"
    fopenP2=f"{loc}/Separate^S_PARA^{pout}.s1p"
    if exists(fopenP1):
        print("Separation Open")
        open1=True
        raw['OpenP1'] =spt.SnPparse(fopenP1,  Verbose=Verbose)
    if exists(fopenP2): 
        open2=True
        raw['OpenP2'] =spt.SnPparse(fopenP2,  Verbose=Verbose)
    # on sub open
    fopenP1=f"{loc}/{SubCal} Open (On Sub)^S_PARA^{pin}.s1p"
    fopenP2=f"{loc}/{SubCal} Open (On Sub)^S_PARA^{pout}.s1p"
    if exists(fopenP1):
        if open1:
            print("WARNING, Both separation and on sub open found for Port1. Latter is saved")
        print(f"On substrate Open")
        raw['OpenP1'] =spt.SnPparse(fopenP1,  Verbose=Verbose)
    if exists(fopenP2): 
        if open1:
            print("WARNING, Both separation and on sub open found for Port2. Latter is saved")
        raw['OpenP2'] =spt.SnPparse(fopenP2,  Verbose=Verbose)
    # short, may not be there in case of LRM
    fshort1=f"{loc}/{SubCal} Short^S_PARA^{pin}.s1p"
    fshort2=f"{loc}/{SubCal} Short^S_PARA^{pout}.s1p"
    if exists(fshort1): raw['ShortP1']=spt.SnPparse(fshort1, Verbose=Verbose)
    if exists(fshort2): raw['ShortP2']=spt.SnPparse(fshort2, Verbose=Verbose)
    # load
    raw['LoadP1']=spt.SnPparse(f"{loc}/{SubCal} Load^S_PARA^{pin}.s1p", Verbose=Verbose)
    raw['LoadP2']=spt.SnPparse(f"{loc}/{SubCal} Load^S_PARA^{pout}.s1p", Verbose=Verbose) 

    # thru and frequency sweep
    if ReadFrequency:
        raw['F'], raw['Thru']  =spt.SnPparse(f"{loc}/{SubCal} Thru^S_PARA^{rev}.s2p", Verbose=Verbose, ReadFrequency=True) 
    else:
        raw['Thru']  =spt.SnPparse(f"{loc}/{SubCal} Thru^S_PARA^{rev}.s2p", Verbose=Verbose)      

    # switch terms
    switchTerms=f"{loc}/{SubCal} Thru^SWITCH_GAMMA^{rev}.s2p"
    if exists(switchTerms):
        print(f"SwitchTerms found")
        SwT=spt.SnPparse(f"{loc}/{SubCal} Thru^S_PARA^{rev}.s2p", Verbose=Verbose) 
        raw['Gamma12']=SwT[:,0,1].copy()
        raw['Gamma21']=SwT[:,1,0].copy()
    return raw

# Procedure for rotation of cartesion coordinates by radian angle
def rotate(x,y,phi):
    xp = x*np.cos(-phi)-y*np.sin(-phi)
    yp = x*np.sin(-phi)+y*np.cos(-phi)
    return xp,yp  

# Create Procedures to Calculate Reflection Coeff of the Open and Short
# For the CMT S911 cal kit, 50 ohm load is assumed to be perfect
def calcshortRho(freq, L0=0.0, L1=0.0, L2=0.0, L3=0.0, Offset=0.0):
    #L0 = 2.0765e-12  # H
    #L1 = -108.54e-24  # H/Hz
    #L2 = 2.1705e-33  # H/Hz^2
    #L3 = -0.01e-42  # H/Hz^3
   # Offs = 31.785e-12	# pS
    ind = L0+L1*freq+L2*freq**2+L3*freq**3
    XL = 2*np.pi*freq*ind*1j
    Rho = (XL-50)/(XL+50)
    phi = 2*np.pi*Offset*freq
    x = Rho.real
    y = Rho.imag
    xp,yp = rotate(x,y,phi)
    Rho = np.array([complex(xv,yv) for xv,yv in zip(xp,yp)])
    return Rho
  
def calcopenRho(freq, C0=0.0, C1=0.0, C2=0.0, C3=0.0, Offset=0.0):
    #C0 = 49.433e-15  # F
    #C1 = -310.13e-27  # F/Hz
    #C2 = 23.168e-36  # F/Hz^2
    #C3 = -0.15966e-45  # F/Hz^3
    #Offs = 29.243e-12	# pS
    cap = C0+C1*freq+C2*freq**2+C3*freq**3
    XC = 1/(2*np.pi*freq*cap*1j)
    Rho = (XC-50)/(XC+50)
    phi = 2*np.pi*Offset*freq
    x = Rho.real
    y = Rho.imag
    xp,yp = rotate(x,y,phi)
    Rho = np.array([complex(xv,yv) for xv,yv in zip(xp,yp)])
    return Rho

def calcloadRho(freq, Z0=50.0, L0=0.0, L1=0.0, L2=0.0, L3=0.0, Offset=0.0):
    ind = L0+L1*freq+L2*freq**2+L3*freq**3
    XL = 2*np.pi*freq*ind*1j    
    Z=Z0+XL
    Rho = (Z-50)/(Z+50)
    phi = 2*np.pi*Offset*freq
    x = Rho.real
    y = Rho.imag
    xp,yp = rotate(x,y,phi)
    Rho = np.array([complex(xv,yv) for xv,yv in zip(xp,yp)])
    return Rho



def OnePortCorrection(freq, st1, st2, st3, m1, m2, m3): # Create the matrices
    """
    v0.1 from Bryan Walker tuto
    ex.:
       sED,  sES,  sERT  = OnePortCorrection(freq, openRho, shortRho, loadRho, iss14['OpenP1'], iss14['ShortP1'], iss14['LoadP1'])
    """
    C  = np.zeros((3,3),dtype=complex)
    CH = np.zeros((3,3),dtype=complex)
    V  = np.zeros((3,1),dtype=complex)
    
    numPoints=len(freq)
    # Form the result matrices Directivity, Reflection tracking and Source Match
    E = np.zeros((3,1),dtype=complex)

    D = np.zeros((numPoints),dtype=complex)
    R = np.zeros((numPoints),dtype=complex)
    S = np.zeros((numPoints),dtype=complex)

    openRho=st1
    mopenRho=m1
    shortRho=st2
    mshortRho=m2
    loadRho=st3
    mloadRho=m3
    
    # Fill them with actual and measured values and do the calculation

    for i in range(numPoints):
        C[0,0] = openRho[i]
        C[0,1] = 1
        C[0,2] = openRho[i]*mopenRho[i]
        C[1,0] = shortRho[i]
        C[1,1] = 1
        C[1,2] = shortRho[i]*mshortRho[i]
        C[2,0] = loadRho[i]
        C[2,1] = 1
        C[2,2] = loadRho[i]*mloadRho[i]
        V[0] = mopenRho[i]
        V[1] = mshortRho[i]
        V[2] = mloadRho[i]
    
        CT = np.transpose(C)
        CH = np.conjugate(CT) #CH is the Hermetian of C
        E1 = CH.dot(C)
        E2 = np.linalg.inv(E1)
        E3 = CH.dot(V)
        E = E2.dot(E3)
    
        D[i] = E[1]
        S[i] = E[2]
        R[i] = E[0]+E[1]*E[2]
    return D, S, R


def TwoPortCorrection(freq, mrefl, mtrans, D, S, R ):
        
    numPoints=len(freq)
    # Form the result matrices Directivity, Reflection tracking and Source Match

    L = np.zeros((numPoints),dtype=complex) # load match
    TT = np.zeros((numPoints),dtype=complex) # transmission tracking
    
    # Fill them with actual and measured values and do the calculation

    for i in range(numPoints):
        De    = D[i]*S[i]-R[i]
        isol   = 0.0 # we neglect isolation, would be here
        L[i]  = (mrefl[i]-D[i])/(mrefl[i]*S[i]-De)
        TT[i] = (mtrans[i]-isol)*(1-S[i]*L[i])
    return L, TT


def SOLT_infinity(iss, Verbose=False, onSub=False):
    err=dict()
    freq=iss['F']
    err['F']=freq.copy()
    shortRho = calcshortRho(freq,L0=3.3e-12)
    if onSub:
        if Verbose: print("Using GSG 100um ON SUBSTRATE OPEN (3.6fF)")
        openRho  = calcopenRho(freq,C0=3.6e-15, C1=0.0, C2=0 , C3=0)
    else:
        if Verbose: print("Using GSG 100um ON SEPARATION OPEN (-6.5fF)")
        openRho  = calcopenRho(freq,C0=-6.5e-15, C1=0.0, C2=0 , C3=0)
    loadRho  = calcloadRho(freq,L0=-0.4e-12)
    err['ED'],   err['ES'],  err['ERT']  = OnePortCorrection(freq, openRho, shortRho, loadRho, iss['OpenP1'], iss['ShortP1'], iss['LoadP1'])
    err['EDr'],  err['ESr'], err['ERTr'] = OnePortCorrection(freq, openRho, shortRho, loadRho, iss['OpenP2'], iss['ShortP2'], iss['LoadP2'])
    err['EL'], err['ETT']   = TwoPortCorrection(freq, iss['Thru'][:,0,0], iss['Thru'][:,1,0],  err['ED'],   err['ES'],  err['ERT'] )
    err['ELr'], err['ETTr'] = TwoPortCorrection(freq, iss['Thru'][:,1,1], iss['Thru'][:,0,1],  err['EDr'],  err['ESr'], err['ERTr']) 
    err['EX']=err['EXr']=np.zeros(np.shape(err['ED']), dtype=np.cdouble)
    return err

def Apply12TermEC(sp, ET):
    """ 
    wrapper around the function Apply12TermCorrection.
    """
    s11=sp[:,0,0]
    s12=sp[:,0,1]
    s21=sp[:,1,0]
    s22=sp[:,1,1]    
    s11a, s12a, s21a, s22a=Apply12TermCorrection(s11, s12, s21, s22, ET['ED'], ET['EL'], ET['ES'], ET['EX'], ET['ETT'], ET['ERT'], ET['EDr'], ET['ELr'], ET['ESr'], ET['EXr'], ET['ETTr'], ET['ERTr'])
    sp2=np.empty(np.shape(sp), dtype=np.cdouble)
    sp2[:,0,0]=s11a
    sp2[:,0,1]=s12a
    sp2[:,1,0]=s21a
    sp2[:,1,1]=s22a    
    return sp2

def Apply12TermCorrection(s11, s12, s21, s22, ED, EL, ES, EX, ETT, ERT, EDr, ELr, ESr, EXr, ETTr, ERTr): 
    num=(s11-ED)/ERT*(1.0+(s22-EDr)*ESr/ERTr)
    num=num-EL*(s21-EX)*(s12-EXr)/ETT/ETTr
    den=(1.0+(s11-ED)/ERT*ES)*(1.0+(s22-EDr)/ERTr*ESr)
    den=den-ELr*EL*(s21-EX)*(s12-EXr)/ETT/ETTr

    s11a=num/den

    num=(s22-EDr)/ERTr*(1.0+(s11-ED)*ES/ERT)
    num=num-ELr*(s21-EX)*(s12-EXr)/ETT/ETTr

    s22a=num/den

    num=(s21-EX)/ETT*(1.0+(s22-EDr)*(ESr-EL)/ERTr)

    s21a=num/den

    num=(s12-EXr)/ETTr*(1.0+(s11-ED)*(ES-ELr)/ERT)

    s12a=num/den
    return s11a, s12a, s21a, s22a

def compare12ETsets(set1, set2, skipgamma=True, skipiso=True, CleanAngle=True,figsize=(20, 4)):

    func1= lambda x: sdb(x)
    lab1='[dB]' 
    if CleanAngle:
        func2= lambda x: CleanAngleDeg(deg(x))
    else:
        func2= lambda x: deg(x)
    lab2='deg'

    fw1=1e-9*set1['F']
    fw2=1e-9*set1['F']

    if np.all(fw1==fw2):
        cando=True
    else:
        print(f"frequency sweeps are not identical, cannot compare")
    for what in whatlist:
        if what=='EX' and skipiso: continue
        if what=='EXr' and skipiso: continue
        if what=='Gamma21' and skipgamma: continue
        if what=='Gamma12' and skipgamma: continue
        #fig= plt.figure(figsize=figsize)
        fig, (ax,bx,cx,dx) = plt.subplots(1, 4, figsize=figsize) # grid_minor_enable=True, grid_major_enable=True, grid_major_fancy=False, grid_minor_fancy=False)
        fig.suptitle(f"{what} - {whatname[what]}",fontsize=12)
        fig.subplots_adjust(wspace=0.35)
    
        ax.plot(fw1, func1(set1[what]),  c='r', marker='', linestyle=':', label='set1')   
        ax.plot(fw2, func1(set2[what]),  c='g', marker='', linestyle=':', label='set2') 
        ax.legend(); ax.set_xlabel(' Frequency [GHz]'); ax.set_ylabel(lab1);
        if cando:
            bx.plot(fw1, func1(set1[what])-func1(set2[what]),  c='b', marker='', linestyle=':',) 
            bx.set_xlabel(' Frequency [GHz]'); bx.set_ylabel(f" diff {lab1}");
    
        cx.plot(fw1, func2(set1[what]),  c='r', marker='', linestyle=':', label='set1')   
        cx.plot(fw2, func2(set2[what]),  c='g', marker='', linestyle=':', label='set2') 
        cx.legend(); cx.set_xlabel(' Frequency [GHz]'); cx.set_ylabel(lab2);

        if cando:
            dx.plot(fw1, func2(set1[what])-func2(set2[what]),  c='b', marker='', linestyle=':',) 
            dx.set_xlabel(' Frequency [GHz]'); dx.set_ylabel(f" diff {lab2}");

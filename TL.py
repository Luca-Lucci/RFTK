#!/usr/bin/env python
# coding: utf-8

# # TL

# ** DO NOT EDIT THIS .py FILE. EDIT THE JUPYTER NOTEBOOK **

# library with routines to manipulate Sparams of transmission lines.
# 
# To actually create the module, export as an executable python script

# Changelog

# - v 0.2: fix ax.grid(b=True) issue for python > 3.9
# - v 0.1: porting of first macros, debugged in  "http://localhost:8888/lab/tree/Investigations/211119_OnWaferCal/20_debug_TL_on_ISS_and_ADS_TL_object.ipynb"

# In[ ]:


# this line is useless, will be added by export filter
# -*- coding: utf-8 -*-

# import re   as re 
import numpy as np  # so far is the only dependable module to be included
import uW as uW
from common import *
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


# In[ ]:


class TL(object):
    def __init__(self, Freq, Sm):
        self.S=Sm
        self.F=np.array(Freq) # Hz
        self.FHz=self.F
        self.FGHz=self.F/1e9 # GHz
        self.omega=2*np.pi*self.F
        


# ADS uses a unwrap internal function that is not well documented (it is not implemented in AEL).
# At the beginning I tested the np.unwrap function but it is working in radiants and the period (at least in 
# this version) is limited to 2pi.
# Dalla funzione atanh i dati sono riportati solo con parte immaginaria nell'internallo -pi/2, pi/2,
# quindi ho bisogno di aggiungere una periodicita di pi o 180.
# sembra che in una nouva versione di matplotlib  il period possa essere specificato. 
# per il momento mi accontento di usare questa versione fast implemented
def unwrapDeg(n,jump=90):
    for i in range(len(n)-1):
        diff=np.abs(n[i+1]-n[i])
        if diff<jump: continue
        for j in range(i+1,len(n)):n[j]+=180
    return n
# In[2]:


def unwrapDeg(n,jump=90,stop=10):
    for i in range(len(n)-1):
        diff=np.abs(n[i+1]-n[i])
        while diff>jump and stop>0: 
            for j in range(i+1,len(n)):n[j]+=180
            diff=np.abs(n[i+1]-n[i])
            stop=stop-1
    return n


# ripoff of a ADS DDs setup to extract basic TL properties

# In[ ]:


class ADS(object):
    # this class implements functions and plot stolen from a DDS
    # template from ADS: Em_Product -> Z0RLGC_FromSParams
    # to get the Zc of the line, you do not have to specify
    # the Z0 of the system, like in Eo paper, since 50 Ohms are
    # given for granted ...
    # length is the line length in meters
    # vel is the estimated velocity in um/ps (i.e. 130 um/ps per il ISS cascade)
    # che serve a convertire la length in delay
    def __init__(self, Freq, Sm, length=1e-3,vel=130,name='TL'):
        self.S=Sm
        self.velocity=vel*1e-6/1e-12 # um/ps to m/s
        self.target_delay=length/self.velocity
        self.F=np.array(Freq)
        self.FGHz=self.F/1e9
        self.Y=uW.S2Y(self.S)
        self.Z=uW.S2Z(self.S)
        self.S11=self.S[:,0,0]
        self.S21=self.S[:,1,0]
        self.Z11=self.Z[:,0,0]
        self.Y11=self.Y[:,0,0]
        self.omega=2*np.pi*self.F
        self.l=length
        self.name=f"{name}"   
    @property
    def Zc(self):
        return np.sqrt(self.Z11/self.Y11)
    @property
    def GammaWrapped(self):
        # arctanh is a multivalued function: for each x there are infinitely many numbers z such that tanh(z) = x. 
        # The convention is to return the z whose imaginary part lies in [-pi/2, pi/2].
        # this causes a lot of problems because in extraction of Beta we have discontinuities
        # here beta is inside -pi/2, pi/2,
        # from Sp  compute  Y,Z matrix (in 50 ohm system then)
        # from Y,Z compute Zc, characteristic impedance.
        # from Y,Zc and line length l compute the Gamma in wrapped form (using atanh)
        # real part is ok, call it alpha
        # take imaginary part of GammaWrapped call it betaWrapped, all values are inside -pi/2,pi/2
        # unwrap so *it is monotonic ?* test with ADS if needed, and call it beta
        # call Gamma=alpha+j beta
        # compute R,L,C,G ....
        return 1/self.l*np.arctanh(1/(self.Zc*self.Y11))
    @property
    def alpha(self):
        # alphais no issue
        return np.real(self.GammaWrapped)
    @property
    def betaWrapped(self):
        # as said, here beta is inside -pi/2, pi/2,
        return np.imag(self.GammaWrapped)
    @property
    def beta(self):
        # ADS moltiplica beta per l e converte in degrees
        #return np.unwrap(self.betaWrapped*self.l*180/np.pi)*np.pi/180/self.l
        #return np.unwrap(self.betaWrapped*self.l)/self.l
        return unwrapDeg(self.betaWrapped*self.l*180/np.pi)*np.pi/180/self.l
    @property
    def Gamma(self):
        return self.alpha+1j*self.beta
    @property
    def R(self):
        # Ohm/m
        return np.real(self.Gamma*self.Zc)
    @property
    def L(self):
        # H/m
        return np.imag(self.Gamma*self.Zc)/self.omega
    @property
    def G(self):
        # S/m
        return np.real(self.Gamma/self.Zc)
    @property
    def C(self):
        # F/m
        return np.imag(self.Gamma/self.Zc)/self.omega
    @property
    def att_dB_per_m(self):
        return 20.0*self.alpha*np.log10(np.exp(1))
    @property
    def att_dB(self):
        return self.att_dB_per_m*self.l
    @property
    def phaseShift(self):
        return self.beta*self.l*180/np.pi
    @property
    def velocity_m_per_s(self):
        bm=self.omega>1e-5
        t1=self.omega/self.beta
        t2=1/self.beta
        return np.where(bm, t1, t2)
    @property
    def delay_ps_per_m(self):
        return 1e12/self.velocity_m_per_s
    @property
    def delay_ps(self):
        return self.l*self.delay_ps_per_m
    @property
    def eff_dielectric_const(self):
        c0=299_792_458 # m / s
        return (c0/self.velocity_m_per_s)**2
    @property
    def delay(self):
        return self.delay_ps*1e-12
    @property
    def v(self):
        return self.velocity_m_per_s
    @property
    def eps(self):
        return self.eff_dielectric_const
    def plotR(self,cmp=[],ylim="auto"):
        fig, ax =plt.subplots(ncols=1,nrows=1,figsize=(5,4))
        ax.plot(self.FGHz, self.R, label=self.name)
        for TL in cmp:
            ax.plot(TL.FGHz, TL.R, label=TL.name)
        ax.minorticks_on()
        ax.grid(which='major',ls='-')
        ax.grid(which='minor',ls=':')
        ax.set_xlabel(f"Freq [GHz]");
        ax.set_ylabel(f"R [Ohm/m]"); 
        if ylim=="auto": pass
        else: ax.set_ylim(ylim)
        ax.legend()
        return fig, ax        
    def plotC(self,cmp=[],ylim="auto"):
        fig, ax =plt.subplots(ncols=1,nrows=1,figsize=(5,4))
        ax.plot(self.FGHz, self.C, label=self.name)
        for TL in cmp:
            ax.plot(TL.FGHz, TL.C, label=TL.name)
        ax.minorticks_on()
        ax.grid(which='major',ls='-')
        ax.grid(which='minor',ls=':')
        ax.set_xlabel(f"Freq [GHz]");
        ax.set_ylabel(f"C [F/m]");        
        if ylim=="auto": pass
        else: ax.set_ylim(ylim)       
        ax.legend()
        return fig, ax
    def plotL(self,cmp=[],ylim="auto"):
        fig, ax =plt.subplots(ncols=1,nrows=1,figsize=(5,4))
        ax.plot(self.FGHz, self.L, label=self.name)
        for TL in cmp:
            ax.plot(TL.FGHz, TL.L, label=TL.name)
        ax.minorticks_on()
        ax.grid(which='major',ls='-')
        ax.grid(which='minor',ls=':')
        ax.set_xlabel(f"Freq [GHz]");
        ax.set_ylabel(f"L [H/m]");        
        if ylim=="auto": pass
        else: ax.set_ylim(ylim)        
        ax.legend()
        return fig, ax        
    def plotG(self,cmp=[],ylim="auto"):
        fig, ax =plt.subplots(ncols=1,nrows=1,figsize=(5,4))
        ax.plot(self.FGHz, self.G, label=self.name)
        for TL in cmp:
            ax.plot(TL.FGHz, TL.G, label=TL.name)
        ax.minorticks_on()
        ax.grid(which='major',ls='-')
        ax.grid(which='minor',ls=':')
        ax.set_xlabel(f"Freq [GHz]");
        ax.set_ylabel(f"G [S/m]");        
        if ylim=="auto": pass
        else: ax.set_ylim(ylim)
        ax.legend()
        return fig, ax        
    def plotRLGC(self,cmp=[]):
        self.plotR(cmp=cmp)
        self.plotL(cmp=cmp)
        self.plotC(cmp=cmp)
        self.plotG(cmp=cmp)
    def plotZc(self,cmp=[], line50=True,ylim=(45,55)):
        fig, ax =plt.subplots(ncols=1,nrows=1,figsize=(5,4))
        ax.plot(self.FGHz, np.abs(self.Zc), label=self.name)
        for TL in cmp:
            ax.plot(TL.FGHz, np.abs(TL.Zc), label=TL.name)
        ax.minorticks_on()
        ax.grid(which='major',ls='-')
        ax.grid(which='minor',ls=':')
        if line50: ax.axhline(y=50, color='r', linestyle='-')
        ax.set_xlabel(f"Freq [GHz]");
        ax.set_ylabel(f"abs(Zc) [Ohm]");
        if ylim=="auto": pass
        else: ax.set_ylim(ylim)        
        ax.legend()
        return fig, ax
    def plotAtt(self,cmp=[],norm=True):
        fig, ax =plt.subplots(ncols=1,nrows=1,figsize=(5,4))
        if norm: ax.plot(self.FGHz, self.att_dB_per_m, label=self.name)
        else: ax.plot(self.FGHz, self.att_dB, label=self.name)
        for TL in cmp:
            if norm: ax.plot(TL.FGHz, TL.att_dB_per_m, label=TL.name)
            else: ax.plot(TL.FGHz, TL.att_dB, label=TL.name)
        ax.minorticks_on()
        ax.grid(which='major',ls='-')
        ax.grid(which='minor',ls=':')
        ax.set_xlabel(f"Freq [GHz]");
        if norm: ax.set_ylabel(f"att [dB/m]")
        else:    ax.set_ylabel(f"att [dB]")
        if cmp==[]:
            ax.set_title(f"for line length {self.l*1e3} mm");
        ax.legend()
        return fig, ax
    def plotDelay(self,cmp=[],norm=False,line=True):
        fig, ax =plt.subplots(ncols=1,nrows=1,figsize=(5,4))
        if norm:
            ax.plot(self.FGHz, self.delay/self.target_delay, 
                        label= self.name if cmp==[] else f"{self.name} l={self.l*1e3} mm")
        else: 
            ax.plot(self.FGHz, self.delay_ps, label=self.name)
            if line: ax.axhline(y=self.target_delay*1e12, color='r', linestyle=':') # ps
                
        for TL in cmp:
            if norm:
                ax.plot(TL.FGHz, TL.delay/TL.target_delay, label=f"{TL.name} l={TL.l*1e3} mm")
            else:
                ax.plot(TL.FGHz, TL.delay_ps, label=TL.name)
                if line: ax.axhline(y=TL.target_delay*1e12, color='r', linestyle=':')
        ax.minorticks_on()
        ax.grid(which='major',ls='-')
        ax.grid(which='minor',ls=':')
        ax.set_xlabel(f"Freq [GHz]");
        if norm:
            ax.set_ylabel(f"delay/target delay []")
        else:
            ax.set_ylabel(f"delay [ps]")   
        if cmp==[]:
            ax.set_title(f"for line length {self.l*1e3} mm");
        ax.legend()
        return fig, ax
    
    def plotEps(self,cmp=[]):
        fig, ax =plt.subplots(ncols=1,nrows=1,figsize=(5,4))
        ax.plot(self.FGHz, self.eps, label=self.name)
        for TL in cmp:
            ax.plot(TL.FGHz, TL.eps, label=TL.name)
        ax.minorticks_on()
        ax.grid(which='major',ls='-')
        ax.grid(which='minor',ls=':')
        ax.set_xlabel(f"Freq [GHz]");
        ax.set_ylabel(f"effective dielectric constant []")
        ax.legend()
        return fig, ax
    def plotPhaseShift(self,cmp=[]):
        fig, ax =plt.subplots(ncols=1,nrows=1,figsize=(5,4))
        ax.plot(self.FGHz, self.phaseShift, label=self.name)
        for TL in cmp:
            ax.plot(TL.FGHz, TL.phaseShift, label=TL.name)
        ax.minorticks_on()
        ax.grid(which='major',ls='-')
        ax.grid(which='minor',ls=':')
        ax.set_xlabel(f"Freq [GHz]");
        ax.set_ylabel(f" phase shift [deg]")
        ax.legend()
        return fig, ax
    def plots(self,cmp=[]):
        self.plotZc(cmp=cmp)
        self.plotRLGC(cmp=cmp)
        self.plotAtt(cmp=cmp)
        self.plotDelay(cmp=cmp)
        self.plotEps(cmp=cmp)
        self.plotPhaseShift(cmp=cmp)
    def plotSxy(self, title="self.name", CleanAngle=False, logx=False, GHz=False):
        if title=="self.name":
            title=self.name
        plotSxy(self.F,self.S, title=title, CleanAngle=CleanAngle, logx=logx, GHz=GHz)


# this implementation of the Eo extraction (the original paper wher ethe extraction were proposed) is more sperimental.
# I kinds of matches till a certain frequency with the ADS extraction. Then problem in phase or value of Gamma emerges and some aspects of the lines are competely crazy.

# ![image.png](attachment:2c710608-86c2-4327-ae02-308597a3c20e.png)

# In[ ]:


class Eo(ADS):
    # this class implements functions and plot stolen from a DDS
    # template from ADS: Em_Product -> Z0RLGC_FromSParams
    # but implements Zc and Gamma from original Eo paper
  
    def __init__(self, Freq, Sm, Zref=50, length=1e-3,name='TL', forceGammaPositive=True):
        self.S=Sm
        self.Z0=Zref
        self.F=np.array(Freq)
        self.FGHz=self.F/1e9
        self.Y=uW.S2Y(self.S)
        self.Z=uW.S2Z(self.S)
        self.S11=self.S[:,0,0]
        self.S21=self.S[:,1,0]
        self.S22=self.S[:,1,1]
        self.Z11=self.Z[:,0,0]
        self.Y11=self.Y[:,0,0]
        self.omega=2*np.pi*self.F
        self.l=length
        self.name=f"{name}"   
        self.PositiveGamma=forceGammaPositive
    @property
    def Zc(self):
        # return np.sqrt(self.Z11/self.Y11)
        #     Z=Z0*np.sqrt( ( (1+s11)**2-(s21)**2 ) / ( (1-s11)**2-(s21)**2 ) )
        return self.Z0*np.sqrt( ( (1+self.S11)**2-(self.S21)**2 ) / ( (1-self.S11)**2-(self.S21)**2 ) )
    @property
    def GammaWrapped(self):
    #def GammaWrapped(S, l=100e-6, forcePositive=True):
        gamma=np.zeros(len(self.F),dtype=np.complex_)   
        for n in range(len(self.F)):
            s11=self.S11[n]
            s22=self.S22[n]
            s21=self.S21[n]
            k=np.sqrt( ( (s22**2-s21**2+1)**2-(2*s11)**2 )/( 2*s21)**2 )
            tmp=((1-s11**2+s21**2)/(2*s21))+k
            tmp=-np.log(1/tmp)
            gamma[n]=tmp/self.l
            if np.real(gamma[n]<0) and self.PositiveGamma:
                tmp=((1-s11**2+s21**2)/(2*s21))-k
                tmp=-np.log(1/tmp)
                gamma[n]=tmp/self.l
            if np.imag(gamma[n])<0 and self.PositiveGamma:
                gamma[n]=np.real(gamma[n])-np.imag(gamma[n])*1j                
        return gamma


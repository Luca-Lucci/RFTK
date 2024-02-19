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
import uW
import ICCAP
import calion as cal
import rftools_noise as rftn 
import rfdataanalysis as rfda
import rfplots as rfp
from smithplot import SmithAxes

#  from psfascii import psfascii ### linux only !?!?
#  import skrf as rf
#%matplotlib qt
import matplotlib as mpl
mpl.rcParams['savefig.dpi']=200  # stp 75


ceacolors=[(184/256, 20/256, 32/256), 
           (115/256, 190/256, 75/256),
           (72/256,172/256,240/256), 
           (160/256, 26/256, 125/256), 
           (229/256, 149/256, 0), 
           (23/256, 26/256, 33/256)]

marker=list('o+x^v<>sdD*phD12348')
linestyle=ls=['-','--','-.',':','']

figsz=(6,6)

av1=dict(color='red',
         #alpha=0.9,
         #marker=marker[1],
         linestyle=linestyle[0],
         linewidth=2);
av2=dict(color='blue',
         #alpha=0.9,
         #marker=marker[1],
         linestyle=linestyle[0],
         linewidth=2);

av3=dict(color='green',
         #alpha=0.9,
         #marker=marker[1],
         linestyle=linestyle[0],
         linewidth=2);


av4=dict(color='purple',
         #alpha=0.9,
         #marker=marker[1],
         linestyle=linestyle[0],
         linewidth=2);

av1b=dict(color='red',
         #alpha=0.9,
         #marker=marker[1],
         linestyle=linestyle[1],
         linewidth=1);
av2b=dict(color='blue', 
         #alpha=0.9,
         #marker=marker[1],
         linestyle=linestyle[1],
         linewidth=1);
av3b=dict(color='green', 
         #alpha=0.9,
         #marker=marker[1],
         linestyle=linestyle[1],
         linewidth=1);

av4b=dict(color='purple', 
         #alpha=0.9,
         #marker=marker[1],
         linestyle=linestyle[1],
         linewidth=1);

av1c=dict(color='red',
         #alpha=0.6,
         #marker=marker[1],
         linestyle=linestyle[3],
         linewidth=1);
av2c=dict(color='blue', 
         #alpha=0.6,
         #marker=marker[1],
         linestyle=linestyle[3],
         linewidth=1);

av3c=dict(color='green', 
         #alpha=0.6,
         #marker=marker[1],
         linestyle=linestyle[3],
         linewidth=1);

av4c=dict(color='purple', 
         #alpha=0.6,
         #marker=marker[1],
         linestyle=linestyle[3],
         linewidth=1);

av1d=dict(color='green',
         #alpha=0.9,
         #marker=marker[1],
         linestyle=linestyle[3],
         linewidth=2);
av2d=dict(color='green', 
         #alpha=0.9,
         #marker=marker[1],
         linestyle=linestyle[3],
         linewidth=2);

refline=dict(linewidth=2, color='g',ls='--',alpha=0.9)
simline=dict(color='red', linestyle=linestyle[1], alpha=1,   linewidth=1);
simlineg=dict(color='green', linestyle=linestyle[1], alpha=1,   linewidth=1);
simlineb=dict(color='blue', linestyle=linestyle[1], alpha=1,   linewidth=1);




# nifty utilities
d2r   = lambda x: np.deg2rad(x)
r2d   = lambda x: np.rad2deg(x)
polar = lambda m,r: m*np.exp(1j*d2r(r))
ralop = lambda z: (np.abs(z),r2d(np.angle(z)))
cart  = lambda z: (np.abs(z),np.angle(z))
conj  = lambda z: (np.conjugate(z))
sdb   = lambda z: (20.0*np.log10(np.abs(z)))
deg   = lambda z: (np.angle(z, deg=True))

def CleanAngleDeg(n):
    for i in range(len(n)-1):
        diff=np.abs(n[i+1]-n[i])
        if diff<50: continue
        if np.abs(n[i+1]-n[i]+360)<diff:
            for j in range(i+1,len(n)):n[j]+=360
        if np.abs(n[i+1]-n[i]-360)<diff:
            for j in range(i+1,len(n)):n[j]-=360
    return n
   


## common simple 2port definitions
# from Franz Sischa pag 10
def C10(Y,f):
    return np.imag(Y[:,0,0]+Y[:,0,1])/(2*np.pi*f)

def C20(Y,f):
    return np.imag(Y[:,1,1]+Y[:,0,1])/(2*np.pi*f)

def L10(Y,f):
    return np.imag((Y[:,0,0]+Y[:,0,1])**-1)/(2*np.pi*f)

def L20(Y,f):
    return np.imag((Y[:,1,1]+Y[:,0,1])**-1)/(2*np.pi*f)

def R10(Y,f):
    return np.real((Y[:,0,0]+Y[:,0,1])**-1)

def R20(Y,f):
    return np.real((Y[:,1,1]+Y[:,0,1])**-1)

def Thru_Delay(S,f):
    return -np.angle(S[:,1,0])/(2*np.pi*f)

def Thru_Delay_rev(S,f):
    return -np.angle(S[:,0,1])/(2*np.pi*f)


def Thru_Z0(S):
    A=rft.StoABCD(S)
    return np.sqrt(A[:,0,1]/A[:,1,0]) 

###### from ideal circuits
def C11(Y,f):
    return -np.imag((Y[:,0,0]+Y[:,0,1])**-1)**-1/(2*np.pi*f)


def C22(Y,f):
    return -np.imag((Y[:,1,1]+Y[:,0,1])**-1)**-1/(2*np.pi*f)

def C12(Y,f):
    return np.imag((Y[:,0,1])**-1)**-1/(2*np.pi*f)

def C21(Y,f):
    return np.imag((Y[:,1,0])**-1)**-1/(2*np.pi*f)

def R1(Y):
    return np.real(rftools.YtoZ_2P(Y)[:,0,0]-rftools.YtoZ_2P(Y)[:,0,1])

def R2(Y):
    return np.real(rftools.YtoZ_2P(Y)[:,1,1]-rftools.YtoZ_2P(Y)[:,1,0])

def L1(Y,f):
    return np.imag(rftools.YtoZ_2P(Y)[:,0,0]-rftools.YtoZ_2P(Y)[:,0,1])/(2*np.pi*f)

def L2(Y,f):
    return np.imag(rftools.YtoZ_2P(Y)[:,1,1]-rftools.YtoZ_2P(Y)[:,1,0])/(2*np.pi*f)



def plotSxy(FsweepHz, Sp, title="", CleanAngle=False, logx=False, GHz=False):
    fig = plt.figure(figsize=(16, 12))
    fig.suptitle(title,fontsize=24)
    ax  = plt.subplot(2, 3, 1, projection='smith', grid_minor_enable=True, grid_major_enable=True, grid_major_fancy=False,
                     grid_minor_fancy=False,grid_major_fancy_threshold=(50, 50))
    bx  = plt.subplot(2, 3, 2)
    cx  = plt.subplot(2, 3, 3)
    dx  = plt.subplot(2, 3, 5)
    ex  = plt.subplot(2, 3, 6)
    
    if GHz:
        F=FsweepHz*1e-9
    else:
        F=FsweepHz
    
    if logx:
        bxplot=bx.semilogx
        cxplot=cx.semilogx
        dxplot=dx.semilogx
        explot=ex.semilogx
    else:
        bxplot=bx.plot
        cxplot=cx.plot
        dxplot=dx.plot    
        explot=ex.plot            

    ax.plot(Sp[:,0,0], linestyle='-',marker="", label="s11", datatype=SmithAxes.S_PARAMETER)
    ax.plot(Sp[:,1,1], linestyle='-',marker="", label="s22", datatype=SmithAxes.S_PARAMETER)

    bxplot(F,sdb(Sp[:,0,0]), linestyle='-',marker="", markevery=10, label='s11')
    bxplot(F,sdb(Sp[:,1,1]), linestyle='-',marker="",  label="s22")

    cxplot(F,sdb(Sp[:,1,0]), linestyle='-',marker="", markevery=10, label='s21')
    cxplot(F,sdb(Sp[:,0,1]), linestyle='-',marker="",  label="s12")

    phase11=np.angle(Sp[:,0,0],deg=True)
    phase22=np.angle(Sp[:,1,1],deg=True)    
    if CleanAngle:
        phase11=CleanAngleDeg(phase11)
        phase22=CleanAngleDeg(phase22)
    dxplot(F, phase11,  linestyle='-',marker="", label='s11')
    dxplot(F,phase22, linestyle='-',marker="",  label="s22")

    phase21=np.angle(Sp[:,1,0],deg=True)
    phase12=np.angle(Sp[:,0,1],deg=True)
    if CleanAngle:
        phase21=CleanAngleDeg(phase21)
        phase12=CleanAngleDeg(phase12)
    explot(F,phase21, linestyle='-',marker="",  label='s21')
    explot(F,phase12, linestyle='-',marker="",  label="s12")
    #ax.legend(loc="best", fontsize=12); 
    bx.set_title("logMag S11/S22")
    bx.set_ylabel("dB")
    dx.set_title("phase(deg) S11/S22")
    dx.set_ylabel("deg")
    cx.set_title("logMag S21/S12")
    cx.set_ylabel("dB")
    ex.set_title("phase(deg) S21/S12")
    ex.set_ylabel("deg")
    
    for zx in [bx, cx, dx, ex]:
        zx.legend(loc="best", fontsize=12)
        zx.minorticks_on()
        zx.grid(which='major',ls='-')
        zx.grid(which='minor',ls=':')
        if GHz:
            zx.set_xlabel("F [GHz]");
        else:
            zx.set_xlabel("F [Hz]");
    fig.subplots_adjust(wspace=0.25) #cx.legend(loc='best')
    return fig

def plotSxySxy(FsweepHz, Sp, Spp, title="", CleanAngle=False, logx=False, GHz=False):
    c1='tab:blue'
    c2='tab:orange'
    fig = plt.figure(figsize=(16, 12))
    fig.suptitle(title,fontsize=24)
    ax  = plt.subplot(2, 3, 1, projection='smith', grid_minor_enable=True, grid_major_enable=True, grid_major_fancy=False,
                     grid_minor_fancy=False,grid_major_fancy_threshold=(50, 50))
    bx  = plt.subplot(2, 3, 2)
    cx  = plt.subplot(2, 3, 3)
    dx  = plt.subplot(2, 3, 5)
    ex  = plt.subplot(2, 3, 6)
    
    if GHz:
        F=FsweepHz*1e-9
    else:
        F=FsweepHz
    
    if logx:
        bxplot=bx.semilogx
        cxplot=cx.semilogx
        dxplot=dx.semilogx
        explot=ex.semilogx
    else:
        bxplot=bx.plot
        cxplot=cx.plot
        dxplot=dx.plot    
        explot=ex.plot            

    ax.plot(Sp[:,0,0], linestyle='-',marker="", color=c1, label="s11", datatype=SmithAxes.S_PARAMETER)
    ax.plot(Sp[:,1,1], linestyle='-',marker="", color=c2, label="s22", datatype=SmithAxes.S_PARAMETER)
    ax.plot(Spp[:,0,0], linestyle=':',marker="", color=c1, label="", datatype=SmithAxes.S_PARAMETER)
    ax.plot(Spp[:,1,1], linestyle=':',marker="", color=c2, label="", datatype=SmithAxes.S_PARAMETER)

    bxplot(F,sdb(Sp[:,0,0]), linestyle='-',marker="", color=c1,  label='s11')
    bxplot(F,sdb(Sp[:,1,1]), linestyle='-',marker="", color=c2,  label="s22")
    bxplot(F,sdb(Spp[:,0,0]), linestyle=':',marker="", color=c1, label='')
    bxplot(F,sdb(Spp[:,1,1]), linestyle=':',marker="", color=c2,  label="")

    cxplot(F,sdb(Sp[:,1,0]), linestyle='-',marker="", color=c1,  label='s21')
    cxplot(F,sdb(Sp[:,0,1]), linestyle='-',marker="", color=c2,  label="s12")
    cxplot(F,sdb(Spp[:,1,0]), linestyle=':',marker="", color=c1,  label='')
    cxplot(F,sdb(Spp[:,0,1]), linestyle=':',marker="", color=c2,  label="")

    phase11=np.angle(Sp[:,0,0],deg=True)
    phase22=np.angle(Sp[:,1,1],deg=True)    
    phase11p=np.angle(Spp[:,0,0],deg=True)
    phase22p=np.angle(Spp[:,1,1],deg=True)    
    if CleanAngle:
        phase11=CleanAngleDeg(phase11)
        phase22=CleanAngleDeg(phase22)
        phase11p=CleanAngleDeg(phase11p)
        phase22p=CleanAngleDeg(phase22p)
    dxplot(F,phase11,  linestyle='-',marker="", color=c1, label='s11')
    dxplot(F,phase22, linestyle='-',marker="", color=c2,  label="s22")
    dxplot(F,phase11p,  linestyle=':',marker="", color=c1, label='')
    dxplot(F,phase22p, linestyle=':',marker="", color=c2,  label="")

    phase21=np.angle(Sp[:,1,0],deg=True)
    phase12=np.angle(Sp[:,0,1],deg=True)
    phase21p=np.angle(Spp[:,1,0],deg=True)
    phase12p=np.angle(Spp[:,0,1],deg=True)
    if CleanAngle:
        phase21=CleanAngleDeg(phase21)
        phase12=CleanAngleDeg(phase12)
        phase21p=CleanAngleDeg(phase21p)
        phase12p=CleanAngleDeg(phase12p)
    explot(F,phase21, linestyle='-',marker="", color=c1,  label='s21')
    explot(F,phase12, linestyle='-',marker="", color=c2,  label="s12")
    explot(F,phase21p, linestyle=':',marker="", color=c1,  label='')
    explot(F,phase12p, linestyle=':',marker="", color=c2,  label="")
    #ax.legend(loc="best", fontsize=12); 
    bx.set_title("logMag S11/S22")
    bx.set_ylabel("dB")
    dx.set_title("phase(deg) S11/S22")
    dx.set_ylabel("deg")
    cx.set_title("logMag S21/S12")
    cx.set_ylabel("dB")
    ex.set_title("phase(deg) S21/S12")
    ex.set_ylabel("deg")
    
    for zx in [bx, cx, dx, ex]:
        zx.legend(loc="best", fontsize=12)
        zx.minorticks_on()
        zx.grid(which='major',ls='-')
        zx.grid(which='minor',ls=':')
        if GHz:
            zx.set_xlabel("F [GHz]");
        else:
            zx.set_xlabel("F [Hz]");
    fig.subplots_adjust(wspace=0.25) #cx.legend(loc='best')
    return fig

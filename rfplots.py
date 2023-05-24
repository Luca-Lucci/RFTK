# -*- coding: utf-8 -*-
"""
Created on Sun Apr  8 22:54:48 2018
"""

from statistics import mean
# import skrf as rf
import re   as re
import numpy as np
import scipy as sp
from scipy import stats
import pickle 
import gzip
import subprocess
import copy
import matplotlib.pyplot as plt
import rftools
import smithplot
import matplotlib as mpl

siz=14
mpl.rc('xtick', labelsize=siz) 
mpl.rc('ytick', labelsize=siz)



def trace_linear(dev, x_axis, y_axis, curve_label, xlabel, ylabel):
    fig, ax = plt.subplots(figsize=(7,5));
    line2d=ax.plot(x_axis,y_axis,label=curve_label)
    ax.set_title(dev)
    ax.legend(loc='best',framealpha=0.7, prop={'size':'medium'},numpoints=1) #   title='Vg',    numpoints=1
    ax.grid(b=True, which='both', axis='both') #, **kwargs)
    ax.set_xlabel(xlabel, fontweight='bold', size=siz); ax.set_ylabel(ylabel, fontweight='bold', size=siz)
    #fig.savefig("101_FigsDambrine/01_UFM2_FG.png",dpi=175,bbox_inches='tight')
    return line2d

def trace_semilogy(dev, x_axis, y_axis, curve_label, xlabel, ylabel):
    fig, ax = plt.subplots(figsize=(7,5));
    line2d=ax.semilogy(x_axis,y_axis,label=curve_label)
    ax.set_title(dev)
    ax.legend(loc='best',framealpha=0.7, prop={'size':'medium'},numpoints=1) #   title='Vg',    numpoints=1
    ax.grid(b=True, which='both', axis='both') #, **kwargs)
    ax.set_xlabel(xlabel, fontweight='bold', size=siz); ax.set_ylabel(ylabel, fontweight='bold', size=siz)
    #fig.savefig("101_FigsDambrine/01_UFM2_FG.png",dpi=175,bbox_inches='tight')
    return line2d

def trace_Fts (dev, Vg, Ftgm_Id, Ft_currentgain, Ft_calc_RF, Ft_calc_DC, Ft_complex_RF):
    fig, cx = plt.subplots(figsize=(7,13));
    cx.plot(Vg,Ft_currentgain,label='Ft_current_gain')
    cx.plot(Vg,Ft_calc_RF,label='Ft_calc_RF')
    cx.plot(Vg,Ft_calc_DC,label='Ft_calc_DC')
    cx.plot(Vg,Ft_complex_RF,label='Ft_complex_RF')
    cx.set_title(dev)
    cx.legend(loc='best',framealpha=0.7, prop={'size':'medium'},numpoints=1) #   title='Vg',    numpoints=1
    cx.grid(b=True, which='both', axis='both') #, **kwargs)
    cx.set_xlabel('Vg[V]', fontweight='bold', size=siz); cx.set_ylabel('Ft [GHz]', fontweight='bold', size=siz)
    
def trace_Fmaxs (dev, Vg, Fmax_powergain, Fmax_calc_RF, Fmax_complex_RF):   
    fig, ex = plt.subplots(figsize=(7,13));
    ex.plot(Vg,Fmax_powergain,label='Fmax_power_gain')
    ex.plot(Vg,Fmax_calc_RF,label='Fmax_calc_RF')
    ex.plot(Vg,Fmax_complex_RF,label='Fmax_complex_RF')
    ex.set_title(dev)
    ex.legend(loc='best',framealpha=0.7, prop={'size':'medium'},numpoints=1) #   title='Vg',    numpoints=1
    ex.grid(b=True, which='both', axis='both') #, **kwargs)
    ex.set_xlabel('Vg [V]', fontweight='bold', size=siz); ex.set_ylabel('Fmax [GHz]', fontweight='bold', size=siz)
    
def trace_Cgd (dev, Vg, Cgd):    
    fig,  bx = plt.subplots(figsize=(7,13));
    bx.plot(Vg,Cgd,label='Cgd')
    bx.set_title(dev)
    bx.legend(loc='best',framealpha=0.7, prop={'size':'medium'},numpoints=1) #   title='Vg',    numpoints=1
    bx.grid(b=True, which='both', axis='both') #, **kwargs)
    bx.set_xlabel('Vg[V]', fontweight='bold', size=siz); bx.set_ylabel('Cgd [F]', fontweight='bold', size=siz)
    
def trace_Cgs (dev, Vg, Cgs):
    
    fig, ax = plt.subplots(figsize=(7,13));
    ax.plot (Vg,Cgs,label='Cgs')
    ax.set_title(dev)
    ax.legend(loc='best',framealpha=0.7, prop={'size':'medium'},numpoints=1) #   title='Vg',    numpoints=1
    ax.grid(b=True, which='both', axis='both') #, **kwargs)
    ax.set_xlabel('Vg[V]', fontweight='bold', size=siz); ax.set_ylabel('Cgs [F]', fontweight='bold', size=siz)
    
def trace_FtGm_Id (dev, Id, Ftgm_Td):
    
    fig, dx = plt.subplots(figsize=(7,13));
    dx.semilogx(Id,Ftgm_Id,label='Ftgm_Id')
    dx.set_title(dev)
    dx.legend(loc='best',framealpha=0.7, prop={'size':'medium'},numpoints=1) #   title='Vg',    numpoints=1
    dx.grid(b=True, which='both', axis='both') #, **kwargs)
    dx.set_xlabel('Id', fontweight='bold', size=siz); dx.set_ylabel('Ft [GHz]', fontweight='bold', size=siz)
   

def trace_resistances_bracale_method(dev, Denom, Rg, Rd, Rs, ShowZeroOnX=False):
    fig, ax, = plt.subplots(figsize=(3,4));
    ax.plot(Denom,Rg,"ro",linestyle='dashed', label='Re(Z11-Z12)')
    ax.plot(Denom,Rd,"bo",label='Re(Z22-Z12)')
    ax.plot(Denom,Rs,"go",label='Re(Z12)')
    #ax.set_title(dev)
    ax.minorticks_on()
    ax.grid(b=True, which='major',ls='-')
    ax.grid(b=True, which='minor',ls=':')
    if ShowZeroOnX: ax.set_xlim(left=-0.1)
    ax.legend(loc='best',framealpha=0.7, prop={'size':'small'},numpoints=1) #   title='Vg',    numpoints=1
    #ax.grid(b=True, which='both', axis='both') #, **kwargs)
    ax.set_xlabel('1/(Vgs-Vth) [1/V]') #, fontweight='bold', size=siz); 
    ax.set_ylabel('Resistances [ohms]') # , fontweight='bold', size=siz)
    #fig.savefig("figs/06_Res_Bracale.png",dpi=175,bbox_inches='tight')
    

def trace_RgVgs (dev, Vg, RgEnz, RgJen, RgDormieu):
    fig, bx = plt.subplots(figsize=(7,5));    
    bx.plot(Vg,RgEnz,label='Rg_Enz')
    bx.plot(Vg,RgJen,label='Rg_Jen')
    bx.plot(Vg,RgDormieu,label='Rg_Dormieu')
    bx.set_title(dev)
    bx.legend(loc='best',framealpha=0.7, prop={'size':'medium'},numpoints=1) #   title='Vg',    numpoints=1
    bx.grid(b=True, which='both', axis='both') #, **kwargs)
    bx.set_xlabel('Vg [V]', fontweight='bold', size=siz); bx.set_ylabel('Resistances [ohms]', fontweight='bold', size=siz)
    
def trace_inductances (dev, Denom, Lg, Ld, Ls,  ShowZeroOnX=False, pH=False):
    fig, ax = plt.subplots(figsize=(7,5));
    if pH: 
        fac=1e12
    else:
        fac=1
    ax.plot(Denom,fac*np.array(Lg),"ro",label='Imag(Z11-Z12)')
    ax.plot(Denom,fac*np.array(Ld),"bo",label='Imag(Z22-Z12) ')
    ax.plot(Denom,fac*np.array(Ls),"go",label='Imag(Z12) ')
    ax.set_title(dev)
    if ShowZeroOnX: ax.set_xlim(left=-0.1)
    ax.legend(loc='best',framealpha=0.7, prop={'size':'medium'},numpoints=1) #   title='Vg',    numpoints=1
    ax.grid(b=True, which='both', axis='both') #, **kwargs)
    ax.set_xlabel('1/(Vgs-Vth)² (1/V²)', fontweight='bold', size=siz); 
    if pH:
        label='Inductance (pH)'
    else:
        label='Inductance (H)'
    ax.set_ylabel(label, fontweight='bold', size=siz)
    #fig.savefig("101_FigsDambrine/01_UFM2_FG.png",dpi=175,bbox_inches='tight')


def trace_inductances_bracale_method(dev, Denom, Lg, Ld, Ls,  ShowZeroOnX=False, pH=False):
    fig, ax = plt.subplots(figsize=(3,4));
    if pH: 
        fac=1e12
    else:
        fac=1
    ax.plot(Denom,fac*np.array(Lg),"ro",label='Im(Z11-Z12)')
    ax.plot(Denom,fac*np.array(Ld),"bo",label='Im(Z22-Z12) ')
    ax.plot(Denom,fac*np.array(Ls),"go",label='Im(Z12) ')
    #ax.set_title(dev)
    if ShowZeroOnX: ax.set_xlim(left=-0.1)
    ax.legend(loc='best',framealpha=0.7, prop={'size':'small'},numpoints=1) #   title='Vg',    numpoints=1
    #ax.grid(b=True, which='both', axis='both') #, **kwargs)
    ax.minorticks_on()
    ax.grid(b=True, which='major',ls='-')
    ax.grid(b=True, which='minor',ls=':')
    ax.set_xlabel('1/(Vgs-Vth)$^2$ [1/V$^2$]') # , fontweight='bold', size=siz); 
    if pH:
        label='Inductance (pH)'
    else:
        label='Inductance (H)'
    ax.set_ylabel(label) #, fontweight='bold', size=siz)
    fig.savefig("figs/06_Ind_Bracale.png",dpi=175,bbox_inches='tight')
    #fig.savefig("101_FigsDambrine/01_UFM2_FG.png",dpi=175,bbox_inches='tight'

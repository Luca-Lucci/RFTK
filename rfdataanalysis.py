from statistics import mean
# import skrf as rf # it is not used in this module
import os   as os
import sys  as sys
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
from rftools import findHz
import smithplot
import matplotlib as mpl

def H21_Ft (Spsweep, bias):
    
    """ The if else below is made in case the user wants to analyse the raw data
    or the de-embedded data.
    Illustration : Not having Y parameters on the dataset most probably means 
    that the data was not de-embedded and you want to see the results
    with the raw data. But having Y parameters does not necessrily mean
    that the S parameters on the dataset are de-embedded. For example,
    the dataset of the 28nm devices from ST has de-embedded Y parameters
    but the S parameters are not. So the else below makes sure that I am
    walways working with de-embedded data for Y  and S parameters. """

    if 'Y' not in Spsweep.keys():
        S=Spsweep['S'][bias]
        Y=rftools.StoY(S)
    else:
        Y=Spsweep['Y'][bias]
        S=rftools.YtoS(Y)     
    F=Spsweep['F']
    Y21=Y[:,1,0]
    Y11=Y[:,0,0]
    Ft, Ftv = rftools.Y4Ptoft(Y, F, 10e9, 30e9)
    currentGain=np.abs(Y21/Y11)
    
    return currentGain, Ft, F

def Mason_Fmax (Spsweep, bias):
    if 'Y' not in Spsweep.keys(): #Refer to the H21_Ft module to understand why this if else was defined.
        S=Spsweep['S'][bias]
        Y=rftools.StoY(S)
    else:
        Y=Spsweep['Y'][bias]
        S=rftools.YtoS(Y)
    
    F=Spsweep['F']

    Fmax, Fmaxv = rftools.Y4Ptofmax(Y, F, 40e9, 50e9)
    masonU=rftools.S2PUMasondB(S)
    
    return masonU, Fmax, F

def IdVg(Spsweep, biases, Lg, Wg) :
    bias_points_inside=len(Spsweep['S'])
    Id=[]
    IdW=[]
    Vg=[]

    for bias in list(range(bias_points_inside)):
        if biases[bias]:
            Id.apend(Spsweep['id'][bias][0])  # we choose 1 point because Id and Vg does not change in function of frequency
            IdW.append((Spsweep['id'][bias][0]*Lg)/Wg)
            Vg.append(Spsweep['vg'][bias])
        
    return Vg, Id, IdW

def estimateVt(IdVg, Vg, W, L, IdTg=300e-9 ):
    Id=np.array(IdVg)
    Idlw=np.log(Id*L/W)
    vg=np.array(Vg)
    f=sp.interpolate.interp1d(Idlw, vg)
    Vth=f( np.log(IdTg) )  
    return float(Vth)

def trivial_extraction(Spsweep, biases,Yplots, Rg=5, Rs=1):  
    
    """ This module is for the extraction of the small signal /
    parameters but without removing the contribution of the /
    Zseries parasitics from the de-embedded S parameters. """
    
    bias_points_inside=len(Spsweep['S'])
    extr=dict()
    extr['Cgs']=[]
    extr['Cgd']=[]
    extr['Fmax_powergain']= []
    extr['Ft_currentgain'] = []
    extr['Ft_complex_RF'] = []
    extr['Fmax_complex_RF'] = []
    extr['Fmax_calc_RF'] = []
    extr['Ft_calc_RF']=[]  # Ft calculated using the formula of Lee or Gilmore
    #extr['Ftgm_Id']=[]
    extr['denomFt']=[]
    extr['Vg']=[]
    extr['Vd']=[]
    extr['Gm']=[]
    extr['Gds']=[]
    extr['gain']=[]
    #extr['Id']=[]
   
    freq=Spsweep['F']

    for bias in list(range(bias_points_inside)):
        if biases[bias]:
            #Idval=Spsweep['id'][bias][0]
            #extr['Id'].append(Idval)  # we choose 1 point because Id and Vg does not change in function of frequency
            if 'VG' in Spsweep.keys():
                extr['Vg'].append(Spsweep['VG'][bias])
                extr['Vd'].append(Spsweep['VD'][bias])
            else :
                extr['Vg'].append(Spsweep['vg'][bias])
                extr['Vd'].append(Spsweep['vd'][bias])
            if 'Y' not in Spsweep.keys(): #Refer to the H21_Ft module to understand why this if else was defined.
                S=Spsweep['S'][bias]
                Y1=rftools.StoY(S)
                Y=substract_Yplots(Y1, Yplots)
            else:
                Y1=Spsweep['Y'][bias]
                Y=substract_Yplots(Y1, Yplots)
                S=rftools.YtoS(Y)

            #Rg, Rgv = rftools.Y2PtoRgDormieuAvg(Y, freq, 20e9, 40e9)
            gm, gmv = rftools.Y2PtoGm(Y, freq, 10e9, 40e9)
            extr['Gm'].append(gm) 
            gd, gdv = rftools.Y2PtoGds(Y, freq, 5e9, 15e9)
            extr['Gds'].append(gd)
            extr['gain'].append(gm/gd)
            cg, cgv = rftools.Y2PtoCgs(Y, freq, 20e9, 40e9)
            extr['Cgs'].append(cg)
            cgd, cgdv = rftools.Y2PtoCgd(Y, freq, 20e9, 40e9)
            extr['Cgd'].append(cgd)
            ft, ftv = rftools.Y4Ptoft(Y, freq, 30e9, 50e9)
            extr['Ft_currentgain'].append(ft*1e9)
            fmax, fmaxv = rftools.Y4Ptofmax(Y, freq, 40e9, 50e9)
            extr['Fmax_powergain'].append(fmax*1e9)
            #extr['Fmax_calc_RF'].append((1/(4*np.pi))*np.sqrt((gm/(cg+cgd))/(Rg*cgd)))
            extr['Fmax_calc_RF'].append((gm/(cg+cgd))/(4*np.pi*np.sqrt((gd*(Rg+Rs))+((gm/(cg+cgd))*Rg*cgd))))
            # Ft_complex_RF.append((gm/(cg*2*np.pi))/((1+(cgd/cg))+(Rs+Rd)*(((cgd/cg)*(gm+gd))+gd)))
            extr['Fmax_complex_RF'].append((gm/(cg*2*np.pi))/(2*(1+(cgd/cg))*np.sqrt((gd*(Rs+Rg))+(0.5*(cgd/cg)*((Rs*gm)+(cgd/cg))))))
            extr["denomFt"].append((cg+cgd)*2*np.pi)
            extr['Ft_calc_RF'].append(gm/((cg+cgd)*2*np.pi))
            #extr["Ftgm_Id"].append((gm*gm)/((cg+cgd)*2*np.pi*Idval))
    return extr

def substract_Zseries (Y, Zserie):
    
    """ This module is defined for the removal of the parasitics /
    from the de-embedded S parameters in order to have the intrinsic /
    elements of the small signal circuit """
    
    if Zserie is not None :
        Z1=rftools.YtoZ_2P(Y)
        Z3=Z1-Zserie
        Y3=rftools.ZtoY_2P(Z3)

    else :
        print ("The Zserie of this device are not known")
        Y3=data[dev][type_param]['Y'][i]
    return Y3
# Danneville didn't mention removing Zseries for the extraction
# of Ft/Fmax frequencies though. He's done it only for the 
# SSEC extraction.
    
def substract_Yplots (Y, Yplots):
    
    """ This module is defined for the removal of the parasitics /
    from the de-embedded S parameters in order to have the intrinsic /
    elements of the small signal circuit """
    
    if Yplots is not None :
        Y3=Y-Yplots

    else :
        print ("The Zserie of this device are not known")
        Y3=data[dev][type_param]['Y'][i]
    return Y3

def Bracale_extraction(Yparams, Yplots, Zserie):
    """ 
    This routine expects as input:
       Yparams[bias,freq,i,j] # bias dependent
       Yplots[freq,i,j]       # Y matrx
       Zseries[freq,i,j]      # Z matrix
       Yparams'=Yparams-Yplots
       Zparams'=inv(Yparams')
       Zparams''=Zparams'-Zseries
       return inv(Zparams'')   
    """
    
 
    Yintr=np.empty(np.shape(Yparams), dtype=np.complex128)
    for bias in range(len(Yparams)):
        Ytmp =Yparams[bias]-Yplots
        Ztmp =rftools.YtoZ_2P(Ytmp)
        Zintr=Ztmp-Zserie
        Yintr[bias]=rftools.ZtoY_2P(Zintr)
    return Yintr

def bracale_extraction(Spsweep, biases, Yplots, Zserie):  #Zserie here is the Zserie of just one device. Ex Zseries["SGL10_D30_C5"]
    
    """ This module is for the extraction of the intrinsic parameters /
    of the small signal circuit after the removal of the contribution /
    of the Zseries parasitics from the de-embedded S parameters. """
    
    bias_points_inside=len(Spsweep['S'])
    extr=dict()
    extr['Cgs']=[]
    extr['Cgd']=[]
    extr['Fmax_powergain']= []
    extr['Ft_currentgain'] = []
    extr['Ft_complex_RF'] = []
    extr['Fmax_complex_RF'] = []
    extr['Fmax_calc_RF'] = []
    extr['Ft_calc_RF']=[]  # Ft calculated using the formula of Lee or Gilmore
    extr['Ft_calc_DC']=[]
    #extr['Ftgm_Id']=[]
    extr['denomFt']=[]
    extr['Vg']=[]
    extr['Vd']=[]
    extr['Gm']=[]
    extr['Gds']=[]
    extr['gain']=[]
    extr['Ri']=[]
    extr['Tau']=[]
    extr['Rgd']=[]
    #extr['Id']=[]
   
    freq=Spsweep['F']

    for bias in list(range(bias_points_inside)):
        if biases[bias]:
            #Idval=Spsweep['id'][bias][0]
            #extr['Id'].append(Idval)  # we choose 1 point because Id and Vg does not change in function of frequency
            if 'VG' in Spsweep.keys():
                extr['Vg'].append(Spsweep['VG'][bias])
                extr['Vd'].append(Spsweep['VD'][bias])
            else :
                extr['Vg'].append(Spsweep['vg'][bias])
                extr['Vd'].append(Spsweep['vd'][bias])
            if 'Y' not in Spsweep.keys(): #Refer to the H21_Ft module to understand why this if else was defined.
                S=Spsweep['S'][bias]
                Y1=rftools.StoY(S)
                Ytot=substract_Yplots(Y1, Yplots)
            else:
                Y1=Spsweep['Y'][bias]
                Ytot=substract_Yplots(Y1, Yplots)
                S=rftools.YtoS(Ytot)
            Y=substract_Zseries(Ytot, Zserie)

            gm, gmv = rftools.Y2PtoGm(Y, freq, 20e9, 40e9)
            extr['Gm'].append(gm)
            Ri, Riv = Y2PtoRi(Y, freq, 25e9, 50e9)
            extr['Ri'].append(Ri)
            Tau, Tauv = Y2Ptotau(Y, freq, 10e9, 35e9)
            extr['Tau'].append(Tau)
            Rgd, Rgdv = Y2PtoRgd(Y, freq, 25e9, 40e9)
            extr['Rgd'].append(Rgd)
            gd, gdv = rftools.Y2PtoGds(Y, freq, 5e9, 15e9)
            extr['Gds'].append(gd)
            extr['gain'].append(gm/gd)
            cg, cgv = rftools.Y2PtoCgs(Y, freq, 20e9, 40e9)
            extr['Cgs'].append(cg)
            cgd, cgdv = rftools.Y2PtoCgd(Y, freq, 20e9, 40e9)
            extr['Cgd'].append(cgd)
            ft, ftv = rftools.Y4Ptoft(Y, freq, 20e9, 40e9)
            extr['Ft_currentgain'].append(ft*1e9)
            fmax, fmaxv = rftools.Y4Ptofmax(Y, freq, 40e9, 50e9)
            extr['Fmax_powergain'].append(fmax*1e9)
            # Fmax_calc_RF.append((1/(4*np.pi))*np.sqrt((gm/(cg+cgd))/(Rg*cgd)))
            # Ft_complex_RF.append((gm/(cg*2*np.pi))/((1+(cgd/cg))+(Rs+Rd)*(((cgd/cg)*(gm+gd))+gd)))
            # Fmax_complex_RF.append((gm/(cg*2*np.pi))/(2*(1+(cgd/cg))*np.sqrt((gd*(Rs+Rg))+(0.5*(cgd/cg)*((Rs*gm)+(cgd/cg))))))
            extr["denomFt"].append((cg+cgd)*2*np.pi)
            extr['Ft_calc_RF'].append(gm/((cg+cgd)*2*np.pi))
            #extr["Ftgm_Id"].append((gm*gm)/((cg+cgd)*2*np.pi*Idval))
    return extr
    
def resistances_values (Denom, Rg, Rd, Rs):
    slopeRg, interceptRg, r_valueRg, p_valueRg, std_errRg = stats.linregress(Denom,Rg)
    #print ("Rg = ", interceptRg, " ohms")
    slopeRd, interceptRd, r_valueRd, p_valueRd, std_errRd = stats.linregress(Denom,Rd)
    #print ("Rd = ", interceptRd, " ohms")
    slopeRs, interceptRs, r_valueRs, p_valueRs, std_errRs = stats.linregress(Denom,Rs)
    #print ("Rs = ", interceptRs, " ohms")
    return interceptRg, interceptRd, interceptRs

def extract_RgVg (Spsweep,  biases, Yplots):
    
    """ This module has been defined in order to verify the /
    values of the gate resistance extracted using Bracale's /
    method. """
    
    bias_points_inside=len(Spsweep['Y'])
    Rg_Enz=[]
    Rg_Jen=[]
    Rg_Dormieu=[]
    Vgs=[]

    for bias in list(range(bias_points_inside)):
        if biases[bias]:
            Vg=Spsweep['vg'][bias]
            Vgs.append(Vg)      
            if 'Y' not in Spsweep.keys(): #Refer to the H21_Ft module to understand why this if else was defined.
                S=Spsweep['S'][bias]
                Y1=rftools.StoY(S)
                Y=substract_Yplots(Y1, Yplots)
            else:
                Y1=Spsweep['Y'][bias]
                Y=substract_Yplots(Y1, Yplots)
                S=rftools.YtoS(Y)
            F=Spsweep['F']

            Rg1, Rg1v = rftools.Y2PtoRgEnzAvg(Y, F, 20e9, 40e9)
            Rg_Enz.append(Rg1)
            Rg2, Rg2v = rftools.Y2PtoRgJenAvg(Y, F, 20e9, 40e9)
            Rg_Jen.append(Rg2)
            Rg3, Rg3v = rftools.Y2PtoRgDormieuAvg(Y, F, 20e9, 40e9)
            Rg_Dormieu.append(Rg3)
        
    return Vgs, Rg_Enz, Rg_Jen, Rg_Dormieu

def extract_Rgfreq (Spsweep, bias, Yplots):
    
    """ It does the same as the extract_RgVg module but as a function /
    of frequency in stead of gate voltage (Vg)"""
    
    Rg_Enz=[]
    Rg_Jen=[]
    Rg_Dormieu=[]
    
    if 'Y' not in Spsweep.keys(): #Refer to the H21_Ft module to understand why this if else was defined.
        S=Spsweep['S'][bias]
        Y1=rftools.StoY(S)
        Y=substract_Yplots(Y1, Yplots)
    else:
        Y1=Spsweep['Y'][bias]
        Y=substract_Yplots(Y1, Yplots)
        S=rftools.YtoS(Y)
    F=Spsweep['F']
    Rg_Enz, Rg_Enzv = rftools.Y2PtoRgEnz(Y, F)
    Rg_Jen, Rg_Jenv = rftools.Y2PtoRgJen(Y, F)
    Rg_Dormieu, Rg_Dormieuv = rftools.Y2PtoRgDormieu(Y, F)
        
    return F, Rg_Enz, Rg_Jen, Rg_Dormieu

def extract_seriesresistances (Spsweep, biases, Yplots, Vth, Vgoverdrive=0.2, fstart=20e9, fstop=40e9):
    Rgextract=[]
    Rdextract=[]
    Rsextract=[]
    Vgs=[]
    Denom=[]

    for bias,go in enumerate(biases):
        if go:
            if 'vg' in Spsweep.keys():
                Vg=Spsweep['vg'][bias]
            else:
                Vg=Spsweep['VG'][bias]
            if Vg < 0 :   # This is for P_type devices
                Vg=Vg*-1
            Vgs.append(Vg)     
            if 'Y' not in Spsweep.keys(): #Refer to the H21_Ft module to understand why this if else was defined.
                S=Spsweep['S'][bias]
                Y1=rftools.StoY(S)
                Y=substract_Yplots(Y1, Yplots)
            else:
                Y1=Spsweep['Y'][bias]
                Y=substract_Yplots(Y1, Yplots)
                S=rftools.YtoS(Y)
            F=Spsweep['F']
        
            if  Vg-Vth >Vgoverdrive :
                Denom.append(1/(Vg-Vth))
                Rg, Rgv = rftools.Y2PtoRgBracale(Y, F, fstart, fstop)
                Rd, Rdv = rftools.Y2PtoRd(Y, F, fstart, fstop)
                Rs, Rsv = rftools.Y2PtoRs(Y, F, fstart, fstop)
                Rgextract.append(Rg)
                Rdextract.append(Rd)
                Rsextract.append(Rs)
            
    return Denom, Rgextract, Rdextract, Rsextract

def extract_seriesinductances (Spsweep, biases, Yplots, Vth, Vgoverdrive=0.2, fstart=20e9, fstop=40e9, fstartlg=25e9, fstoplg=1000e9):
    
    #bias_points_inside=len(Spsweep['S'])
    Lgextract=[]
    Ldextract=[]
    Lsextract=[]
    Vgs=[]
    Denom=[]

    for bias,go in enumerate(biases):
        if go:
            if 'vg' in Spsweep.keys():
                Vg=Spsweep['vg'][bias]
            else :
                Vg=Spsweep['VG'][bias]
            if Vg < 0 :   # This is for P_type devices
                Vg=Vg*-1
            Vgs.append(Vg)      
            if 'Y' not in Spsweep.keys(): #Refer to the H21_Ft module to understand why this if else was defined.
                S=Spsweep['S'][bias]
                Y1=rftools.StoY(S)
                Y=substract_Yplots(Y1, Yplots)
            else:
                Y1=Spsweep['Y'][bias]
                Y=substract_Yplots(Y1, Yplots)
                S=rftools.YtoS(Y)
            #Y=substract_Zseries(Ytot, Zserie)
            F=Spsweep['F']
 
            startlg, stoplg=findHz(F, fstartlg, fstoplg)
            if  Vg-Vth >Vgoverdrive :
                Denom.append(1/((Vg-Vth)**2))
                Lg, Lgv, Fl = rftools.Y2PtoLg(Y, F)
                slopeLg1, interceptLg1, r_valueLg1, p_valueLg1, std_errLg1 = stats.linregress(Fl[startlg:stoplg],Lgv[startlg:stoplg])
                Ld, Ldv = rftools.Y2PtoLdAvg(Y, F, fstart, fstop)
                Ls, Lsv = rftools.Y2PtoLsAvg(Y, F, fstart, fstop)
                Lgextract.append(interceptLg1)
                Ldextract.append(Ld)
                Lsextract.append(Ls)
            
    return Denom, Lgextract, Ldextract, Lsextract

def inductances_values (Denom, Lg, Ld, Ls):
    slopeLg, interceptLg, r_valueLg, p_valueLg, std_errLg = stats.linregress(Denom,Lg)
    #print ("Lg = ", interceptLg, " H")
    slopeLd, interceptLd, r_valueLd, p_valueLd, std_errLd = stats.linregress(Denom,Ld)
    #print ("Ld = ", interceptLd, " H")
    slopeLs, interceptLs, r_valueLs, p_valueLs, std_errLs = stats.linregress(Denom,Ls)
    #print ("Ls = ", interceptLs, " H")
    return interceptLg, interceptLd, interceptLs

def extract_OpenShort (Yopen, Yshort, freq):
    
    """ This module is for the extraction of the Cpg, Cpd, Ccoupling..."""
    extr={}
    extr["Cop1"]=[]
    extr["Cop2"]=[]
    extr["C12"]=[]
    extr["C21"]=[]
    extr["Rs1"]=[]
    extr["Rs2"]=[]
    
    extr["Cop1"], Copv = rftools.Y2PtoCopen1(Yopen, freq)
    extr["Cop2"], Copv2 = rftools.Y2PtoCopen2(Yopen, freq)
    extr["C12"], C12v = rftools.Y2PtoC12(Yopen, freq)
    extr["C21"], C21v = rftools.Y2PtoC21(Yopen, freq)
    Ysh= Yshort-Yopen
    extr["Rs1"], Rs1v = rftools.Y2PtoRs1(Ysh, freq)
    extr["Rs2"], Rs2v = rftools.Y2PtoRs2(Ysh, freq)
    
    return extr

def extract_NFminAvg (NFsweep, bias, fstart=10, fstop=1000e9):
    fHz=NFsweep['F']
    start, stop = findHz(fHz, fstart, fstop)
    NFminv=NFsweep['v1_NFmin'][bias]
    NFmin=sum(NFminv[start:stop])/len(NFminv[start:stop])
    return NFmin, NFminv

def extract_RnAvg (NFsweep, bias, fstart=10, fstop=1000e9):
    fHz=NFsweep['F']
    start, stop = findHz(fHz, fstart, fstop)
    Rnv=NFsweep['v1_Rn'][bias]
    Rn=sum(Rnv[start:stop])/len(Rnv[start:stop])
    return Rn, Rnv

def Y2PtoRi (Y, fHz, fstart=10, fstop=1000e9):
    """ This is for the extraction of the channel or intrinsic resistance"""
    start, stop = findHz(fHz, fstart, fstop)
    Xv=np.real(1/(Y[:,0,0]+Y[:,0,1]))
    #Xv=np.real(Y[:,0,0]) / np.sqr(np.imag(Y[:,O,0]+Y[:,0,1]))
    X=sum(Xv[start:stop])/len(Xv[start:stop])
    return X, Xv

def Y2PtoRgd (Y, fHz, fstart=10, fstop=1000e9):
    """ This is for the extraction of the channel or intrinsic resistance"""
    start, stop = findHz(fHz, fstart, fstop)
    Xv=-np.real(1/(Y[:,0,1]))
    X=sum(Xv[start:stop])/len(Xv[start:stop])
    return X, Xv

def Y2Ptotau (Y, fHz, fstart=10, fstop=1000e9):
    """This module is for the determination of the transit time (tau) when the series 
    admittances are not removed using Bracale. So this applies directly
    to the de-embedded data"""
    start, stop = findHz(fHz, fstart, fstop)
    Xv=(-1/(fHz*2*np.pi))*(np.angle((Y[:,1,0]-Y[:,0,1])/(Y[:,0,0]+Y[:,0,1]), deg=False)+(np.pi/2))
    #The following line contains a mistake. I used Y[:,1,1] instead of Y[:,0,0].
    #Xv=(-1/(fHz*2*np.pi))*(np.angle((Y[:,1,0]-Y[:,0,1])/(Y[:,1,1]+Y[:,0,1]), deg=False)+(np.pi/2))
    X=sum(Xv[start:stop])/len(Xv[start:stop])
    return X, Xv

def Y2PtotauBracale (Y, fHz, fstart=10, fstop=1000e9):
    """This module is for the determination of the transit time (tau) when the series 
    admittances are removed using Bracale. So this applies directly
    to the de-embedded data"""
    start, stop = findHz(fHz, fstart, fstop)
    Xv=(-1/(fHz*2*np.pi))*np.angle(Y[:,1,0]-Y[:,0,1], deg=False)
    # np.angle receives Z param not Y
    X=sum(Xv[start:stop])/len(Xv[start:stop])
    return X, Xv

def Y2PtoRgChen(Y, fHz, fstart=10, fstop=1000e9):
    """ This extraction uses the new formula proposed by Chen (McMaster Uni) in his thesis
    NOISE CHARACTERIZATION AND MODELING OF NANOSCALE MOSFETS (2017)"""
    start, stop = findHz(fHz, fstart, fstop)
    ZZ=rftools.YtoZ_2P(Y)
    Xv=np.real(ZZ[:,0,0])-np.real(ZZ[:,1,1]/4)
    X=sum(Xv[start:stop])/len(Xv[start:stop])
    return X, Xv

def Y2PtoCcal1(Y, fHz, fstart=10, fstop=1000e9):
    start, stop = findHz(fHz, fstart, fstop)
    Copenv=(np.imag(Y[:,0,0]+Y[:,0,1]))/(fHz*2*np.pi)
    Copen=sum(Copenv[start:stop])/len(Copenv[start:stop])
    return Copen, Copenv

def Y2PtoCcal2(Y, fHz, fstart=10, fstop=1000e9):
    start, stop = findHz(fHz, fstart, fstop)
    Copenv=(np.imag(Y[:,1,1]+Y[:,0,1]))/(fHz*2*np.pi)
    Copen=sum(Copenv[start:stop])/len(Copenv[start:stop])
    return Copen, Copenv

def Y2PtoLcal1(Y, fHz, fstart=10, fstop=1000e9):
    start, stop = findHz(fHz, fstart, fstop)
    Lshortv=np.imag(1/(Y[:,0,0]+Y[:,0,1]))/(fHz*2*np.pi)
    Lshort=sum(Lshortv[start:stop])/len(Lshortv[start:stop])
    return Lshort, Lshortv

def Y2PtoLcal2(Y, fHz, fstart=10, fstop=1000e9):
    start, stop = findHz(fHz, fstart, fstop)
    Lshortv=np.imag(1/(Y[:,1,1]+Y[:,0,1]))/(fHz*2*np.pi)
    Lshort=sum(Lshortv[start:stop])/len(Lshortv[start:stop])
    return Lshort, Lshortv

def Y2PtoRcal1 (Y, fHz, fstart=10, fstop=1000e9):
    start, stop = findHz(fHz, fstart, fstop)
    Xv=np.real(1/(Y[:,0,0]+Y[:,0,1]))
    X=sum(Xv[start:stop])/len(Xv[start:stop])
    return X, Xv

def Y2PtoRcal2 (Y, fHz, fstart=10, fstop=1000e9):
    start, stop = findHz(fHz, fstart, fstop)
    Xv=np.real(1/(Y[:,1,1]+Y[:,0,1]))
    X=sum(Xv[start:stop])/len(Xv[start:stop])
    return X, Xv

def Y2Ptodelaytime (S, fHz, fstart=10, fstop=1000e9):
    start, stop = findHz(fHz, fstart, fstop)
    Xv=(-1/(fHz*2*np.pi))*np.angle(S[:,1,0], deg=False)
    X=sum(Xv[start:stop])/len(Xv[start:stop])
    return X, Xv

def Ytofmax2(Y, fHz, fstart=10, fstop=1000e9):
    start, stop = findHz(fHz, fstart, fstop)
    S=rftools.YtoS(Y)
    num=S[:,1,0]
    den=S[:,0,1]
    fmaxv=np.abs(num/den)  
    fmax=sum(fmaxv[start:stop])/1e9/len(fmaxv[start:stop])
    return fmax, fmaxv

def reject_outliers_2(data, freq, m = 5.):
    # eliminates undesirable points from the noise data
    d = np.abs(data - np.median(data))
    mdev = np.median(d)
    s = d/(mdev if mdev else 1.)
    return data[s<m], freq[s<m]

def getSparamfromS2Pfile (files):
    f=files.readlines()
    for counter, value in enumerate(f): # This loop is to get the frequency range
        vect=value.strip().split(' ')
        tupl=vect
        while '' in tupl:
            tupl.remove('') 
        tupl2=[]
        for i in tupl:
            tupl2.append(i)
        if tupl2==['!']: #this list corresponds to the line after which the first frequency
                         # starts.
            Sstart=counter
        if tupl2==['!', 'Noise', 'params']:  #this list corresponds to the line before which
                                             #the last frequency is.
            Sstop=counter
            
    numOfFreqPoints=Sstop-Sstart-1
    
    sm=np.zeros((numOfFreqPoints,2,2),dtype=np.complex128)
    
    for counter, value in enumerate(f):
        if counter>Sstart  and counter<Sstop:
            vect=value.strip().split(' ')
            tupl=vect
            while '' in tupl:
                tupl.remove('')  
            tupl2=[]
            for i in tupl:
                tupl2.append(float(i))
            S11=tupl2[1]+tupl2[2]*1j
            S21=tupl2[3]+tupl2[4]*1j
            S12=tupl2[5]+tupl2[6]*1j
            S22=tupl2[7]+tupl2[8]*1j
            freqpoint=counter-5
            sm[freqpoint,0,0]=S11
            sm[freqpoint,1,0]=S21
            sm[freqpoint,0,1]=S12
            sm[freqpoint,1,1]=S22
    return sm

def getSparamfromS2P (file):
    f=file.readlines()
    Sstop=np.size(f)
    test=False
    for counter, value in enumerate(f): # This loop is to get the frequency range
        vect=value.strip().split(' ')
        tupl=vect
        while '' in tupl:
            tupl.remove('') 
        tupl2=[]
        for i in tupl:
            tupl2.append(i)
        mode=str.isnumeric((tupl2[0][0]))
        if mode==True and test==False:
            Sstart=counter
            test=True
    numOfFreqPoints=Sstop-Sstart
    freq=np.zeros((numOfFreqPoints),dtype=np.float_)
    sm=np.zeros((numOfFreqPoints,2,2),dtype=np.complex_)
    for counter, value in enumerate(f):
        if counter>=Sstart  and counter<=Sstop:
            vect=value.strip().split(' ')
            tupl=vect
            while '' in tupl:
                tupl.remove('')  
            tupl2=[]
            for i in tupl:
                tupl2.append(float(i))
            F=tupl2[0]
            S11=tupl2[1]+tupl2[2]*1j
            S12=tupl2[3]+tupl2[4]*1j
            S21=tupl2[5]+tupl2[6]*1j
            S22=tupl2[7]+tupl2[8]*1j
            freqpoint=counter-Sstart
            sm[freqpoint,0,0]=S11
            sm[freqpoint,1,0]=S21
            sm[freqpoint,0,1]=S12
            sm[freqpoint,1,1]=S22 
            freq[freqpoint]=F
    return freq, sm

def getIcorfromS2P (file):
    f=file.readlines()
    Sstop=np.size(f)
    test=False
    for counter, value in enumerate(f): # This loop is to get the frequency range
        vect=value.strip().split(' ')
        tupl=vect
        while '' in tupl:
            tupl.remove('') 
        tupl2=[]
        for i in tupl:
            tupl2.append(i)
        mode=str.isnumeric((tupl2[0][0]))
        if mode==True and test==False:
            Sstart=counter
            test=True
    numOfFreqPoints=Sstop-Sstart
    freq=np.zeros((numOfFreqPoints),dtype=np.complex128)
    icor=np.zeros((numOfFreqPoints,2,2),dtype=np.complex128)
    for counter, value in enumerate(f):
        if counter>=Sstart  and counter<=Sstop:
            vect=value.strip().split(' ')
            tupl=vect
            while '' in tupl:
                tupl.remove('')  
            tupl2=[]
            for i in tupl:
                tupl2.append(float(i))
            F=tupl2[0]
            icor11=tupl2[1]+tupl2[2]*1j
            icor12=tupl2[3]+tupl2[4]*1j
            icor21=tupl2[5]+tupl2[6]*1j
            icor22=tupl2[7]+tupl2[8]*1j
            freqpoint=counter-Sstart
            icor[freqpoint,0,0]=icor11
            icor[freqpoint,1,0]=icor21
            icor[freqpoint,0,1]=icor12
            icor[freqpoint,1,1]=icor22 
            freq[freqpoint]=F
    return freq, icor

def getSoptfromADS(files):
    f=files.readlines()
    for counter, value in enumerate(f): # This loop is to get the frequency range
        vect=value.strip().split(' ')
        tupl=vect
        while '' in tupl:
            tupl.remove('')
        #print (tupl[0])
        tupl2=[]
        for i in tupl:
            tupl2.append(i)
        #print (tupl2)
        if tupl2==['freq\tSopt']:
            Sstart=counter
            #print (Sstart)
        if tupl2==[]: 
            Sstop=counter
            #print (Sstop)
            numOfFreqPoints=Sstop-Sstart-1
            #print (numOfFreqPoints)
            break    

    sm=np.zeros(numOfFreqPoints, dtype=np.complex128)
    sm2=np.zeros(numOfFreqPoints, dtype=np.complex128)
    i=0
    for counter, value in enumerate(f): # This loop is to get the frequency range
        if counter>Sstart  and counter<Sstop:
            vect=value.strip().split(' ')
            tupl=vect
            while '' in tupl:
                tupl.remove('')
            tupl2=value.split('\t', 2)
            tupl3=tupl2[-1].split('/', 2)
            tupl4=tupl3[-1].split('\n')
            #tupl3=value.split('/')
            #print (tupl2[0], tupl3[0], tupl4[0])
            Sopt=float(tupl3[0])+float(tupl4[0])*1j
            sm[i]=Sopt
            i+=1
    #conversion from angles and absolutes to a complex number
    sm2=np.array([np.complex(np.real(a)*np.cos(math.radians(np.imag(a))),np.real(a)*np.sin(math.radians(np.imag(a)))) for a in sm])
            
    return sm2

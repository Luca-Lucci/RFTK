# -*- coding: latin1 -*-
# ver. 0.7 2023-09 bug resolved in SS_STM
# ver. 0.6 2023-07 update graphs to py 3.91 (delete b=True from grid)
#                  allow SnPparse from list (i.e. zipfiles)
# ver. 0.5 2023-03 improved export of Touchstone files
# ver. 0.4 2023-03 fixed SnPparse that was ignoring GHz and onther f-specs
# ver. 0.3 added some properties to SStrivial
# ver. 0.2
# added deftable to SS_Trivial to make it easier to mode SS definitions from 
# devices with different f sweeps

import subprocess
import uW 
import rfdataanalysis
import rftools as rft
import os
import numpy as np
import time
import re as re
import matplotlib.pyplot as plt

windows=True

class ADSsim :
   status = 0 # init has been performed
   # path to run ADS into (save logs and data without explicit path redir
   PathADS="C:/ADStmp/"
   # for the time being this is hardcoded, may relax this one day ...
   version="2020.00"
   version="2021"
   version="2021myrev210507"
   fHz = None
   # will store environment variables to run ADS
   my_env= None
   netlist_Hparam_old="""
; old netlist with Ri and Cgs invertite
; Top Design: "SSC_lib:MW2216_s4d4w2_PSO:sgslvtn_s4d4w2"
; Netlisted using Hierarchy Policy: "Standard"

Options ResourceUsage=yes UseNutmegFormat=no EnableOptim=no TopDesignName="SSC_lib:MW2216_s4d4w2_PSO:sgslvtn_s4d4w2_Vg=0.8V" DcopOutputNodeVoltages=yes DcopOutputPinCurrents=yes DcopOutputAllSweepPoints=no DcopOutputDcopType=0
#uselib "ckt" , "VCCS"
S_Param:SP1 CalcS=yes CalcY=yes CalcZ=no GroupDelayAperture=1e-4 FreqConversion=no FreqConversionPort=1 StatusLevel=2 CalcNoise=yes SortNoise=0 BandwidthForNoise=100 Hz DevOpPtLevel=0 \\
SweepVar="freq" SweepPlan="SP1_stim" OutputPlan="SP1_Output" 

;; SweepPlan: SP1_stim Start=100 MHz Stop=50 GHz Step=100 MHz 

SweepPlan: SP1_stim {SWEEP}

OutputPlan:SP1_Output \\
      Type="Output" \\
      UseEquationNestLevel=yes \\
      EquationNestLevel=2 \\
      UseSavedEquationNestLevel=yes \\
      SavedEquationNestLevel=2

#load "python","LinearCollapse"
Component Module="LinearCollapse" Type="ModelExtractor" NetworkRepresentation=2
C:Cgs  N__15 N__4 C={Cgs} F
R:Ri  N__9 N__15 R={Ri} Ohm Noise=no 
R:Gds  N__8 N__4 R={Gds} Ohm Noise=no 
C:Cgd  N__9 N__16 C={Cgd} F
R:Rgd  N__16 N__8 R={Rgd} Ohm Noise=no 
VCCS:Gm  N__9 N__4 N__8 N__4 G={Gm} S R1=1e100 Ohm R2=1e100 Ohm F=0.0 GHz 

Options:Options1 Temp=16.85 Tnom=16.85 TopologyCheck=yes ForceS_Params=yes GiveAllWarnings=yes  \\
MaxWarnings=10 ForceM_Params=yes InitialGuessAnnotation=0 TopologyCheckMessages=no doDeltaAC=no  \\
ReduceSPortRatio=0.5 WarnSOA=yes MaxWarnSOA=5 Census=no MinNodalR=1e-2 Ohm  \\
LargeCapThreshold=1e-3 F EnableAssertChecks=no EnableSpareRemoval=yes 

Port:Term1  N__11 0 Num=1 Z=50 Ohm Noise=yes Temp=27 
Port:Term2  N__3 0 Num=2 Z=50 Ohm Noise=yes Temp=27 
R:Rd  N__8 N__3 R={Rd} Ohm Noise=no 
R:Rg  N__11 N__9 R={Rg} Ohm Noise=no 
R:Rs  N__4 0 R={Rs} Ohm Noise=no 
aele H=stoh(S);


aele x=write_snp("{dev}_Sparam.s2p",S,"S-par simulation data","Hz","MA",50,50, 1, 9, "  ");

aele temporary2=write_var("{dev}_Sparam2.s2p", "W", "## S11.r S11.i S12.r S12.i S21.r S21.i S22.r S22.i ", " ", "s", 10, \\
      freq, S(1,1) , S(1,2) ,S(2,1) ,S(2,2)  ) ;
"""


   netlist_CpgCpd="""
; Top Design: "SpiceTools_lib:DevCpdCpg:schematic"

Options ResourceUsage=yes UseNutmegFormat=no EnableOptim=no TopDesignName="SpiceTools_lib:DevCpdCpg:schematic" \\
        DcopOutputNodeVoltages=yes DcopOutputPinCurrents=yes DcopOutputAllSweepPoints=no DcopOutputDcopType=0

S_Param:SP1 CalcS=yes CalcY=yes CalcZ=no GroupDelayAperture=1e-4 FreqConversion=no FreqConversionPort=1 StatusLevel=2 CalcNoise=yes SortNoise=0 BandwidthForNoise=100 Hz DevOpPtLevel=0 \\
SweepVar="freq" SweepPlan="SP1_stim" OutputPlan="SP1_Output" 

SweepPlan: SP1_stim {SWEEP}

OutputPlan:SP1_Output \\
      Type="Output" \\
      UseEquationNestLevel=yes \\
      EquationNestLevel=2 \\
      UseSavedEquationNestLevel=yes \\
      SavedEquationNestLevel=2

Port:TermG1  Nin  0 Num=1 Z=50 Ohm Noise=yes 
Port:TermG2  Nout 0 Num=2 Z=50 Ohm Noise=yes 
R:Ri   Nintgs Nsource R={Ri}  Noise=no 
R:RGD  Nintgd Ndrain  R={Rgd} Noise=no
C:CGD  Ngate  Nintgd  C={Cgd} F
C:CGS  Ngate  Nintgs  C={Cgs} F 
R:Rds  Ndrain Nsource R={Rds} Noise=no
#uselib "ckt" , "VCCS"
VCCS:GM  Ngate Nsource Ndrain Nsource G={Gm} S T=0.0 nsec R1=1e100 Ohm R2=1e100 Ohm F=0.0 GHz 

C:Cpg  Nin       0     C={Cpg} F
C:Cpd  Nout      0     C={Cpd} F 
L:Lg   Nin       Ngate L={Lg} R={Rg} Noise=no  
L:Ls   Nsource   0     L={Ls} R={Rs} Noise=no
L:Ld   Ndrain    Nout  L={Ld} R={Rd} Noise=no 

Options:Options1 Temp=16.85 Tnom=16.85 TopologyCheck=yes ForceS_Params=yes GiveAllWarnings=yes  \\
MaxWarnings=10 ForceM_Params=yes InitialGuessAnnotation=0 TopologyCheckMessages=no doDeltaAC=no  \\
ReduceSPortRatio=0.5 WarnSOA=yes MaxWarnSOA=5 Census=no MinNodalR=1e-2 Ohm  \\
LargeCapThreshold=1e-3 F EnableAssertChecks=no EnableSpareRemoval=yes 

aele H=stoh(S);
aele x=write_snp("{dev}_CpgCpd_Sparam.s2p",S,"S-par simulation data","Hz","MA",50,50, 1, 9, "  ");
aele temporary2=write_var("{dev}_CpgCpd_Sparam2.s2p", "W", "## S11.r S11.i S12.r S12.i S21.r S21.i S22.r S22.i ", " ", "s", 10, \\
      freq, S(1,1) , S(1,2) ,S(2,1) ,S(2,2)  ) ;
"""
   netlist_Hparam="""
; Top Design: "SSC_lib:MW2216_s4d4w2_PSO:sgslvtn_s4d4w2"
; Netlisted using Hierarchy Policy: "Standard"

Options ResourceUsage=yes UseNutmegFormat=no EnableOptim=no TopDesignName="SSC_lib:MW2216_s4d4w2_PSO:sgslvtn_s4d4w2_Vg=0.8V" \\
        DcopOutputNodeVoltages=yes DcopOutputPinCurrents=yes DcopOutputAllSweepPoints=no DcopOutputDcopType=0
#uselib "ckt" , "VCCS"
S_Param:SP1 CalcS=yes CalcY=yes CalcZ=no GroupDelayAperture=1e-4 FreqConversion=no FreqConversionPort=1 StatusLevel=2 CalcNoise=yes SortNoise=0 BandwidthForNoise=100 Hz DevOpPtLevel=0 \\
SweepVar="freq" SweepPlan="SP1_stim" OutputPlan="SP1_Output" 

;; SweepPlan: SP1_stim Start=100 MHz Stop=50 GHz Step=100 MHz 

SweepPlan: SP1_stim {SWEEP}

OutputPlan:SP1_Output \\
      Type="Output" \\
      UseEquationNestLevel=yes \\
      EquationNestLevel=2 \\
      UseSavedEquationNestLevel=yes \\
      SavedEquationNestLevel=2

#load "python","LinearCollapse"
;; Component Module="LinearCollapse" Type="ModelExtractor" NetworkRepresentation=2
C:Cgs  N__9 N__15 C={Cgs} F
R:Ri  N__15 N__4 R={Ri} Ohm Noise=no 
R:Gds  N__8 N__4 R={Gds} Ohm Noise=no 
C:Cgd  N__9 N__16 C={Cgd} F
R:Rgd  N__16 N__8 R={Rgd} Ohm Noise=no 
VCCS:Gm  N__9 N__4 N__8 N__4 G={Gm} S R1=1e100 Ohm R2=1e100 Ohm F=0.0 GHz 

Options:Options1 Temp=16.85 Tnom=16.85 TopologyCheck=yes ForceS_Params=yes GiveAllWarnings=yes  \\
MaxWarnings=10 ForceM_Params=yes InitialGuessAnnotation=0 TopologyCheckMessages=no doDeltaAC=no  \\
ReduceSPortRatio=0.5 WarnSOA=yes MaxWarnSOA=5 Census=no MinNodalR=1e-2 Ohm  \\
LargeCapThreshold=1e-3 F EnableAssertChecks=no EnableSpareRemoval=yes 

Port:Term1  N__11 0 Num=1 Z=50 Ohm Noise=yes Temp=27 
Port:Term2  N__3 0 Num=2 Z=50 Ohm Noise=yes Temp=27 
R:Rd  N__8 N__3 R={Rd} Ohm Noise=no 
R:Rg  N__11 N__9 R={Rg} Ohm Noise=no 
R:Rs  N__4 0 R={Rs} Ohm Noise=no 

aele H=stoh(S);
aele x=write_snp("{dev}_Sparam.s2p",S,"S-par simulation data","Hz","MA",50,50, 1, 9, "  ");
aele temporary2=write_var("{dev}_Sparam2.s2p", "W", "## S11.r S11.i S12.r S12.i S21.r S21.i S22.r S22.i ", " ", "s", 10, \\
      freq, S(1,1) , S(1,2) ,S(2,1) ,S(2,2)  ) ;
"""

   netlist_NoisePSO="""
; Top Design: "SpiceTools_lib:NoiseOS_vs_genericMDIF:schematic_PSO" 
; Netlist for PSO simultations with noise

Options ResourceUsage=yes UseNutmegFormat=no EnableOptim=no \\
        TopDesignName="SpiceTools_lib:NoiseOS_vs_genericMDIF:schematic_PSO" DcopOutputNodeVoltages=yes DcopOutputPinCurrents=yes\\
        DcopOutputAllSweepPoints=no DcopOutputDcopType=0

Port:TermG1  Port1 0 Num=1 Z=50 Ohm Noise=yes Temp=27 
Port:TermG2  Port2 0 Num=2 Z=50 Ohm Noise=no Temp=27 

S_Param:SP1 CalcS=yes CalcY=yes CalcZ=no GroupDelayAperture=1e-4 FreqConversion=no FreqConversionPort=1 StatusLevel=2 CalcNoise=yes SortNoise=0 BandwidthForNoise={DF} Hz DevOpPtLevel=0 NoiseInputPort=1 NoiseOutputPort=2 \\
SweepVar="freq" SweepPlan="SP1_stim" OutputPlan="SP1_Output" 

SweepPlan: SP1_stim {SWEEP}

OutputPlan:SP1_Output \\
      Type="Output" \\
      UseEquationNestLevel=yes \\
      EquationNestLevel=2 \\
      UseSavedEquationNestLevel=yes \\
      SavedEquationNestLevel=2 

C:Cpad2   Port2 0     C={Cpad2}
C:Cpad1   Port1 0     C={Cpad1}
C:Cpad12  Port2 Port1 C={Cpad12}

L:Lshort1  Port1          Nin   L={Lsh1} R={Rsh1} Noise=yes  
L:Lshort2  Nout           Port2 L={Lsh2} R={Rsh2} Noise=yes  
L:Lshort3  GNDbeforeShort 0     L={Lsh3} R={Rsh3} Noise=yes  

C:Copen1   Nin  GNDbeforeShort C={Copen1}
C:Copen2   Nout GNDbeforeShort C={Copen2}
C:Copen12  Nin  Nout           C={Copen12}

C:Cpg  Nin  GNDbeforeShort C={Cpg}
C:Cpd  Nout GNDbeforeShort C={Cpd}

L:Lg  Nin        Nint           L={Lg}  Noise=no  
R:Rg  Nintgate   Ngate          R={Rg}  Noise=no 
L:Ld  Noutdrain  Nout           L={Ld}  Noise=no  
R:Rd  Ndrain     Noutdrain      R={Rd}  Noise=no 
L:Ls  Nintsource GNDbeforeShort L={Ls}  Noise=no  
R:Rs  Nsource    Nintsource     R={Rs}  Noise=no 

R:Ri   Nintgs Nsource R={Ri}  Noise=no 
R:RGD  Nintgd Ndrain  R={Rgd} Noise=no 
C:CGD  Ngate  Nintgd  C={Cgd}
C:CGS  Ngate  Nintgs  C={Cgs}
R:Rds  Ndrain Nsource R={Rds} Noise=no 
#uselib "ckt" , "VCCS"
VCCS:GM  Ngate Nsource Ndrain Nsource G={Gm} T=0.0 nsec R1=1e100 Ohm R2=1e100 Ohm F=0.0 GHz 

kb=1.38e-23
DF={DF}
Tout={Tout}
Tin=290
H22={H22}
H11={H11}
vg_noise=sqrt(4*kb*Tin*DF*H11)
id_noise=sqrt(4*kb*Tout*DF*H22)

V_Source:noiseVg  Nint Nintgate Type="V_Noise" V_Noise=vg_noise SaveCurrent=1 
NoiseCorr2Port:SRC3  CorrCoeff=0.0 Source1="noiseVg" Source2="noiseId" 
I_Source:noiseId  Noutdrain Nintsource Type="I_Noise" I_Noise=id_noise 

Options:Options1 Temp=16.85 Tnom=25 TopologyCheck=yes ForceS_Params=yes V_RelTol=1e-6  \\
I_RelTol=1e-6 GiveAllWarnings=yes MaxWarnings=10 ForceM_Params=yes InitialGuessAnnotation=0  \\
TopologyCheckMessages=no doDeltaAC=yes ReduceSPortRatio=0.5 WarnSOA=yes MaxWarnSOA=5  \\
Census=no MinNodalR=1e-2 Ohm LargeCapThreshold=1e-3 F EnableAssertChecks=no EnableSpareRemoval=yes 


aele x=write_snp("{dev}_PSO_Sparam.s2p",S,"S-par simulation data","Hz","MA",50,50, 1, 9, "  ");
aele temporary1=write_var("{dev}_PSO_NFdata.dat", "W", "## freq, NFmin, Rn, Sopt, nf(1), nf(2), te(1), te(2) ", " ", "s", 10, \\
      freq, NFmin, Rn, Sopt, nf(1), nf(2), te(1), te(2) );
aele temporary2=write_var("{dev}_PSO_Sparam2.dat", "W", "## S11.r S11.i S12.r S12.i S21.r S21.i S22.r S22.i ", " ", "s", 10, \\
      freq, S(1,1) , S(1,2) ,S(2,1) ,S(2,2)  ) ;
aele temporary3=write_var("{dev}_PSO_Cydata.dat", "W", "## icor11.r icor11.i icor12.r icor12.i icor21.r icor21.i icor22.r icor22.i ", " ", "s", 10, \\
      freq, icor(1,1) , icor(1,2) ,icor(2,1) ,icor(2,2)  ) ;
"""

   netlist_NoisePSOopt="""
; Top Design: "SpiceTools_lib:NoiseOS_vs_genericMDIF:schematic_PSO_wOPT" 
; Netlist for PSO simultations with noise

Options ResourceUsage=yes UseNutmegFormat=no EnableOptim=yes \\
        TopDesignName="SpiceTools_lib:NoiseOS_vs_genericMDIF:schematic_PSO" DcopOutputNodeVoltages=yes DcopOutputPinCurrents=yes\\
        DcopOutputAllSweepPoints=no DcopOutputDcopType=0

Port:TermG1  Port1 0 Num=1 Z=50 Ohm Noise=yes Temp=27 
Port:TermG2  Port2 0 Num=2 Z=50 Ohm Noise=no Temp=27 

S_Param:SP1 CalcS=yes CalcY=yes CalcZ=no GroupDelayAperture=1e-4 FreqConversion=no FreqConversionPort=1 StatusLevel=2 CalcNoise=yes SortNoise=0 BandwidthForNoise={DF} Hz DevOpPtLevel=0 NoiseInputPort=1 NoiseOutputPort=2 \\
SweepVar="freq" SweepPlan="SP1_stim" OutputPlan="SP1_Output" 

SweepPlan: SP1_stim {SWEEP}

OutputPlan:SP1_Output \\
      Type="Output" \\
      UseEquationNestLevel=yes \\
      EquationNestLevel=2 \\
      UseSavedEquationNestLevel=yes \\
      SavedEquationNestLevel=2 

C:Cpad1   Port1 0     C={Cpad1}
C:Cpad2   Port2 0     C={Cpad2}
C:Cpad12  Port2 Port1 C={Cpad12}

L:Lshort1  Port1          Nin   L={Lsh1} R={Rsh1} Noise=yes  
L:Lshort2  Nout           Port2 L={Lsh2} R={Rsh2} Noise=yes  
L:Lshort3  GNDbeforeShort 0     L={Lsh3} R={Rsh3} Noise=yes  

C:Copen1   Nin  GNDbeforeShort C={Copen1}
C:Copen2   Nout GNDbeforeShort C={Copen2}
C:Copen12  Nin  Nout           C={Copen12}

C:Cpg  Nin  GNDbeforeShort C={Cpg}
C:Cpd  Nout GNDbeforeShort C={Cpd}

L:Lg  Nin        Nint           L={Lg}  Noise=no  
R:Rg  Nintgate   Ngate          R={Rg}  Noise=no 
L:Ld  Noutdrain  Nout           L={Ld}  Noise=no  
R:Rd  Ndrain     Noutdrain      R={Rd}  Noise=no 
L:Ls  Nintsource GNDbeforeShort L={Ls}  Noise=no  
R:Rs  Nsource    Nintsource     R={Rs}  Noise=no 

R:Ri   Nintgs Nsource R={Ri}  Noise=no 
R:RGD  Nintgd Ndrain  R={Rgd} Noise=no 
C:CGD  Ngate  Nintgd  C={Cgd}
C:CGS  Ngate  Nintgs  C={Cgs}
R:Rds  Ndrain Nsource R={Rds} Noise=no 
#uselib "ckt" , "VCCS"
VCCS:GM  Ngate Nsource Ndrain Nsource G={Gm} T=0.0 nsec R1=1e100 Ohm R2=1e100 Ohm F=0.0 GHz 

kb=1.38e-23
DF={DF}
Tout={Tout} opt{{ {ToutMin} to {ToutMax} }}
Tin=290
H22={H22}
H11={H11}
vg_noise=sqrt(4*kb*Tin*DF*H11)
id_noise=sqrt(4*kb*Tout*DF*H22)

V_Source:noiseVg  Nint Nintgate Type="V_Noise" V_Noise=vg_noise SaveCurrent=1 
NoiseCorr2Port:SRC3  CorrCoeff=0.0 Source1="noiseVg" Source2="noiseId" 
I_Source:noiseId  Noutdrain Nintsource Type="I_Noise" I_Noise=id_noise 

Options:Options1 Temp=16.85 Tnom=25 TopologyCheck=yes ForceS_Params=yes V_RelTol=1e-6  \\
I_RelTol=1e-6 GiveAllWarnings=yes MaxWarnings=10 ForceM_Params=yes InitialGuessAnnotation=0  \\
TopologyCheckMessages=no doDeltaAC=yes ReduceSPortRatio=0.5 WarnSOA=yes MaxWarnSOA=5  \\
Census=no MinNodalR=1e-2 Ohm LargeCapThreshold=1e-3 F EnableAssertChecks=no EnableSpareRemoval=yes 

#uselib "ckt" , "DAC"
DAC:DACmeas  File="{noisefile}" Type="gmdif" InterpMode="spline" InterpDom="ri" ExtrapMode="interpMode" iVar1="freq" iVal1=freq iVar2="Vg" iVal2={Vg} iVar3="Vd" iVal3={Vd} 
NF50=file{{DACmeas, "NF50"}}
aele NF50meas=NF50;NF50sim=nf(2);NF50diff=NF50meas-NF50sim;

OptimGoal:OptimGoal1 Expr="NF50diff" SimInstanceName="SP1" Weight=1.0 \\
       IndepVar[1]="freq" SpecLimitLine[1]="OptimGoal1_limit1" 
SpecLimitLine:"OptimGoal1_limit1" Type="Inside" Min=-0.1 Max=0.1 Weight=1.0 \\
       IndepVar[1]="freq" IndepMin[1]={IndepMin} IndepMax[1]={IndepMax} 

Optim:"Optim1" OptimType="hpVMO" ErrorForm="L2" MaxIters=25 P=2 DesiredError=0.0 StatusLevel=4 FinalAnalysis="None"              \\
      NormalizeGoals=yes SetBestValues=yes SaveSolns=yes SaveGoals=no SaveOptimVars=no UpdateDataset=no SaveNominal=no           \\
      SaveAllIterations=no UseAllOptVars=yes UseAllGoals=yes SaveCurrentEF=no InitialTempControlMode=0 NumShootsPerIter=20       \\
      EnableCockpit=no NormalizeError=yes SaveAllTrials=no UseAdvTermCriteria=yes CostRelativeTol=1.0e-8 LimitOfContSmallImprovement=5 

aele ExportNoiseIEEE=write_var("{dev}_PSO_NF_data.dat", WriteMode, " freq NFmin Rn Sopt nf(1) nf(2) te(1) te(2)", Delimiter, NumberFormat, precision,  freq[0,::], NFmin[0,::], Rn[0,::], Sopt[0,::], nf(1)[0,::], nf(2)[0,::], te(1)[0,::], te(2)[0,::] )
aele ExportSpar=write_var("{dev}_PSO_Sparam.dat", WriteMode, " freq S11 S12 S21 S22", Delimiter, NumberFormat, precision,  freq[0,::], S(1,1)[0,::], S(1,2)[0,::], S(2,1)[0,::], S(2,2)[0,::])
aele ExportNoiseCy=write_var("{dev}_PSO_Cy_data.dat", WriteMode, " freq icorr11 icorr12 icorr21 icorr22", Delimiter, NumberFormat, precision,  freq[0,::], icor(1,1)[0,::], icor(1,2)[0,::], icor(2,1)[0,::], icor(2,2)[0,::])
aele ExportTout=write_var("{dev}_PSO_Tout.dat", WriteMode, " Optimat Tout", Delimiter, NumberFormat, precision, Tout[0]);
aele NumberFormat="s";
aele precision=10;
aele WriteMode="W";
aele Delimiter=" ";
"""


   netlist_NoiseOS="""
Options ResourceUsage=yes UseNutmegFormat=no EnableOptim=no \\
        TopDesignName="SpiceTools_lib:NoiseCpdCpg:schematic" DcopOutputNodeVoltages=yes DcopOutputPinCurrents=yes\\
        DcopOutputAllSweepPoints=no DcopOutputDcopType=0

Port:TermG1  Port1 0 Num=1 Z=50 Ohm Noise=yes Temp=27 
Port:TermG2  Port2 0 Num=2 Z=50 Ohm Noise=yes Temp=27 

S_Param:SP1 CalcS=yes CalcY=yes CalcZ=no GroupDelayAperture=1e-4 FreqConversion=no FreqConversionPort=1 StatusLevel=2 CalcNoise=yes SortNoise=0 BandwidthForNoise={DF} Hz DevOpPtLevel=0 \\
SweepVar="freq" SweepPlan="SP1_stim" OutputPlan="SP1_Output" 

SweepPlan: SP1_stim {SWEEP}

OutputPlan:SP1_Output \\
      Type="Output" \\
      UseEquationNestLevel=yes \\
      EquationNestLevel=2 \\
      UseSavedEquationNestLevel=yes \\
      SavedEquationNestLevel=2


C:Cport1   Port1 0     C={Cport1} 
C:Cport2   Port2 0     C={Cport2} 
C:Cport12  Port1 Port2 C={Cport12} 
L:Lshort1  Port1 Nin   L={Lsh1} R={Rsh1} Noise=yes  
L:Lshort2  Nout  Port2 L={Lsh2} R={Rsh2} Noise=yes  

C:Cpg      Nin   0     C={Cpg} 
C:Cpd      Nout  0     C={Cpd}

L:Lg  Nin        Nint       L={Lg} Noise=no  
R:Rd  Ndrain     Noutdrain  R={Rd} Noise=no 
L:Ls  Nintsource 0          L={Ls} Noise=no  
R:Rs  Nsource    Nintsource R={Rs} Noise=no 
R:Rg  Nintgate   Ngate      R={Rg} Noise=no 
L:Ld  Noutdrain  Nout       L={Ld} Noise=no  

R:Ri   Nintgs    Nsource    R={Ri}  Noise=no 
R:RGD  Nintgd    Ndrain     R={Rgd} Noise=no 
C:CGD  Ngate     Nintgd     C={Cgd} 
C:CGS  Ngate     Nintgs     C={Cgs} 
R:Rds  Ndrain    Nsource    R={Rds} Noise=no 
#uselib "ckt" , "VCCS"
VCCS:GM  Ngate Nsource Ndrain Nsource G={Gm} T=0.0 nsec R1=1e100 Ohm R2=1e100 Ohm F=0.0 GHz 

kb=1.38e-23
DF={DF}
Tout={Tout}
Tin=290
H22={H22}
H11={H11}
vg_noise=sqrt(4*kb*Tin*DF*H11)
id_noise=sqrt(4*kb*Tout*DF*H22)

V_Source:noiseVg   Nint      Nintgate   Type="V_Noise" V_Noise=vg_noise SaveCurrent=1 
NoiseCorr2Port:SRC3  CorrCoeff=0.0 Source1="noiseVg" Source2="noiseId" 
I_Source:noiseId   Noutdrain Nintsource Type="I_Noise" I_Noise=id_noise 

 Options:Options1 Temp=16.85 Tnom=16.85 TopologyCheck=yes ForceS_Params=yes GiveAllWarnings=yes  \\
MaxWarnings=10 ForceM_Params=yes InitialGuessAnnotation=0 TopologyCheckMessages=no doDeltaAC=no  \\
ReduceSPortRatio=0.5 WarnSOA=yes MaxWarnSOA=5 Census=no MinNodalR=1e-2 Ohm  \\
LargeCapThreshold=1e-3 F EnableAssertChecks=no EnableSpareRemoval=yes

aele x=write_snp("{dev}_OS_Sparam.s2p",S,"S-par simulation data","Hz","MA",50,50, 1, 9, "  ");
aele temporary1=write_var("{dev}_OS_NFdata.dat", "W", "## freq, NFmin, Rn, Sopt, nf(1), nf(2), te(1), te(2) ", " ", "s", 10, \\
      freq, NFmin, Rn, Sopt, nf(1), nf(2), te(1), te(2) );
aele temporary2=write_var("{dev}_OS_Sparam2.dat", "W", "## S11.r S11.i S12.r S12.i S21.r S21.i S22.r S22.i ", " ", "s", 10, \\
      freq, S(1,1) , S(1,2) ,S(2,1) ,S(2,2)  ) ;
aele temporary3=write_var("{dev}_OS_Cydata.dat", "W", "## icor11.r icor11.i icor12.r icor12.i icor21.r icor21.i icor22.r icor22.i ", " ", "s", 10, \\
      freq, icor(1,1) , icor(1,2) ,icor(2,1) ,icor(2,2)  ) ;
"""

   netlist_NoiseOSopt="""
; Top Design: "SpiceTools_lib:NoiseOS_vs_genericMDIF:schematic_wOPT"
; Netlisted using Hierarchy Policy: "Standard"
Options ResourceUsage=yes UseNutmegFormat=no EnableOptim=yes \\
        TopDesignName="SpiceTools_lib:NoiseOS_vs_genericMDIF:schematic_wOPT" DcopOutputNodeVoltages=yes DcopOutputPinCurrents=yes \\
        DcopOutputAllSweepPoints=no DcopOutputDcopType=0

Port:TermG1  Port1 0 Num=1 Z=50 Ohm Noise=yes Temp=27 
Port:TermG2  Port2 0 Num=2 Z=50 Ohm Noise=yes Temp=27 

S_Param:SP1 CalcS=yes CalcY=yes CalcZ=no GroupDelayAperture=1e-4 FreqConversion=no FreqConversionPort=1 StatusLevel=2 CalcNoise=yes SortNoise=0 \\
          BandwidthForNoise={DF} Hz DevOpPtLevel=0 NoiseInputPort=1 NoiseOutputPort=2 \\
          SweepVar="freq" SweepPlan="SP1_stim" OutputPlan="SP1_Output" 

SweepPlan: SP1_stim {SWEEP}

OutputPlan:SP1_Output \\
      Type="Output" \\
      UseEquationNestLevel=yes \\
      EquationNestLevel=2 \\
      UseSavedEquationNestLevel=yes \\
      SavedEquationNestLevel=2

C:Cport1   Port1 0     C={Cport1} 
C:Cport2   Port2 0     C={Cport2} 
C:Cport12  Port1 Port2 C={Cport12} 
L:Lshort1  Port1 Nin   L={Lsh1} R={Rsh1} Noise=yes  
L:Lshort2  Nout  Port2 L={Lsh2} R={Rsh2} Noise=yes  

C:Cpg      Nin   0     C={Cpg} 
C:Cpd      Nout  0     C={Cpd}

L:Lg  Nin        Nint       L={Lg} Noise=no  
R:Rd  Ndrain     Noutdrain  R={Rd} Noise=no 
L:Ls  Nintsource 0          L={Ls} Noise=no  
R:Rs  Nsource    Nintsource R={Rs} Noise=no 
R:Rg  Nintgate   Ngate      R={Rg} Noise=no 
L:Ld  Noutdrain  Nout       L={Ld} Noise=no  

R:Ri   Nintgs    Nsource    R={Ri}  Noise=no 
R:RGD  Nintgd    Ndrain     R={Rgd} Noise=no 
C:CGD  Ngate     Nintgd     C={Cgd} 
C:CGS  Ngate     Nintgs     C={Cgs} 
R:Rds  Ndrain    Nsource    R={Rds} Noise=no 
#uselib "ckt" , "VCCS"
VCCS:GM  Ngate Nsource Ndrain Nsource G={Gm} T=0.0 nsec R1=1e100 Ohm R2=1e100 Ohm F=0.0 GHz 

kb=1.38e-23
DF={DF}
Tout={Tout} opt{{ {ToutMin} to {ToutMax}  }}
Tin=290
H22={H22}
H11={H11}
vg_noise=sqrt(4*kb*Tin*DF*H11)
id_noise=sqrt(4*kb*Tout*DF*H22)

V_Source:noiseVg   Nint      Nintgate   Type="V_Noise" V_Noise=vg_noise SaveCurrent=1 
NoiseCorr2Port:SRC3  CorrCoeff=0.0 Source1="noiseVg" Source2="noiseId" 
I_Source:noiseId   Noutdrain Nintsource Type="I_Noise" I_Noise=id_noise 

 Options:Options1 Temp=16.85 Tnom=16.85 TopologyCheck=yes ForceS_Params=yes GiveAllWarnings=yes  \\
MaxWarnings=10 ForceM_Params=yes InitialGuessAnnotation=0 TopologyCheckMessages=no doDeltaAC=no  \\
ReduceSPortRatio=0.5 WarnSOA=yes MaxWarnSOA=5 Census=no MinNodalR=1e-2 Ohm  \\
LargeCapThreshold=1e-3 F EnableAssertChecks=no EnableSpareRemoval=yes

#uselib "ckt" , "DAC"
DAC:DACmeas  File="{noisefile}" Type="gmdif" InterpMode="spline" InterpDom="ri" ExtrapMode="interpMode" iVar1="freq" iVal1=freq iVar2="Vg" iVal2={Vg} iVar3="Vd" iVal3={Vd} 
NF50=file{{DACmeas, "NF50"}}
aele NF50meas=NF50;NF50sim=nf(2);NF50diff=NF50meas-NF50sim;

OptimGoal:OptimGoal1 Expr="NF50diff" SimInstanceName="SP1" Weight=1.0 \\
       IndepVar[1]="freq" SpecLimitLine[1]="OptimGoal1_limit1" 
SpecLimitLine:"OptimGoal1_limit1" Type="Inside" Min=-0.1 Max=0.1 Weight=1.0 \\
       IndepVar[1]="freq" IndepMin[1]={IndepMin} IndepMax[1]={IndepMax} 

Optim:"Optim1" OptimType="hpVMO" ErrorForm="L2" MaxIters=25 P=2 DesiredError=0.0 StatusLevel=4 FinalAnalysis="None"              \\
      NormalizeGoals=yes SetBestValues=yes SaveSolns=yes SaveGoals=no SaveOptimVars=no UpdateDataset=no SaveNominal=no           \\
      SaveAllIterations=no UseAllOptVars=yes UseAllGoals=yes SaveCurrentEF=no InitialTempControlMode=0 NumShootsPerIter=20       \\
      EnableCockpit=no NormalizeError=yes SaveAllTrials=no UseAdvTermCriteria=yes CostRelativeTol=1.0e-8 LimitOfContSmallImprovement=5 

aele ExportNoiseIEEE=write_var("{dev}_OS_NF_data.dat", WriteMode, " freq NFmin Rn Sopt nf(1) nf(2) te(1) te(2)", Delimiter, NumberFormat, precision,  freq[0,::], NFmin[0,::], Rn[0,::], Sopt[0,::], nf(1)[0,::], nf(2)[0,::], te(1)[0,::], te(2)[0,::] )
aele ExportSpar=write_var("{dev}_OS_Sparam.dat", WriteMode, " freq S11 S12 S21 S22", Delimiter, NumberFormat, precision,  freq[0,::], S(1,1)[0,::], S(1,2)[0,::], S(2,1)[0,::], S(2,2)[0,::])
aele ExportNoiseCy=write_var("{dev}_OS_Cy_data.dat", WriteMode, " freq icorr11 icorr12 icorr21 icorr22", Delimiter, NumberFormat, precision,  freq[0,::], icor(1,1)[0,::], icor(1,2)[0,::], icor(2,1)[0,::], icor(2,2)[0,::])
aele ExportTout=write_var("{dev}_OS_Tout.dat", WriteMode, " Optimat Tout", Delimiter, NumberFormat, precision, Tout[0]);
aele NumberFormat="s";
aele precision=10;
aele WriteMode="W";
aele Delimiter=" ";  
"""

   netlist_NoiseCpgCpd="""
Options ResourceUsage=yes UseNutmegFormat=no EnableOptim=no \\
        TopDesignName="SpiceTools_lib:NoiseCpdCpg:schematic" DcopOutputNodeVoltages=yes DcopOutputPinCurrents=yes\\
        DcopOutputAllSweepPoints=no DcopOutputDcopType=0
Port:TermG1  Nin 0 Num=1 Z=50 Ohm Noise=yes Temp=27 
Port:TermG2  Nout 0 Num=2 Z=50 Ohm Noise=yes Temp=27 

S_Param:SP1 CalcS=yes CalcY=yes CalcZ=no GroupDelayAperture=1e-4 FreqConversion=no FreqConversionPort=1 StatusLevel=2 CalcNoise=yes SortNoise=0 BandwidthForNoise={DF} Hz DevOpPtLevel=0 \\
SweepVar="freq" SweepPlan="SP1_stim" OutputPlan="SP1_Output" 

SweepPlan: SP1_stim {SWEEP}

OutputPlan:SP1_Output \\
      Type="Output" \\
      UseEquationNestLevel=yes \\
      EquationNestLevel=2 \\
      UseSavedEquationNestLevel=yes \\
      SavedEquationNestLevel=2

L:Lg  Nin Nint L={Lg} Noise=no  
R:Rg  Nintgate Ngate R={Rg} Noise=no 
R:Rd  Ndrain Noutdrain R={Rd} Noise=no 
L:Ld  Noutdrain Nout L={Ld} Noise=no  
L:Ls  Nintsource 0 L={Ls}  Noise=no  
R:Rs  Nsource Nintsource R={Rs} Noise=no
  
R:Ri  Nintgs Nsource R={Ri} Noise=no 

R:RGD  Nintgd Ndrain R={Rgd} Noise=no 
C:CGD  Ngate Nintgd C={Cgd} 
C:CGS  Ngate Nintgs C={Cgs} 
R:Rds  Ndrain Nsource R={Rds} Noise=no 

#uselib "ckt" , "VCCS"
VCCS:GM  Ngate Nsource Ndrain Nsource G={Gm} T=0.0 nsec R1=1e100 Ohm R2=1e100 Ohm F=0.0 GHz 

C:Cpg  Nint 0 C={Cpg} F
C:Cpd  Nout 0 C={Cpd} F

V_Source:noiseVg  Nin Nintgate Type="V_Noise" V_Noise=vg_noise SaveCurrent=1 
NoiseCorr2Port:SRC3  CorrCoeff=0.0 Source1="noiseVg" Source2="noiseId" 
I_Source:noiseId  Noutdrain Nintsource Type="I_Noise" I_Noise=id_noise

kb=1.38e-23
DF={DF}
Tout={Tout}
Tin=290
H22={H22}
H11={H11}
vg_noise=sqrt(4*kb*Tin*DF*H11)
id_noise=sqrt(4*kb*Tout*DF*H22)
 
 Options:Options1 Temp=16.85 Tnom=16.85 TopologyCheck=yes ForceS_Params=yes GiveAllWarnings=yes  \\
MaxWarnings=10 ForceM_Params=yes InitialGuessAnnotation=0 TopologyCheckMessages=no doDeltaAC=no  \\
ReduceSPortRatio=0.5 WarnSOA=yes MaxWarnSOA=5 Census=no MinNodalR=1e-2 Ohm  \\
LargeCapThreshold=1e-3 F EnableAssertChecks=no EnableSpareRemoval=yes

aele x=write_snp("{dev}_CpgCpd_Sparam.s2p",S,"S-par simulation data","Hz","MA",50,50, 1, 9, "  ");
aele temporary1=write_var("{dev}_CpgCpd_NFdata.dat", "W", "## freq, NFmin, Rn, Sopt, nf(1), nf(2), te(1), te(2) ", " ", "s", 10, \\
      freq, NFmin, Rn, Sopt, nf(1), nf(2), te(1), te(2) );
aele temporary2=write_var("{dev}_CpgCpd_Sparam2.dat", "W", "## S11.r S11.i S12.r S12.i S21.r S21.i S22.r S22.i ", " ", "s", 10, \\
      freq, S(1,1) , S(1,2) ,S(2,1) ,S(2,2)  ) ;
aele temporary3=write_var("{dev}_CpgCpd_Cydata.dat", "W", "## icor11.r icor11.i icor12.r icor12.i icor21.r icor21.i icor22.r icor22.i ", " ", "s", 10, \\
      freq, icor(1,1) , icor(1,2) ,icor(2,1) ,icor(2,2)  ) ;
"""
   # nestlist_sys implements the SSEC with the all BEOL, PADS included
   netlist_sys="""
; Top Design: "SSC_lib:MW2216_s4d4w2_PSO:sgslvtn_s4d4w2"
; Netlisted using Hierarchy Policy: "Standard"

Options ResourceUsage=yes UseNutmegFormat=no EnableOptim=no TopDesignName="SSC_lib:MW2216_s4d4w2_PSO:sgslvtn_s4d4w2_Vg=0.8V" DcopOutputNodeVoltages=yes DcopOutputPinCurrents=yes DcopOutputAllSweepPoints=no DcopOutputDcopType=0
C:Ypad1  N__7 0 C={Ypad1} F
C:C2  N__7 N__5 C={Ypad3} F
#uselib "ckt" , "VCCS"
VCCS:SRC1  N__9 N__4 N__8 N__4 G={Gm} S R1=1e100 Ohm R2=1e100 Ohm F=0.0 GHz
V_Source:SRC3  N__12 N__9 Type="V_Noise" V_Noise=vg_noise SaveCurrent=1
S_Param:SP1 CalcS=yes CalcY=yes CalcZ=no GroupDelayAperture=1e-4 FreqConversion=no FreqConversionPort=1 StatusLevel=2 CalcNoise=yes SortNoise=0 BandwidthForNoise={DF} Hz DevOpPtLevel=0 \
SweepVar="freq" SweepPlan="SP1_stim" OutputPlan="SP1_Output"

SweepPlan: SP1_stim Start=100 MHz Stop=50 GHz Step=100 MHz

OutputPlan:SP1_Output \\
      Type="Output" \\
      UseEquationNestLevel=yes \\
      EquationNestLevel=2 \\
      UseSavedEquationNestLevel=yes \\
      SavedEquationNestLevel=2

#load "python","LinearCollapse"
Component Module="LinearCollapse" Type="ModelExtractor" NetworkRepresentation=2
NoiseCorr2Port:SRC6  CorrCoeff=0 Source1="SRC3" Source2="SRC4"
I_Source:SRC4  N__8 N__3 Type="I_Noise" I_Noise=id_noise
R:Rd  N__8 N__6 R={Rd} Ohm Noise=no
R:Rg  N__11 N__12 R={Rg} Ohm Noise=no
C:Yvias  N__2 N__3 C={Yvia1} F
C:Cgs  N__15 N__4 C={Cgs} F
R:Ri  N__9 N__15 R={Ri} Ohm Noise=no 
C:Yvia2  N__0 N__3 C={Yvia2} F
R:Rs  N__4 N__10 R={Rs} Ohm Noise=no
R:Gds  N__8 N__4 R={Gds} Ohm Noise=no
C:C1  N__2 N__0 C={Yvia3} F
C:Cgd  N__9 N__16 C={Cgd} F
R:Rgd  N__16 N__8 R={Rgd} Ohm Noise=no 
C:Cpd  N__0 N__3 C={Cpd} F
C:Cpg  N__2 N__3 C={Cpg} F
C:Ypad2  N__5 0 C={Ypad2} F
L:L3  N__6 N__0 L={Ld} H Noise=yes
L:L1  N__2 N__11 L={Lg} H Noise=yes
L:L2  N__3 N__10 L={Ls} H Noise=yes
L:L4  N__3 0 L={Ls3} H R={Rs3} Ohm Temp=27 Noise=yes
L:Zs1  N__7 N__2 L={Ls1} H R={Rs1} Ohm Temp=27 Noise=no
L:Zs2  N__0 N__5 L={Ls2} H R={Rs2} Ohm Temp=27 Noise=no
;;aele H=stoh(S);

kb=1.38e-23
DF={DF}
Tout={Tout}
id_noise=sqrt(4*kb*Tout*DF*H22ext)
vg_noise=sqrt(4*kb*Tin*DF*H11ext)
H11ext={H11}
H22ext={H22}
Tin=290
#uselib "ckt" , "DAC"
;;DAC:DAC1  File="/home/310.3.2-Theses/OK252244/ADS/SSC_wrk/data/Hparam.ds" Type="dscr" InterpMode="linear" InterpDom="ri" ExtrapMode="interpMode"

Options:Options1 Temp=16.85 Tnom=16.85 TopologyCheck=yes ForceS_Params=yes GiveAllWarnings=yes  \\
MaxWarnings=10 ForceM_Params=yes InitialGuessAnnotation=0 TopologyCheckMessages=no doDeltaAC=no  \\
ReduceSPortRatio=0.5 WarnSOA=yes MaxWarnSOA=5 Census=no MinNodalR=1e-2 Ohm  \\
LargeCapThreshold=1e-3 F EnableAssertChecks=no EnableSpareRemoval=yes

Port:Term1  N__7 0 Num=1 Z=50 Ohm Noise=yes Temp=27
Port:Term2  N__5 0 Num=2 Z=50 Ohm Noise=yes Temp=27

aele x=write_snp("{dev}_Sparam_sys.s2p",S,"S-par simulation data","Hz","MA",50,50, 1, 9, "  ");

aele temporary1=write_var("{dev}_NFdata_sys.dat", "W", "## freq, NFmin, Rn, Sopt, nf(1), nf(2), te(1), te(2) ", " ", "s", 10, \\
      freq, NFmin, Rn, Sopt, nf(1), nf(2), te(1), te(2) );
aele temporary2=write_var("{dev}_Sparam2_sys.dat", "W", "## S11.r S11.i S12.r S12.i S21.r S21.i S22.r S22.i ", " ", "s", 10, \\
      freq, S(1,1) , S(1,2) ,S(2,1) ,S(2,2)  ) ;
aele temporary3=write_var("{dev}_Cydata_sys.dat", "W", "## icor11.r icor11.i icor12.r icor12.i icor21.r icor21.i icor22.r icor22.i ", " ", "s", 10, \\
      freq, icor(1,1) , icor(1,2) ,icor(2,1) ,icor(2,2)  ) ;
"""
   # netlist_deembedded implements the de-embedded SSEC, i.e. including Cpg and Cpd
   netlist_deembedded="""; Top Design: "SSC_lib:MW2216_s4d4w2_PSO:sgslvtn_s4d4w2"
; Netlisted using Hierarchy Policy: "Standard"

Options ResourceUsage=yes UseNutmegFormat=no EnableOptim=no TopDesignName="SSC_lib:MW2216_s4d4w2_PSO:sgslvtn_s4d4w2_Vg=0.8V" DcopOutputNodeVoltages=yes DcopOutputPinCurrents=yes DcopOutputAllSweepPoints=no DcopOutputDcopType=0
#uselib "ckt" , "VCCS"
S_Param:SP1 CalcS=yes CalcY=yes CalcZ=no GroupDelayAperture=1e-4 FreqConversion=no FreqConversionPort=1 StatusLevel=2 CalcNoise=yes SortNoise=0 BandwidthForNoise={DF} Hz DevOpPtLevel=0 \
SweepVar="freq" SweepPlan="SP1_stim" OutputPlan="SP1_Output" 

SweepPlan: SP1_stim Start=100 MHz Stop=50 GHz Step=100 MHz 

OutputPlan:SP1_Output \\
      Type="Output" \\
      UseEquationNestLevel=yes \\
      EquationNestLevel=2 \\
      UseSavedEquationNestLevel=yes \\
      SavedEquationNestLevel=2

#load "python","LinearCollapse"
Component Module="LinearCollapse" Type="ModelExtractor" NetworkRepresentation=2
C:Cgs  N__15 N__4 C={Cgs} F
R:R1  N__9 N__15 R={Ri} Ohm Noise=no 
R:Gds  N__8 N__4 R={Gds} Ohm Noise=no 
C:Cgd  N__9 N__16 C={Cgd} F
R:Rgd  N__16 N__8 R={Rgd} Ohm Noise=no 
VCCS:Gm  N__9 N__4 N__8 N__4 G={Gm} S R1=1e100 Ohm R2=1e100 Ohm F=0.0 GHz 

Options:Options1 Temp=16.85 Tnom=16.85 TopologyCheck=yes ForceS_Params=yes GiveAllWarnings=yes  \\
MaxWarnings=10 ForceM_Params=yes InitialGuessAnnotation=0 TopologyCheckMessages=no doDeltaAC=no  \\
ReduceSPortRatio=0.5 WarnSOA=yes MaxWarnSOA=5 Census=no MinNodalR=1e-2 Ohm  \\
LargeCapThreshold=1e-3 F EnableAssertChecks=no EnableSpareRemoval=yes 

Port:Term1  N__11 0 Num=1 Z=50 Ohm Noise=yes Temp=27 
Port:Term2  N__3 0 Num=2 Z=50 Ohm Noise=yes Temp=27 
R:Rd  N__8 N__3 R={Rd} Ohm Noise=no 
I_Source:SRC4  N__8 0 Type="I_Noise" I_Noise=id_noise 
V_Source:SRC3  N__12 N__9 Type="V_Noise" V_Noise=vg_noise SaveCurrent=1 
R:Rg  N__11 N__12 R={Rg} Ohm Noise=no 
R:Rs  N__4 0 R={Rs} Ohm Noise=no 
aele H=stoh(S);

kb=1.38e-23
DF=1
Tout={Tout}
id_noise=sqrt(4*kb*Tout*DF*H22ext)
vg_noise=sqrt(4*kb*Tin*DF*H11ext)
H11ext={H11}
H22ext={H22}
Tin=290
NoiseCorr2Port:SRC6  CorrCoeff=0 Source1="SRC3" Source2="SRC4" 

aele x=write_snp("borrame_demb.s2p",S,"S-par simulation data","Hz","MA",50,50, 1, 9, "  ");
aele temporary1=write_var("NFdata_demb_{dev}.s2p", "W", "## freq, NFmin, Rn, Sopt, nf(1), nf(2), te(1), te(2) ", " ", "s", 10, \\
      freq, NFmin, Rn, Sopt, nf(1), nf(2), te(1), te(2) );
aele temporary2=write_var("SpData_demb_{dev}.s2p", "W", "## S11.r S11.i S12.r S12.i S21.r S21.i S22.r S22.i ", " ", "s", 10, \\
      freq, S(1,1) , S(1,2) ,S(2,1) ,S(2,2)  ) ;
aele temporary3=write_var("CyData_demb_{dev}.s2p", "W", "## icor11.r icor11.i icor12.r icor12.i icor21.r icor21.i icor22.r icor22.i ", " ", "s", 10, \\
      freq, icor(1,1) , icor(1,2) ,icor(2,1) ,icor(2,2)  ) ;
"""
   #netlist_int is the intrinsic device
   netlist_int="""; Top Design: "SSC_lib:MW2216_s4d4w2_PSO:sgslvtn_s4d4w2"
; Netlisted using Hierarchy Policy: "Standard"

Options ResourceUsage=yes UseNutmegFormat=no EnableOptim=no TopDesignName="SSC_lib:MW2216_s4d4w2_PSO:sgslvtn_s4d4w2_Vg=0.8V" DcopOutputNodeVoltages=yes DcopOutputPinCurrents=yes DcopOutputAllSweepPoints=no DcopOutputDcopType=0
#uselib "ckt" , "VCCS"
S_Param:SP1 CalcS=yes CalcY=yes CalcZ=no GroupDelayAperture=1e-4 FreqConversion=no FreqConversionPort=1 StatusLevel=2 CalcNoise=yes SortNoise=0 BandwidthForNoise=100 Hz DevOpPtLevel=0 \
SweepVar="freq" SweepPlan="SP1_stim" OutputPlan="SP1_Output" 

SweepPlan: SP1_stim {SWEEP}

OutputPlan:SP1_Output \
      Type="Output" \
      UseEquationNestLevel=yes \
      EquationNestLevel=2 \
      UseSavedEquationNestLevel=yes \
      SavedEquationNestLevel=2

#load "python","LinearCollapse"
Component Module="LinearCollapse" Type="ModelExtractor" NetworkRepresentation=2
C:Cgs  N__15 0 C={Cgs} F
R:Ri  N__9 N__15 R={Ri} Ohm Noise=no 
R:Gds  N__3 0 R={Rds} Ohm Noise=no 
C:Cgd  N__9 N__16 C={Cgd} F
R:Rgd  N__16 N__3 R={Rgd} Ohm Noise=no 
VCCS:Gm  N__9 0 N__3 0 G={Gm} S R1=1e100 Ohm R2=1e100 Ohm F=0.0 GHz 

Options:Options1 Temp=16.85 Tnom=16.85 TopologyCheck=yes ForceS_Params=yes GiveAllWarnings=yes  \
MaxWarnings=10 ForceM_Params=yes InitialGuessAnnotation=0 TopologyCheckMessages=no doDeltaAC=no  \
ReduceSPortRatio=0.5 WarnSOA=yes MaxWarnSOA=5 Census=no MinNodalR=1e-2 Ohm  \
LargeCapThreshold=1e-3 F EnableAssertChecks=no EnableSpareRemoval=yes 

Port:Term1  N__9 0 Num=1 Z=50 Ohm Noise=yes Temp=27 
Port:Term2  N__3 0 Num=2 Z=50 Ohm Noise=yes Temp=27 

aele x=write_snp("{dev}_Sparam_int.s2p",S,"S-par simulation data","Hz","MA",50,50, 1, 9, "  ");

aele temporary2=write_var("{dev}_Sparam2_int.s2p", "W", "## S11.r S11.i S12.r S12.i S21.r S21.i S22.r S22.i ", " ", "s", 10, \
      freq, S(1,1) , S(1,2) ,S(2,1) ,S(2,2)  ) ;
"""

   def __init__(self):
      self.fHz=np.linspace(0.5,50,100)*1.0e9
      if windows:
              my_env = os.environ.copy()
              my_env['HPEESOF_DIR']="C:\\Program Files\\Keysight\\ADS2022_update1"
              my_env['COMPL_DIR']="C:\\Program Files\\Keysight\\ADS2022_update1"
              my_env['SIMARCH']="win32_64"
              my_env['PATH']=my_env['HPEESOF_DIR']+"\\bin;"+  \
                 my_env['HPEESOF_DIR']+"\\bin\\"+my_env['SIMARCH']+";"+   \
                 my_env['HPEESOF_DIR']+"\\lib\\"+my_env['SIMARCH']+";"+   \
                 my_env['HPEESOF_DIR']+"\\adsptolemy\\lib."+my_env['SIMARCH']+";"+  \
                 my_env['PATH']
              self.my_env=my_env
      else: 
          self.my_env = os.environ.copy()
          self.my_env["SIMARCH"]="linux_x86_64"
          self.my_env["HPEESOF_DIR"]="/home/cao/agilent/ADS2020_00"
          self.my_env["PATH"]="/home/cao/agilent/ADS2020_00/bin/:/home/cao/agilent/ADS2020_00/linux_x86_64/bin:"+self.my_env["PATH"]
          self.my_env["LD_LIBRARY_PATH"]="/home/cao/agilent/ADS2020_00/lib/linux_x86_64:/home/cao/agilent/ADS2020_00/lib:/home/cao/agilent/ADS2020_00/lib/linux_x86_64/xcbSupport://home/cao/agilent/ADS2020_00/adsptolemy/lib.linux_x86_64:.:"+self.my_env['LD_LIBRARY_PATH']

   def parseVarFileNoiseData(self, file):
     f=file.readlines()
     Sstop=np.size(f)
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

   def Hparam_symb(self, Rg, Rd, Rs, Gm, Gds, Cgd, Cgs, Ri=0, Rgd=0, seed=None, debug=False):
      # we need to add explicitly the frequeny sweep that is hardocode in the 
      # ADS netlist in the standard function
      # function should work if all the circuit elements are passed as reals or as arrays.
      # in the Mathematica file I can handle tau in the Gm, not included here
      omega=2*np.pi*self.fHz

      y11num=1j*Cgd*omega+1j*Cgs*omega-Cgd*Cgs*Rgd*omega*omega-Cgd*Cgs*Ri*omega*omega
      y11den=(-1j+Cgd*Rgd*omega)*(-1j+Cgs*Ri*omega)
      y11=-y11num/y11den
      #print(f"ga")    
      y12=(Cgd*omega)/(1j-Cgd*Rgd*omega)

      y21num=Gm*(-1-1j*Cgd*Rgd*omega)+Cgd*omega*(1j-Cgs*Ri*omega)
      y21den=(-1j+Cgd*Rgd*omega)*(-1j+Cgs*Ri*omega)

      y21=y21num/y21den

      y22=Gds+(omega*Cgd)/(-1j+Cgd*Rgd*omega)

      Y=np.empty([len(omega), 2, 2], dtype=np.complex128)
      Z=np.empty([len(omega), 2, 2], dtype=np.complex128)
      for i,om in enumerate(omega):
         Y[i,0,0]=y11[i]
         Y[i,1,1]=y22[i]
         Y[i,1,0]=y21[i]
         Y[i,0,1]=y12[i]
         Z[i,0,0]=Rs+Rg
         Z[i,1,1]=Rs+Rd
         Z[i,1,0]=Rs
         Z[i,0,1]=Rs
      tmp1=rft.YtoZ_2P(Y)
      tmp2=Z+tmp1
      H=rft.ZtoH_2P(tmp2)
      Y2=rft.ZtoY_2P(tmp2)
      #H=rft.ZtoH_2P(Z+rft.YtoZ_2P(Y))
      reH11=np.mean(np.real(H[:,0,0])) 
      reH22=np.mean(np.real(H[:,1,1])) 
 #     if debug:
 #         return reH11,reH22,H,Y2
 #     else:
      return reH11, reH22
      
      
   def Hparam(self, Rg, Rd, Rs, Gm, Gds, Cgd, Cgs, Ri=0, Rgd=0,
              SWEEP="Start=100 MHz Stop=50 GHz Step=100 MHz",
              seed=None, debug=False):
      """
      what is called Gds is in reality 1/Gds, Gds must be passed as Ohms!
      In netlist: 
         Rds 1 2 R=Gds
      and it does not make much sense!
      Non ha senso lo so!
      """
      if seed is None:
         val=repr(np.random.randint(9999999))
      else:
         val=repr(seed)

      netl=self.netlist_Hparam.format(Rg=Rg,Rd=Rd,Rs=Rs, Gm=Gm, Gds=Gds, Cgd=Cgd, Cgs=Cgs, Ri=Ri, Rgd=Rgd, SWEEP=SWEEP, dev=val)
      fname="Hp_"+val+".cir"
      f = open(self.PathADS+fname, 'w')
      f.write(netl)
      f.close()
      #time.sleep(1) # may be needed on NFS or samba filesystems
      if debug : print(f"Ready to process {fname}")
      # this is old, low level syntax
      # os.subprocess.Popen(["adssim",fname],env=my_env,cwd=self.PathADS)
      # os.subprocess.run was added in python 3.5 
      #resu=subprocess.run(["adssim",fname],shell=True,env=self.my_env,capture_output=True,cwd=self.PathADS)
      #resu=subprocess.run(["adssim",fname],env=self.my_env,capture_output=True,cwd=self.PathADS)
      #resu=subprocess.run(["adssim",fname],env=self.my_env,capture_output=True)
      #resu=subprocess.run(["adssim",fname],env=self.my_env,capture_output=True)
      ## This Works
      # Popen looks be low level and still used by the run. So in principle it  is 
      # more portable than run
      # adssim is a script around hpeesofsim
      # there is also hpeesofsess but I cannot find documentation on its use
      if not windows:
          resu=subprocess.Popen("hpeesofsim "+fname,shell=True,env=self.my_env,cwd=self.PathADS,stderr=subprocess.PIPE,stdout=subprocess.PIPE)
          resu.wait()
          stdout, stderr = resu.communicate()
          if debug: print(stdout)
          if len(stderr)>0: print(stderr)
      else:
          command=["hpeesofsim.exe", fname]
          resu=subprocess.run(command, shell=True, env=self.my_env, cwd=self.PathADS, stderr=subprocess.PIPE,
                         stdout=subprocess.PIPE)
          if debug: print(str(resu.stdout, 'utf-8').replace('"','').strip() )
          # leave this print active. If it is not there, there is not enough lag
          # between the subprocess.run and the next line i.e. output=open(self.PathADS+val ...
          # and the system won't find the file because it is in the 
          # process of being saved!!!

#      output = open(self.PathADS+val+"_Sparam.s2p", "r")
      output = open(self.PathADS+val+"_Sparam2.s2p", "r")
      fr, Sdut=rfdataanalysis.getSparamfromS2P(output)
      Hdut=rft.StoH(Sdut)
      Ydut=rft.StoY(Sdut)
      H11=np.mean(np.real(Hdut[:,0,0]))
      H22=np.mean(np.real(Hdut[:,1,1])) 
#      if debug:
#          return H11,H22,Hdut,Ydut, Sdut 
#      else:
      return H11,H22

   def CpgCpd(self, Rg, Rd, Rs, Ls, Ld, Lg, Cpg, Cpd, Gm, Rds, Cgd, Cgs, Ri=0, Rgd=0,
              SWEEP="Start=100 MHz Stop=50 GHz Step=100 MHz",
              seed=None, debug=False):
      """
      what is called Gds is in reality 1/Gds, Gds is passed as Ohms! ->
      this bug has now been corrected, what was called Gds now is called Rds and
      this Rds has to be passed. in all OK flow, h stores Rds in a variable called 'Gds'
      return Sparams!
      """
      if seed is None:
         val=repr(np.random.randint(9999999))
      else:
         val=repr(seed)

      netl=self.netlist_CpgCpd.format(Rg=Rg,Rd=Rd,Rs=Rs,Lg=Lg,Ld=Ld,Ls=Ls, Cpg=Cpg, Cpd=Cpd, Gm=Gm, Rds=Rds, Cgd=Cgd, Cgs=Cgs, Ri=Ri, Rgd=Rgd, SWEEP=SWEEP, dev=val)
      fname="CpgCpd_"+val+".cir"
      f = open(self.PathADS+fname, 'w')
      f.write(netl)
      f.close()
      #time.sleep(1) # may be needed on NFS or samba filesystems
      if debug : print(f"Ready to process {fname}")
      if not windows:
          resu=subprocess.Popen("hpeesofsim "+fname,shell=True,env=self.my_env,cwd=self.PathADS,stderr=subprocess.PIPE,stdout=subprocess.PIPE)
          resu.wait()
          stdout, stderr = resu.communicate()
          if debug: print(stdout)
          if len(stderr)>0: print(stderr)
      else:
          command=["hpeesofsim.exe", fname]
          resu=subprocess.run(command, shell=True, env=self.my_env, cwd=self.PathADS, stderr=subprocess.PIPE,
                         stdout=subprocess.PIPE)
          print(str(resu.stdout, 'utf-8').replace('"','').strip() )

      output = open(self.PathADS+val+"_CpgCpd_Sparam2.s2p", "r")
      fr, Sdut=rfdataanalysis.getSparamfromS2P(output)
      return Sdut 

   def SimInt(self, Gm, Gds, Cgd, Cgs, Ri=0, Rgd=0,
           SWEEP="Start=100 MHz Stop=50 GHz Step=100 MHz", seed=None, debug=False):
      """
      Funzione che che dovrebbe essere sostituita da un calcolo analitico
      """
      if seed is None:
         val=repr(np.random.randint(9999999))
      else:
         val=repr(seed)
      netl=self.netlist_int.format( Gm=Gm, Rds=(1.0/Gds), Cgd=Cgd, Cgs=Cgs, Ri=Ri, Rgd=Rgd, SWEEP=SWEEP, dev=val)
      fname=val+"Spint.cir"
      f = open(self.PathADS+fname, 'w')
      f.write(netl)
      f.close()
      #time.sleep(1) # may be needed on NFS or samba filesystem
      if debug : print(f"Processing {self.PathADS+fname}")
      if not windows:
          resu=subprocess.Popen("hpeesofsim "+fname,shell=True,env=self.my_env,cwd=self.PathADS,stderr=subprocess.PIPE,stdout=subprocess.PIPE)
          resu.wait()
          stdout, stderr = resu.communicate()
      else:
          command=["hpeesofsim.exe", fname]
          resu=subprocess.run(command, shell=True, env=self.my_env, cwd=self.PathADS, stderr=subprocess.PIPE,
                         stdout=subprocess.PIPE)
          stdout = str(resu.stdout, 'utf-8').replace('"','').strip()  
          stderr = str(resu.stderr, 'utf-8').replace('"','').strip()  

      if debug: print(stdout)
      if stderr: print(stderr)

      output = open(self.PathADS+val+"_Sparam2_int.s2p", "r")
      fr, Sdut=rfdataanalysis.getSparamfromS2P(output)
      output.close()
      return Sdut

   def NoiseCpgCpd(self, Rg, Rd, Rs, Ls, Ld, Lg, Cpg, Cpd, Gm, Rds, Cgd, Cgs, H11, H22, DF, Tout, Ri=0, Rgd=0,
              SWEEP="Start=100 MHz Stop=50 GHz Step=100 MHz",
              seed=None, debug=False):
      """
      what is called Gds is in reality 1/Gds, Gds is passed as Ohms! ->
      this bug has now been corrected, what was called Gds now is called Rds and
      this Rds has to be passed. in all OK flow, h stores Rds in a variable called 'Gds'
      return Sparams!
      """
      IEEEnoise=dict()
      if seed is None:
         val=repr(np.random.randint(9999999))
      else:
         val=repr(seed)

      netl=self.netlist_NoiseCpgCpd.format(Rg=Rg,Rd=Rd,Rs=Rs,Lg=Lg,Ld=Ld,Ls=Ls, Cpg=Cpg, Cpd=Cpd, Gm=Gm, Rds=Rds, Cgd=Cgd, Cgs=Cgs, Ri=Ri, Rgd=Rgd, SWEEP=SWEEP, dev=val, DF=DF, Tout=Tout, H11=H11, H22=H22)
      fname=val+"CpgCpd.cir"
      f = open(self.PathADS+fname, 'w')
      f.write(netl)
      f.close()
      #time.sleep(1) # may be needed on NFS or samba filesystems
      if debug : print(f"Ready to process {fname}")
      if not windows:
          resu=subprocess.Popen("hpeesofsim "+fname,shell=True,env=self.my_env,cwd=self.PathADS,stderr=subprocess.PIPE,stdout=subprocess.PIPE)
          resu.wait()
          stdout, stderr = resu.communicate()
          if debug: print(stdout)
          if len(stderr)>0: print(stderr)
      else:
          command=["hpeesofsim.exe", fname]
          resu=subprocess.run(command, shell=True, env=self.my_env, cwd=self.PathADS, stderr=subprocess.PIPE,
                         stdout=subprocess.PIPE)
          stdout = str(resu.stdout, 'utf-8').replace('"','').strip()  
          stderr = str(resu.stderr, 'utf-8').replace('"','').strip() 

      if debug: print(stdout)
      if stderr: print(stderr)
          
      soutput = open(self.PathADS+val+"_CpgCpd_Sparam2.dat", "r")
      noutput = open(self.PathADS+val+"_CpgCpd_NFdata.dat", "r")
      coutput = open(self.PathADS+val+"_CpgCpd_Cydata.dat", "r")
      f, Sp = rfdataanalysis.getSparamfromS2P(soutput)
      NFmin, Rn, Sopt, NF, TE = self.parseVarFileNoiseData(noutput)
      fc, Cy = rfdataanalysis.getSparamfromS2P(coutput)
      soutput.close()
      noutput.close()
      coutput.close()
      # return f, Sp, NFmin, Rn, Sopt, NF, TE, Cy
      # create a dictionary for noise results
      IEEEnoise['F']=f
      IEEEnoise['simSp']=Sp
      IEEEnoise['NFmin']=NFmin
      IEEEnoise['Rn']=Rn
      IEEEnoise['Sopt']=Sopt
      IEEEnoise['simNF50']=NF
      IEEEnoise['simTe']=TE
      IEEEnoise['Cy']=Cy
      IEEEnoise['Tout']=Tout
      return IEEEnoise

   def NoisePSOopt(self, Cpad1, Cpad2, Cpad12, 
               Rsh1, Rsh2, Rsh3, Lsh1, Lsh2,  Lsh3,
               Copen1, Copen2, Copen12, 
               Rg, Rd, Rs, Ls, Ld, Lg, Cpg, Cpd,
               Gm, Rds, Cgd, Cgs, H11, H22, DF, 
               noisefile, Vg, Vd,
               Tout, ToutMin=100, ToutMax=8000, IndepMin=10e9, IndepMax=30e9, 
               Ri=0, Rgd=0,
               SWEEP="Start=100 MHz Stop=50 GHz Step=100 MHz",
               seed=None, debug=False):
      """
      what is called Gds is in reality 1/Gds, Gds is passed as Ohms! ->
      this bug has now been corrected, what was called Gds now is called Rds and
      this Rds has to be passed. in all OK flow, h stores Rds in a variable called 'Gds'
      return Sparams!
      """
      if seed is None:
         val=repr(np.random.randint(9999999))
      else:
         val=repr(seed)

      netl=self.netlist_NoisePSOopt.format(Cpad1=Cpad1, Cpad2=Cpad2, Cpad12=Cpad12,
              Rsh1=Rsh1, Rsh2=Rsh2, Rsh3=Rsh3, Lsh1=Lsh1, Lsh2=Lsh2, Lsh3=Lsh3,
              Copen1=Copen1, Copen2=Copen2, Copen12=Copen12,
              Rg=Rg, Rd=Rd, Rs=Rs, Lg=Lg, Ld=Ld, Ls=Ls, Cpg=Cpg, Cpd=Cpd,
              Gm=Gm, Rds=Rds, Cgd=Cgd, Cgs=Cgs, Ri=Ri, Rgd=Rgd, noisefile=noisefile, Vg=Vg, Vd=Vd, IndepMin=IndepMin, IndepMax=IndepMax,
              SWEEP=SWEEP, dev=val, DF=DF, Tout=Tout, ToutMin=ToutMin, ToutMax=ToutMax, H11=H11, H22=H22)
      IEEEnoise=dict()
      fname=val+"_PSOopt.cir"
      f = open(self.PathADS+fname, 'w')
      f.write(netl)
      f.close()
      #time.sleep(1) # may be needed on NFS or samba filesystems
      if debug : print(f"Ready to process {fname}")
      if not windows:
          resu=subprocess.Popen("hpeesofsim "+fname,shell=True,env=self.my_env,cwd=self.PathADS,stderr=subprocess.PIPE,stdout=subprocess.PIPE)
          resu.wait()
          stdout, stderr = resu.communicate()
      else:
          command=["hpeesofsim.exe", fname]
          resu=subprocess.run(command, shell=True, env=self.my_env, cwd=self.PathADS, stderr=subprocess.PIPE,
                         stdout=subprocess.PIPE)
          stdout = str(resu.stdout, 'utf-8').replace('"','').strip()  
          stderr = str(resu.stderr, 'utf-8').replace('"','').strip()  


      if debug: print(stdout)
      if stderr: print(stderr)
      
      soutput = open(self.PathADS+val+"_PSO_Sparam.dat"  , "r")
      noutput = open(self.PathADS+val+"_PSO_NF_data.dat" , "r")
      coutput = open(self.PathADS+val+"_PSO_Cy_data.dat" , "r")
      toutput = open(self.PathADS+val+"_PSO_Tout.dat"    , "r")
      f, Sp = rfdataanalysis.getSparamfromS2P(soutput)
      NFmin, Rn, Sopt, NF, TE = self.parseVarFileNoiseData(noutput)
      fc, Cy = rfdataanalysis.getSparamfromS2P(coutput)
      # load the temperature from the second line of the file :-(
      tmp = toutput.readlines()
      Tout=float(tmp[1].strip());
      if debug: print(f"Tout is {Tout}")
      soutput.close()
      noutput.close()
      coutput.close()
      toutput.close()
      # return Tout, f, Sp, NFmin, Rn, Sopt, NF, TE, Cy
      IEEEnoise['F']=f
      IEEEnoise['simSp']=Sp
      IEEEnoise['NFmin']=NFmin
      IEEEnoise['Rn']=Rn
      IEEEnoise['Sopt']=Sopt
      IEEEnoise['simNF50']=NF
      IEEEnoise['simTe']=TE
      IEEEnoise['Cy']=Cy
      IEEEnoise['Tout']=Tout
      return IEEEnoise


   def NoiseOSopt(self, Cport1, Cport2,Cport12, Rsh1, Rsh2, Lsh1, Lsh2, 
               Rg, Rd, Rs, Ls, Ld, Lg, Cpg, Cpd,
               Gm, Rds, Cgd, Cgs, H11, H22, DF, 
               noisefile, Vg, Vd,
               Tout, ToutMin=100, ToutMax=8000, IndepMin=10e9, IndepMax=30e9, 
               Ri=0, Rgd=0,
               SWEEP="Start=100 MHz Stop=50 GHz Step=100 MHz",
               seed=None, debug=False):
      """
      what is called Gds is in reality 1/Gds, Gds is passed as Ohms! ->
      this bug has now been corrected, what was called Gds now is called Rds and
      this Rds has to be passed. in all OK flow, h stores Rds in a variable called 'Gds'
      return Sparams!
      """
      IEEEnoise=dict()

      if seed is None:
         val=repr(np.random.randint(9999999))
      else:
         val=repr(seed)

      netl=self.netlist_NoiseOSopt.format(Cport1=Cport1, Cport2=Cport2, 
              Cport12=Cport12, Rsh1=Rsh1, Rsh2=Rsh2, Lsh1=Lsh1, Lsh2=Lsh2,
              Rg=Rg,Rd=Rd,Rs=Rs,Lg=Lg,Ld=Ld,Ls=Ls, Cpg=Cpg, Cpd=Cpd,
              Gm=Gm, Rds=Rds, Cgd=Cgd, Cgs=Cgs, Ri=Ri, Rgd=Rgd, noisefile=noisefile, Vg=Vg, Vd=Vd, IndepMin=IndepMin, IndepMax=IndepMax,
              SWEEP=SWEEP, dev=val, DF=DF, Tout=Tout, ToutMin=ToutMin, ToutMax=ToutMax, H11=H11, H22=H22)
      fname=val+"_OSopt.cir"
      f = open(self.PathADS+fname, 'w')
      f.write(netl)
      f.close()
      #time.sleep(1) # may be needed on NFS or samba filesystems
      if debug : print(f"Ready to process {fname}")
      if not windows:
          resu=subprocess.Popen("hpeesofsim "+fname,shell=True,env=self.my_env,cwd=self.PathADS,stderr=subprocess.PIPE,stdout=subprocess.PIPE)
          resu.wait()
          stdout, stderr = resu.communicate()
      else:
          command=["hpeesofsim.exe", fname]
          resu=subprocess.run(command, shell=True, env=self.my_env, cwd=self.PathADS, stderr=subprocess.PIPE,
                         stdout=subprocess.PIPE)
          stdout = str(resu.stdout, 'utf-8').replace('"','').strip()  
          stderr = str(resu.stderr, 'utf-8').replace('"','').strip()  


      if debug: print(stdout)
      if stderr: print(stderr)
      
      soutput = open(self.PathADS+val+"_OS_Sparam.dat"  , "r")
      noutput = open(self.PathADS+val+"_OS_NF_data.dat" , "r")
      coutput = open(self.PathADS+val+"_OS_Cy_data.dat" , "r")
      toutput = open(self.PathADS+val+"_OS_Tout.dat"    , "r")
      f, Sp = rfdataanalysis.getSparamfromS2P(soutput)
      NFmin, Rn, Sopt, NF, TE = self.parseVarFileNoiseData(noutput)
      fc, Cy = rfdataanalysis.getSparamfromS2P(coutput)
      # load the temperature from the second line of the file :-(
      tmp = toutput.readlines()
      Tout=float(tmp[1].strip());
      if debug: print(f"Tout is {Tout}")
      soutput.close()
      noutput.close()
      coutput.close()
      toutput.close()
      # create a dictionary for noise results
      IEEEnoise['F']=f
      IEEEnoise['simSp']=Sp
      IEEEnoise['NFmin']=NFmin
      IEEEnoise['Rn']=Rn
      IEEEnoise['Sopt']=Sopt
      IEEEnoise['simNF50']=NF
      IEEEnoise['simTe']=TE
      IEEEnoise['Cy']=Cy
      IEEEnoise['Tout']=Tout
      return IEEEnoise

   def NoisePSO(self, Cpad1, Cpad2, Cpad12, Rsh1, Rsh2, Rsh3, Lsh1, Lsh2, Lsh3, Copen1, Copen2, Copen12,
               Rg, Rd, Rs, Ls, Ld, Lg, Cpg, Cpd,
               Gm, Rds, Cgd, Cgs, H11, H22, DF, Tout, Ri=0, Rgd=0,
               SWEEP="Start=100 MHz Stop=50 GHz Step=100 MHz",
               seed=None, debug=False):
      """
      what is called Gds is in reality 1/Gds, Gds is passed as Ohms! ->
      this bug has now been corrected, what was called Gds now is called Rds and
      this Rds has to be passed. in all OK flow, h stores Rds in a variable called 'Gds'
      return Sparams!
      """
      IEEEnoise=dict()
      if seed is None:
         val=repr(np.random.randint(9999999))
      else:
         val=repr(seed)

      netl=self.netlist_NoisePSO.format(Cpad1=Cpad1, Cpad2=Cpad2, Cpad12=Cpad12,
              Rsh1=Rsh1, Rsh2=Rsh2, Rsh3=Rsh3, Lsh1=Lsh1, Lsh2=Lsh2,  Lsh3=Lsh3,
              Copen1=Copen1, Copen2=Copen2, Copen12=Copen12,
              Rg=Rg,Rd=Rd,Rs=Rs,Lg=Lg,Ld=Ld,Ls=Ls, Cpg=Cpg, Cpd=Cpd,
              Gm=Gm, Rds=Rds, Cgd=Cgd, Cgs=Cgs, Ri=Ri, Rgd=Rgd,
              SWEEP=SWEEP, dev=val, DF=DF, Tout=Tout, H11=H11, H22=H22)
      fname="PSO_"+val+".cir"
      f = open(self.PathADS+fname, 'w')
      f.write(netl)
      f.close()
      #time.sleep(1) # may be needed on NFS or samba filesystems
      if debug : print(f"Ready to process {fname}")
      if not windows:
          resu=subprocess.Popen("hpeesofsim "+fname,shell=True,env=self.my_env,cwd=self.PathADS,stderr=subprocess.PIPE,stdout=subprocess.PIPE)
          resu.wait()
          stdout, stderr = resu.communicate()
          if debug: print(stdout)
          if len(stderr)>0: print(stderr)
      else:
          command=["hpeesofsim.exe", fname]
          resu=subprocess.run(command, shell=True, env=self.my_env, cwd=self.PathADS, stderr=subprocess.PIPE,
                         stdout=subprocess.PIPE)
          print(str(resu.stdout, 'utf-8').replace('"','').strip() )
      soutput = open(self.PathADS+val+"_PSO_Sparam2.dat", "r")
      noutput = open(self.PathADS+val+"_PSO_NFdata.dat", "r")
      coutput = open(self.PathADS+val+"_PSO_Cydata.dat", "r")
      f, Sp = rfdataanalysis.getSparamfromS2P(soutput)
      NFmin, Rn, Sopt, NF, TE = self.parseVarFileNoiseData(noutput)
      fc, Cy = rfdataanalysis.getSparamfromS2P(coutput)
      soutput.close()
      noutput.close()
      coutput.close()
      #return f, Sp, NFmin, Rn, Sopt, NF, TE, Cy
      IEEEnoise['F']=f
      IEEEnoise['simSp']=Sp
      IEEEnoise['NFmin']=NFmin
      IEEEnoise['Rn']=Rn
      IEEEnoise['Sopt']=Sopt
      IEEEnoise['simNF50']=NF
      IEEEnoise['simTe']=TE
      IEEEnoise['Cy']=Cy
      IEEEnoise['Tout']=Tout
      return IEEEnoise


   def NoiseOS(self, Cport1, Cport2,Cport12, Rsh1, Rsh2, Lsh1, Lsh2, 
               Rg, Rd, Rs, Ls, Ld, Lg, Cpg, Cpd,
               Gm, Rds, Cgd, Cgs, H11, H22, DF, Tout, Ri=0, Rgd=0,
               SWEEP="Start=100 MHz Stop=50 GHz Step=100 MHz",
               seed=None, debug=False):
      """
      what is called Gds is in reality 1/Gds, Gds is passed as Ohms! ->
      this bug has now been corrected, what was called Gds now is called Rds and
      this Rds has to be passed. in all OK flow, h stores Rds in a variable called 'Gds'
      return Sparams!
      """
      IEEEnoise=dict()
      if seed is None:
         val=repr(np.random.randint(9999999))
      else:
         val=repr(seed)

      netl=self.netlist_NoiseOS.format(Cport1=Cport1, Cport2=Cport2, 
              Cport12=Cport12, Rsh1=Rsh1, Rsh2=Rsh2, Lsh1=Lsh1, Lsh2=Lsh2,
              Rg=Rg,Rd=Rd,Rs=Rs,Lg=Lg,Ld=Ld,Ls=Ls, Cpg=Cpg, Cpd=Cpd,
              Gm=Gm, Rds=Rds, Cgd=Cgd, Cgs=Cgs, Ri=Ri, Rgd=Rgd,
              SWEEP=SWEEP, dev=val, DF=DF, Tout=Tout, H11=H11, H22=H22)
      fname="OS_"+val+".cir"
      f = open(self.PathADS+fname, 'w')
      f.write(netl)
      f.close()
      #time.sleep(1) # may be needed on NFS or samba filesystems
      if debug : print(f"Ready to process {fname}")
      if not windows:
          resu=subprocess.Popen("hpeesofsim "+fname,shell=True,env=self.my_env,cwd=self.PathADS,stderr=subprocess.PIPE,stdout=subprocess.PIPE)
          resu.wait()
          stdout, stderr = resu.communicate()
          if debug: print(stdout)
          if len(stderr)>0: print(stderr)
      else:
          command=["hpeesofsim.exe", fname]
          resu=subprocess.run(command, shell=True, env=self.my_env, cwd=self.PathADS, stderr=subprocess.PIPE,
                         stdout=subprocess.PIPE)
          print(str(resu.stdout, 'utf-8').replace('"','').strip() )
      soutput = open(self.PathADS+val+"_OS_Sparam2.dat", "r")
      noutput = open(self.PathADS+val+"_OS_NFdata.dat", "r")
      coutput = open(self.PathADS+val+"_OS_Cydata.dat", "r")
      f, Sp = rfdataanalysis.getSparamfromS2P(soutput)
      NFmin, Rn, Sopt, NF, TE = self.parseVarFileNoiseData(noutput)
      fc, Cy = rfdataanalysis.getSparamfromS2P(coutput)
      soutput.close()
      noutput.close()
      coutput.close()
      #return f, Sp, NFmin, Rn, Sopt, NF, TE, Cy
      IEEEnoise['F']=f
      IEEEnoise['simSp']=Sp
      IEEEnoise['NFmin']=NFmin
      IEEEnoise['Rn']=Rn
      IEEEnoise['Sopt']=Sopt
      IEEEnoise['simNF50']=NF
      IEEEnoise['simTe']=TE
      IEEEnoise['Cy']=Cy
      IEEEnoise['Tout']=Tout
      return IEEEnoise


   def SimSys(self, Ypad1, Ypad2, Ypad3, Yvia1, Yvia2, Yvia3, Ls1, Ls2, Ls3, Rs1, Rs2, Rs3, 
                    Cpd, Cpg, Rg, Rd, Rs, Lg, Ld, Ls, 
                    Gm, Gds, Cgd, Cgs, H11, H22, Tout, Ri, Rgd, seed=None, debug=False):
      if seed is None:
         val=repr(np.random.randint(9999999))
      else:
         val=repr(seed)
      netl=self.netlist_sys.format( Ypad1=Ypad1, Ypad2=Ypad2, Ypad3=Ypad3, Yvia1=Yvia1, Yvia2=Yvia2,
                   Yvia3=Yvia3, Ls1=Ls1, Ls2=Ls2, Ls3=Ls3, Rs1=Rs1, Rs2=Rs2,
                   Rs3=Rs3, Cpd=Cpd, Cpg=Cpg, Rg=Rg, Rd=Rd, Rs=Rs,Lg=Lg,
                   Ld=Ld, Ls=Ls, Gm=Gm, Gds=Gds, Cgd=Cgd, Cgs=Cgs, H11=H11,
                    H22=H22, Tout=Tout, Ri=Ri, Rgd=Rgd, DF=1, dev=val ) 
      fname="Spsys_"+val+".cir"
      f = open(self.PathADS+fname, 'w')
      f.write(netl)
      f.close()
      if debug : print(fname)
      resu=subprocess.Popen("hpeesofsim "+fname,shell=True,env=self.my_env,cwd=self.PathADS,stderr=subprocess.PIPE,stdout=subprocess.PIPE)
      resu.wait()
      stdout, stderr = resu.communicate()
      if debug: print(stdout)
      if len(stderr)>0: print(stderr)
      soutput = open(self.PathADS+val+"_Sparam2_sys.dat", "r")
      noutput = open(self.PathADS+val+"_NFdata_sys.dat", "r")
      coutput = open(self.PathADS+val+"_Cydata_sys.dat", "r")
      f, Sp = rfdataanalysis.getSparamfromS2P(soutput)
      NFmin, Rn, Sopt, NF, TE = self.parseVarFileNoiseData(noutput)
      fc, Cy = rfdataanalysis.getSparamfromS2P(coutput)
      soutput.close()
      noutput.close()
      coutput.close()
      return f, Sp, NFmin, Rn, Sopt, NF, TE, Cy

################################################################################
################################################################################
################################################################################

def touchstoneS2Psave(Fileandpath, Spar, Freq, comment='', fmt='RI', precision=10, toGHz=False):
    """
    This function save in touchstone format the Sparameters for a 2 Port or coplex impedance.
    Can add a comment for the first line of the file in the 
    optional parameter comment and can save in Hz or GHz.
    defalt fourmat is RI, but MA and DB are supported.
    """
    dig=precision # digits of precision
    col=dig+8 # allows +- sign, dot, and three digits exponent
    # decide for one or two port sweep:
    shape=np.shape(Spar)
    print(f"shape to save is {shape}")
    if len(np.shape(Spar)) ==1:
        TwoPort=False # it is a 1Port sweep
        if not type(Spar[0]) is np.complex_:
            print("WARNING: sweep is not complex values? ")
    elif len(np.shape(Spar)) ==3 and shape[-1]==2 and shape[-2]==2:
        TwoPort=True
    else:
        print(f"Fatal: expecting a complex array or array of 2x2 matrix.")
        return
    
    if fmt.upper()=='RI':
        cola= lambda z: np.real(z)
        colb= lambda z: np.imag(z)
        comm="! f ReS11 ImS11 ReS21 ImS21 ReS12 ImS12 ReS22 ImS22"
    elif fmt.upper()=='MA':
        cola= lambda z: np.abs(z)
        colb= lambda z: np.rad2deg(np.angle(z)) # tochstone is 
        print("Warning, pls double check, this is untested!")
        comm="! f absS11 degS11 absS21 degS21 absS12 degS12 absS22 degS22"
    elif fmt.upper()=='DB':
        cola= lambda z: 20*np.log10(np.abs(z)) # touchstone are 20*log10
        colb= lambda z: np.rad2deg(np.angle(z)) 
        print("Warning, pls double check, this is untested!")
        comm="! f dBS11 degS11 dBS21 degS21 dBS12 degS12 dBS22 degS22 dB is 20log10"
    else:
        cola= lambda z: np.real(z)
        colb= lambda z: np.imag(z)
        print(f"Warning: could not reognized number format {fmt}, defaulting to RI")
    
    
    if not (len(Freq)==len(Spar)):
        print(f"fatal: len(F)={len(Freq)} is not matching len(Spar)={len(Spar)}!")
        return

    with open(Fileandpath, "w") as f:
        if comment:
            print("! "+comment, file=f)
        print(comm, file=f)    
        if toGHz:
            print(f"# GHz S {fmt.upper()} R 50", file=f)
        else:
            print(f"# Hz S {fmt.upper()} R 50", file=f)
        for i, fval in enumerate(Freq):
            if toGHz:
                print(f"{fval*1e-9:<15e}", end=" ", file=f)
            else:
                print(f"{int(fval):<13}", end=" ", file=f)
            if TwoPort:
                print(f"{cola(Spar[i,0,0]):{col}.{dig}e}  {colb(Spar[i,0,0]):{col}.{dig}e}", end="  ", file=f)
                print(f"{cola(Spar[i,1,0]):{col}.{dig}e}  {colb(Spar[i,1,0]):{col}.{dig}e}", end="  ", file=f)
                print(f"{cola(Spar[i,0,1]):{col}.{dig}e}  {colb(Spar[i,0,1]):{col}.{dig}e}", end="  ", file=f)
                print(f"{cola(Spar[i,1,1]):{col}.{dig}e}  {colb(Spar[i,1,1]):{col}.{dig}e}", file=f)
            else:
                print(f"{cola(Spar[i]):{col}.{dig}e}  {colb(Spar[i]):{col}.{dig}e}", file=f)

def GenericMDIFSpar(Fileandpath, Spar, Freq, NF50=[], comment='', blockname='Sparams'):
    """
    Write first, document it later ...
    """
    YesNoise=False
    if len(NF50)> 0: YesNoise=True
    if not (len(Freq)==len(Spar)):
        print(f"fatal: len(F)={len(Freq)} is not matching len(Spar)={len(Spar)}!")
        return
    if YesNoise and not (len(Freq)==len(NF50)):
        print(f"fatal: Len(NF)={len(NF50)} is not matching len(F)={len(Freq)}!")
        return
    if not np.shape(Spar[0])==(2,2): 
        print(f"fatal: Spar doesn't seems to be a 2x2 matrix, but it is a {np.shape(Spar[0])} !")
        return
    with open(Fileandpath, "w") as f:
        if comment:   
            print("! "+comment, file=f)
        print(f"BEGIN {blockname}", file=f)
        if YesNoise:
            print("% Freq(real)  S[1,1](complex) S[2,1](complex) S[1,2](complex) S[2,2](complex) NF50(real)", file=f)
        else:
            print("% Freq(real)  S[1,1](complex) S[2,1](complex) S[1,2](complex) S[2,2](complex)", file=f)
        for i, fval in enumerate(Freq):
            print(f"{int(fval):<13}", end=" ", file=f)
            print(f"{np.real(Spar[i,0,0]):20.17e}  {np.imag(Spar[i,0,0]):20e}", end="  ", file=f)
            print(f"{np.real(Spar[i,1,0]):20.17e}  {np.imag(Spar[i,1,0]):20e}", end="  ", file=f)
            print(f"{np.real(Spar[i,0,1]):20e}  {np.imag(Spar[i,0,1]):20e}", end="  ", file=f)
            if YesNoise:
                print(f"{np.real(Spar[i,1,1]):20e}  {np.imag(Spar[i,1,1]):20e}  {NF50[i]:20e}", file=f)
            else:
                print(f"{np.real(Spar[i,1,1]):20e}  {np.imag(Spar[i,1,1]):20e}", file=f)
        print("END", file=f)

def ADS_GenericMDIF_SSexport(Fileandpath, SSobj, comment='', reH11=[], reH22=[], Tout=[], blockname='SSmodel'):
    """
    This routine saves all extractions that are bias dependent (i.e. the INTRINSIC)
    paramenter of a DUT from an Object SSobj (SS = small signal) and also the parameter necessary
    to compute the noise, i.e. H11, H22 and Tout.
    Saved format is a 'Generic MDIF' file that can be read easily in ADS.

    The routine tries to be smart in the sense that Ri, Tau of Rgd may be present or not, and
    are checked individually. 
    Also H11, H22 and Tout are optionally save. But in these case if either one of H11 and H22 is present, then all three must be 
    there and if Tout is not present it is defaulted to 1000K.
    The extrinsic part is not saved and it is supposed to be included manually
    """
    HasGm =hasattr(SSobj, "Gm")
    HasGds=hasattr(SSobj, "Gds")
    HasVg =hasattr(SSobj, "vg")
    HasCgs=hasattr(SSobj, "Cgs")
    HasCgd=hasattr(SSobj, "Cgd")
    HasCds=hasattr(SSobj, "Cds")
    HasRi =hasattr(SSobj, "Ri")
    HasTau=hasattr(SSobj, "Tau")
    HasRgd=hasattr(SSobj, "Rgd")
    HasVd =hasattr(SSobj, "vd")
    HasH11=True if len(reH11)>0 else False
    HasH22=True if len(reH22)>0 else False
    HasTout=True if len(Tout)>0 else False
    
    if not all([HasGm, HasVd, HasVg, HasGds, HasCgs, HasCgd, HasCds]):
        print(f"FATAL: passed object seems not a SS circuit with bias")
        return False
    # still the case the sweep is assign the NONE
    if SSobj.vd.all()==None:
        HasVd=False
        setVd=['0.0'] ## if SSobj.vd is none Vd defalts to 0.0 

    if HasVd:
        setVd=list(sorted(set(SSobj.vd)))
    setVg=list(sorted(set(SSobj.vg)))

    # if Tout is not passed, but H11 and H22 yes, then Tout defauts to 1000
    # to allow optimization
    if HasH11: 
        if not HasTout:
            Tout=[1000 for i in reH11]
            HasTout=True

    if HasH11: 
        if not all([HasH11, HasH22, HasTout]):
            print(f"FATAL, If one of H11, H22 or Tout is present, the other two must be there also!")
            return False
        if not len(reH11)==len(reH22)==len(Tout)==len(SSobj.vg):
            print(f"FATAL: bias sweep length is not matching!")
            return False
    ##### atually start doing something!

    head="% Vg(real) "
    head=head+" Gm(real) Gds(real) Cgs(real) Cgd(real) Cds(real) "
    if HasRi: head=head+" Ri(real)"
    if HasTau: head=head+" Tau(real)"
    if HasRgd: head=head+" Rgd(real)"
    if HasH11: head=head+" H11(real)  H22(real)  Tout(real)"

    with open(Fileandpath, "w") as f:
        if comment:
            print(f"! {comment}", file=f)
        for vd in setVd:
            if HasVd: print(f"var Vd(real)={vd}", file=f)
            print(f"BEGIN {blockname}", file=f)
            print(f"{head}", file=f)
            for vg in setVg:
                i=SSobj.findbias(vg=vg, vd=vd)
                line=f" {vg:<10}"
                line=line+f" {SSobj.Gm[i]:15} {SSobj.Gds[i]:15}"
                line=line+f" {SSobj.Cgs[i]:15} {SSobj.Cgd[i]:15}"
                line=line+f" {SSobj.Cds[i]:15}"
                if HasRi : line=line+f" {SSobj.Ri[i]:15}"
                if HasTau: line=line+f" {SSobj.Tau[i]:15}"
                if HasRgd: line=line+f" {SSobj.Rgd[i]:15}"
                if HasH11: line=line+f" {reH11[i]:15} {reH22[i]:15} {Tout[i]:15}"  
                print(f"{line}", file=f) 
            print(f"END", file=f)

def ADS_discrete_export(Fileandpath, SSobj, comment=''):
    """
    File is OK, format seems OK. But could not be used from within ADS. Reason not clear
    Write first, document later! ...
    """
    HasGm =hasattr(SSobj, "Gm")
    HasGds=hasattr(SSobj, "Gds")
    HasVg =hasattr(SSobj, "vg")
    HasCgs=hasattr(SSobj, "Cgs")
    HasCgd=hasattr(SSobj, "Cgd")
    HasCds=hasattr(SSobj, "Cds")
    HasRi =hasattr(SSobj, "Ri")
    HasTau=hasattr(SSobj, "Tau")
    HasRgd=hasattr(SSobj, "Rgd")
    HasVd =hasattr(SSobj, "vd")
    
    if not all([HasGm, HasVd, HasVg, HasGds, HasCgs, HasCgd, HasCds]):
        print(f"FATAL: passed object seems not a SS circuit with bias")
        return False
    # still the case the sweep is assign the NONE
    if SSobj.vd.all()==None:
        HasVd=False

    with open(Fileandpath, "w") as f:
        if comment:
            print(f"REM {comment}", file=f)
        ind=1
        print(f"BEGIN DSCRDATA", file=f)
        head="% INDEX Vg "
        if HasVd: head=head+" Vd "
        head=head+" Gm Gds Cgs Cgd Cds "
        if HasRi: head=head+" Ri"
        if HasTau: head=head+" Tau"
        if HasRgd: head=head+" Rgd"
        print(f"{head}", file=f)
        for i, vg in enumerate(SSobj.vg):
            line=""
            line=f"{ind:4}  {vg:10}"
            if HasVd: 
                line=line+f" {SSobj.vd[i]:10}"
            line=line+f" {SSobj.Gm[i]:15} {SSobj.Gds[i]:15}"
            line=line+f" {SSobj.Cgs[i]:15} {SSobj.Cgd[i]:15}"
            line=line+f" {SSobj.Cds[i]:15}"
            if HasRi : line=line+f" {SSobj.Ri[i]:15}"
            if HasTau: line=line+f" {SSobj.Tau[i]:15}"
            if HasRgd: line=line+f" {SSobj.Rgd[i]:15}"
            print(f"{line}", file=f) 
            ind=ind+1
        print(f"END DSCRDATA", file=f)


def GenericMDIFSparVorI(Fileandpath, Spar, Freq, varname, varsweep, NF50=[], comment='', blockname='Sparams'):
    """
    Dumps A Sparameter measurement sweep in a 'GenericMDIF' file format that
    is easily loaded by ADS with a data element
    it is suppose to write Sparams in a sweep, either in voltage or current.
    Fileandpath 
       is a string with the name of the file to be written
    Spar 
       is supposed to be a ndarray of [nb, nf, 2, 2], so 4 dimensions:
       first dimension is the bias sweep, second is the frequency sweeep
       then the 2x2 matrix at that bias and frequency.
    varname
       is the (float) variable name: Vg, Vanode, Ianode ...
    varsweep 
       is a one-dimensional array with the frequency sweeep values
    NF50 
       is (optional) and it may store NF50 measurements values
    """
    YesNoise=False
    if len(NF50)> 0: YesNoise=True
    nbias=len(varsweep)
    sSp=np.shape(Spar)
    sNf=np.shape(NF50)
    # some sanity checks on the parameter passed ...
    if not len(sSp)==4: 
        print(f"fatal: sparams must be a 4 dimesional array [bias][freq][2x2 matrix] !")
        return
    if YesNoise and not len(sNf)==2:
        print(f"fatal, NF50 must be a 2 dimensional array [bias][frequency]!")
        return
    if not nbias==sSp[0]:
        print(f"fatal, bias points are not matching! bias vector has  has {nbias} points, Sparams primary sweep is {sSp[0]} points.")
        return
    if YesNoise and not nbias==sNf[0]:
        print(f"fatal, bias points are not matching! NF50 has {sNf[0]} bias points, bias vector has {nbias} points")
        return
    if not (len(Freq)==sSp[1]):
        print(f"fatal: len(F)={len(Freq)} is not matching frequency sweep in Spar that has {sSp[1]} points!")
        return
    if YesNoise and not (len(Freq)==sNf[1]):
        print(f"fatal: frequency sweep in NF50 that has {sNf[1]} points is not matching len(F)={len(Freq)}!")
        return
    if not (sSp[2], sSp[3])==(2,2): 
        print(f"fatal: Spar doesn't seems to be a 2x2 matrix!")
        return
    # start doing 

    with open(Fileandpath, "w") as f:
        if comment:   
            print("! "+comment, file=f)
        for j, biasp in enumerate(varsweep): 
            print(f"var {varname}(real)={biasp}", file=f)
            print(f"BEGIN {blockname}", file=f)
            if YesNoise:
                print("% Freq(real)  S[1,1](complex) S[2,1](complex) S[1,2](complex) S[2,2](complex) NF50(real)", file=f)
            else:
                print("% Freq(real)  S[1,1](complex) S[2,1](complex) S[1,2](complex) S[2,2](complex)", file=f)
            for i, fval in enumerate(Freq):
                print(f"{int(fval):<13}", end=" ", file=f)
                print(f"{np.real(Spar[j,i,0,0]):15e}  {np.imag(Spar[j,i,0,0]):15e}", end="  ", file=f)
                print(f"{np.real(Spar[j,i,1,0]):15e}  {np.imag(Spar[j,i,1,0]):15e}", end="  ", file=f)
                print(f"{np.real(Spar[j,i,0,1]):15e}  {np.imag(Spar[j,i,0,1]):15e}", end="  ", file=f)
                if YesNoise:
                    print(f"{np.real(Spar[j,i,1,1]):15e}  {np.imag(Spar[j,i,1,1]):15e}  {NF50[j,i]:15e}", file=f)
                else:
                    print(f"{np.real(Spar[j,i,1,1]):15e}  {np.imag(Spar[j,i,1,1]):15e}", file=f)
            print("END", file=f)


def GenericMDIFSparVdVg(Fileandpath, Spar, Freq, vvd, vvg, NF50=[], comment='', blockname='Sparams'):
    """
    Dumps A Sparameter measurement sweep in a 'GenericMDIF' file format that
    is easily loaded by ADS with a data element
    """
    YesNoise=False
    if len(NF50)> 0: YesNoise=True
    nbias=len(vvg)
    if not len(vvg) == nbias:
        print(f"fatal: len(vvg)={len(vvg)} is not matching len(vvd)={len(vvd)}!")
        return
    sSp=np.shape(Spar)
    sNf=np.shape(NF50)
    if not len(sSp)==4: 
        print(f"fatal: sparams must be a 4 dimesional array [bias][freq][2x2 matrix] !")
        return
    if YesNoise and not len(sNf)==2:
        print(f"fatal, NF50 must be a 2 dimensional array [bias][frequency]!")
        return
    if not nbias==len(vvd)==sSp[0]:
        print(f"fatal, bias points are not matching! Vg has {len(vvg)} points, Vd {len(vvd)}, Sparams {sSp[0]}")
        return
    if YesNoise and not nbias==sNf[0]:
        print(f"fatal, bias points are not matching! NF50 has {sNf[0]} bias points, Vg has {len(vvg)}")
        return

    if not (len(Freq)==sSp[1]):
        print(f"fatal: len(F)={len(Freq)} is not matching frequency sweep in Spar that has {sSp[1]} points!")
        return
    if YesNoise and not (len(Freq)==sNf[1]):
        print(f"fatal: frequency sweep in NF50 that has {sNf[1]} points is not matching len(F)={len(Freq)}!")
        return

    if not (sSp[2], sSp[3])==(2,2): 
        print(f"fatal: Spar doesn't seems to be a 2x2 matrix!")
        return


    with open(Fileandpath, "w") as f:
        if comment:   
            print("! "+comment, file=f)

        for j, vg in enumerate(vvg): 
            vd=vvd[j]
            print(f"var Vg(real)={vg}", file=f)
            print(f"var Vd(real)={vd}", file=f)
            print(f"BEGIN {blockname}", file=f)
            if YesNoise:
                print("% Freq(real)  S[1,1](complex) S[2,1](complex) S[1,2](complex) S[2,2](complex) NF50(real)", file=f)
            else:
                print("% Freq(real)  S[1,1](complex) S[2,1](complex) S[1,2](complex) S[2,2](complex)", file=f)
            for i, fval in enumerate(Freq):
                print(f"{int(fval):<13}", end=" ", file=f)
                print(f"{np.real(Spar[j,i,0,0]):15e}  {np.imag(Spar[j,i,0,0]):15e}", end="  ", file=f)
                print(f"{np.real(Spar[j,i,1,0]):15e}  {np.imag(Spar[j,i,1,0]):15e}", end="  ", file=f)
                print(f"{np.real(Spar[j,i,0,1]):15e}  {np.imag(Spar[j,i,0,1]):15e}", end="  ", file=f)
                if YesNoise:
                    print(f"{np.real(Spar[j,i,1,1]):15e}  {np.imag(Spar[j,i,1,1]):15e}  {NF50[j,i]:15e}", file=f)
                else:
                    print(f"{np.real(Spar[j,i,1,1]):15e}  {np.imag(Spar[j,i,1,1]):15e}", file=f)
            print("END", file=f)


def ParseMDFSparams(filename, trans={"Freq" :  "F", "S1_1" :  "s11", "S1_2" :  "s12", "S2_1" :  "s21", "S2_2":"s22"}, debug=False, ReadFrequency=False):
    """
    This routine is meant to load Sparams saved by the PNAX
    Example usage:
       f, Sm = ParseMDFSparams(datadir+"/1stCal/sweep_Z60X_-10dBm.mdf", {"Freq":"F","S1_1":"s11", "S1_2":"s12", "S2_1":"s21", "S2_2":"s22"})
       
    it is possible to implemnt 'on the fly' conversions passing a 'Trans' dictionary: 
      {   "Freq" :  "F",
          "S1_1" :  "s11",
          "S1_2" :  "s12", 
          "S2_1" :  "s21",
          "S2_2":"s22"
          }
    in case of another port is used:

    MdifSp34={ "Freq":"F",
                  "S3_3":"s11", 
                  "S3_4":"s12", 
                  "S4_3":"s21", 
                  "S4_4":"s22",
                } # tradotti per essere processati come dati port1 port2 

    This routine is a wrapper around singleBlackMDFparse.
    That routin is low level and returns a hash. This one is for one line Sparam reading
    trans may be "12" or "34" to use the standard hashes:
    MdifSp12={"Freq":"F",
              "S1_1":"s11", 
              "S1_2":"s12", 
              "S2_1":"s21", 
              "S2_2":"s22"
              } # from  and to variable names
    or
    MdifSp34   ={"Freq":"F","S3_3":"s11", "S3_4":"s12", "S4_3":"s21", "S4_4":"s22"} # from  and to variable names
    tans may be adictionary, in that case it is used as that.                
    """

    MdifSp12={"Freq":"F",
              "S1_1":"s11", 
              "S1_2":"s12", 
              "S2_1":"s21", 
              "S2_2":"s22"
              }
    MdifSp34   ={"Freq":"F",
                "S3_3":"s11", 
                "S3_4":"s12", 
                "S4_3":"s21", 
                "S4_4":"s22"} # from  and to variable names
    if type(trans)==dict:
        uset=trans
    elif type(trans)==str:
        if trans=='12': 
            uset=MdifSp12
        elif trans=='34': 
            uset=MdifSp34
        else:
            print(f"Cannot interpret {trans}")
            uset=MdifSp12
    else:
        print(f"Cannot interpret trans parameter")
        uset=MdifSp12

    tmp= singleBlockMDFparse(filename=filename, trans=uset, debug=debug, SaveComments=False)
    freq=tmp['F']
    s11=tmp['s11']
    s21=tmp['s21']
    s12=tmp['s12']
    s22=tmp['s22']
    S=np.empty([len(freq), 2, 2],dtype=np.complex128)
    f=np.empty(len(freq),dtype=np.float_)
    for i,ff in enumerate(freq):
        f[i]=freq[i]
        S[i,0,0]=s11[i]
        S[i,0,1]=s12[i]
        S[i,1,0]=s21[i]
        S[i,1,1]=s22[i]
    if ReadFrequency:
        return f, S
    else:
        return S


def singleBlockMDFparse(filename, trans=dict(), debug=False, SaveComments=False):
    """
    This routine is meant to load files saved by the PNAX, it is a low level 
    routine that returns a list of ndarrays (reads in everything but does not
    try to reconstruct the matrices. For faster processing use ParseMDFSparams.

    Example usage:
       da = singleBlockMDFparse(datadir+"/1stCal/sweep_Z60X_-10dBm.mdf",
         trans = {"Freq":"F","S1_1":"s11", "S1_2":"s12", "S2_1":"s21", "S2_2":"s22"})
       comm, da = singleBlockMDFparse(datadir+"/1stCal/sweep_Z60X_-10dBm.mdf",
         trans = {"Freq":"F","S1_1":"s11", "S1_2":"s12", "S2_1":"s21", "S2_2":"s22"}, SaveComments=True)

    As shown above, it is possible to implemnt 'on the fly' conversions passing a 'Trans' dictionary: 
      {   "Freq" :  "F",
          "S1_1" :  "s11",
          "S1_2" :  "s12", 
          "S2_1" :  "s21",
          "S2_2":"s22"
          }

    this is how the PNAX is saving the data in the MDFs
    MdifSp12={"Freq":"F",
              "S1_1":"s11", 
              "S1_2":"s12", 
              "S2_1":"s21", 
              "S2_2":"s22"
              } # from  and to variable names
    in case the unration measurements are saved also:

    MdifSp12unr={ "Freq":"F",
                  "S1_1":"s11", 
                  "S1_2":"s12", 
                  "S2_1":"s21", 
                  "S2_2":"s22",
                  "A,1":"a", 
                  "B,2":"b", 
                  "R1,1":"r1", 
                  "R2,2":"r2", 
                  "A,2":"ap", 
                  "B,1":"bp"
                } # ap is a' prime. bp is b'
    Without noise most of the time is better to use port 3 and 4

    MdifSp34unr={ "Freq":"F",
                  "S3_3":"s11", 
                  "S3_4":"s12", 
                  "S4_3":"s21", 
                  "S4_4":"s22",
                  "C_3":"a", 
                  "D,4":"b", 
                  "R3,3":"r1", 
                  "R4,4":"r2",
                  "C,4":"ap", 
                  "D,3":"bp"
                } # tradotti per essere processati come dati port1 port2 
    
    # from and to variable names for 2poer Sp on ports 3 and 4
    MdifSp34   ={"Freq":"F","S3_3":"s11", "S3_4":"s12", "S4_3":"s21", "S4_4":"s22"} 
    # same with unratio measureemnts, ap is a' prime. bp is b'
    MdifSp34unr={"Freq":"F","S3_3":"s11", "S3_4":"s12", "S4_3":"s21", "S4_4":"s22",
                 "C_3":"a", "D,4":"b", "R3,3":"r1", "R4,4":"r2","C,4":"ap", "D,3":"bp"}

    # from and to variable names for 2poer Sp on ports 1 and 2
    MdifSp12   ={"Freq":"F","S1_1":"s11", "S1_2":"s12", "S2_1":"s21", "S2_2":"s22"}
    # Same with unrationed measurements: ap is a' prime. bp is b'
    MdifSp12unr={"Freq":"F","S1_1":"s11", "S1_2":"s12", "S2_1":"s21", "S2_2":"s22",
                 "A,1":"a", "B,2":"b", "R1,1":"r1", "R2,2":"r2", "A,2":"ap", "B,1":"bp"}

    Notes from doc:
    All analyzer models have test port receivers and reference receivers. 
    For 4-port models: 
    R1, R2, R3, and R4 are reference receivers, they measure the signal as it leaves 
                              the analyzer source.
    R1 measures the signal out of Port 1
    R4 measures the signal out of Port 4

    A, B, C, and D are test port receivers. They measure the signal out 
    (or reflecting off ) of the DUT.

    A measures the signal into VNA Port 1
    B measures the signal into VNA Port 2
    C measures the signal into VNA Port 3
    D measures the signal into VNA Port 4
  
    A comma is added to specify the ACTIVE port:
    R1,1 is signal out of port1 when Port1 is active. Usually only R1,1 R2,2 R3,3 and R4,4 are meaningful.

    A,1 is signal received in port 1 when port 1 is active, so S11=A,1/R1,1
    B,1 is signal received in port 2 when port 1 is active, so S21=B,1/R1,1


    """

    data=dict()
    comments=list()
    nocom=list()
    headers=list()
    
    # parse the file in one go
    with open(filename, encoding='Latin1') as fin:
        lines=fin.read().splitlines()
        
    # find comments, start, stop and list of variables in the header
    for i,line in enumerate(lines): 
        #print(line)
        if len(line)<1: continue # drops empty lines
        if line[0]=='!':            # comment lines start with "!"
            comments.append(line)
            continue
        b=re.search("^% (.*)", line) # header line starts with "%"
        if b:  # 
            headers=b.group(1).split()
        start=re.search("^BEGIN (.*)",line)
        stop=re.search("^END(.*)", line)
        if start:
            blockname=start.group(1)
            if debug: print(f"blockname is {blockname}")
            st=i
        if stop:
            sp=i
    if debug: print(f"start {st} and stop {sp}")
        
    # create the hash and parse the variables  
    varnames=list()
    vartypes=list()
    for vname in headers:
        if debug: print(vname)
        test=re.search("^(\S+)\((\w+)\)", vname)
        if test:
            varname=(test.group(1))
            tmp=test.group(2)
            if tmp=="int" or tmp=="0":
                vartype="I"
            elif tmp=="real" or tmp=="1":
                vartype="F"
            elif tmp=="string" or tmp=="2":
                print("string type not implemented yet")
                return False
            elif tmp=="complex" or tmp=="3":
                vartype="C"
            elif tmp=="boolean" or tmp=="4":
                print("boolean type not implemented yet")
                return False
            elif tmp=="binary" or tmp=="5":
                print("binary type not implemented yet")
                return False
            elif tmp=="octal" or tmp=="6":
                print("octal type not implemented yet")
                return False
            elif tmp=="hexadecimal" or tmp=="7":
                print("hex type not implemented yet")
                return False
            elif tmp=="byte16" or tmp=="8":
                print("byte16 type not implemented yet")
                return False
            else:
                print(f"cannot interpret datatype '{tmp}'")
                return False
            if varname=="Freq": vartype="I"
            if varname in trans:
                varname=trans[varname]
            while varname in data:
                varname=varname+"_1"
            data[varname]=list()
            
            varnames.append(varname)
            vartypes.append(vartype)
    if debug: print(varnames)
    if debug: print(vartypes)
        
    # actually parse the data
    for i in range(st+2,sp):       # skipt begin line and header line
        dataline=lines[i].split()
        j=0
        for i,var in enumerate(varnames):
            if vartypes[i]=="I":
                data[var].append(int(float(dataline[j]))) # will do a floor of that data 
                j+=1
            if vartypes[i]=="F":
                data[var].append(float(dataline[j]))
                j+=1
            if vartypes[i]=="C":
                data[var].append(complex(float(dataline[j]), float(dataline[j+1]) ) )
                j+=2

    # last step, convert lists into np.arrays
    for varname in data.keys():
        data[varname]=np.array(data[varname])
    if SaveComments:
        return comments, data
    else:
        print(f"{comments}")
        return data

def MDFparseHeader(headers=list(), trans=dict(), debug=False):
    # create the hash and parse the variables  
    # Header is in the form: 
    # % Freq(real) NF(complex) S1_1(complex) S2_2(complex) S2_1(complex)
    # the % has been dropped and the header has been splitted in a list
    # each variable is varname(vartype)

    varnames=list()
    vartypes=list()
    data=dict()
    for vname in headers:
        if debug: print(vname)
        test=re.search("^(\S+)\((\w+)\)", vname)  # try to split vname in some1(some2)
        if test:
            varname=(test.group(1))
            tmp=test.group(2)
            if tmp=="int" or tmp=="0":
                vartype="I"
            elif tmp=="real" or tmp=="1":
                vartype="F"
            elif tmp=="string" or tmp=="2":
                print("string type not implemented yet")
                return list(), list(), dict()
            elif tmp=="complex" or tmp=="3":
                vartype="C"
            elif tmp=="boolean" or tmp=="4":
                print("boolean type not implemented yet")
                return list(), list(), dict()
            elif tmp=="binary" or tmp=="5":
                print("binary type not implemented yet")
                return list(), list(), dict()
            elif tmp=="octal" or tmp=="6":
                print("octal type not implemented yet")
                return list(), list(), dict()
            elif tmp=="hexadecimal" or tmp=="7":
                print("hex type not implemented yet")
                return list(), list(), dict()
            elif tmp=="byte16" or tmp=="8":
                print("byte16 type not implemented yet")
                return list(), list(), dict()
            else:
                print(f"cannot interpret datatype '{tmp}'")
                return list(), list(), dict()
      
            # fix some idiosincracies"
            if varname=="Freq": vartype="I"
            if varname in trans:
                varname=trans[varname]
            while varname in data:
                varname=varname+"_1"
            data[varname]=list()
            
            varnames.append(varname)
            vartypes.append(vartype)
        else:
            print("WARNING: couldn't parse {vname}!")
    if debug: print(varnames)
    if debug: print(vartypes)
    return varnames, vartypes, data

def SnPparse(filename="", lines=[], SaveComments=False, Verbose=False, ReadFrequency=False):
    """
    quick routine to load/imoprt/read touchstone Touchstone files saved from wincal.
    Won't support all touchstone file format features
    In principle all kind of matrix (S, Y, Z can be stored). Only Sparams are supported
    only whole line comments are allowed
    only 1 2 3 and 4 port are supported. 5 port and above are not supported

    use filename="namefile.s2p" to OPEN the file
    or use lines to pass all lines, e.g. if you unpacked from a zipfile
    """
    
    # parse the file in one go
    if not filename=="": 
    	with open(filename) as fin:
        	lines=fin.read().splitlines()

    header=list()    # usually "# Hz S RI R 50" 
    comments=list()  # usually at the header of the file, line starting with '!'
    datalines=list() # pure data

    # find comments, start, stop and list of variables in the header
    for i,line in enumerate(lines): 
        #print(line)
        if len(line)<1: continue # drops empty lines
        if line[0]=='!':            # comment lines start with "!"
            stripl=line[1:].strip()  # cancella leading '!' and trailing spaces
            if len( stripl )<1: continue # drops empty comments lines
            comments.append( stripl )
            continue
        b=re.search("^# (.*)", line) # header line starts with "#" 
        if b:  # found a header files
            header=b.group(1).upper().split() # Split the header into a list ['HZ', 'S', 'RI', 'R', '50']
            continue
        datalines.append(line)
    
    fact=header[0] # usually 'HZ', has been made uppercase!

    if header[1]!='S':
        print(f"Warning, expectiong 'S' for Sparams, but getting {header[0]}")
    if header[3]!='R':
        print(f"Warning, Header is not conforming to Touchstone!")
    if int(float(header[4]))!=50:
        print(f"Warning, parsing non-50-Ohms data? {header[4]} Ohms!\n{header}")

    if fact=='HZ':
        ff=1.0
    elif fact=='GHZ':
        ff=1e9
    elif fact=='MHZ':
        ff=1e6
    elif fact=='KHZ':
        ff=1e3
    elif fact=='THZ':
        ff=1e12
    else:
        print(f"Warning: cannot make sense of frequency units specified as {fact} !")
        ff=1.1

    if header[2]=='RI':
        string2num = lambda r,i: (float(r)+1j*float(i))
    elif header[2]=='MA':
        string2num = lambda m,a: float(m)*np.exp(1j*np.deg2rad(float(a)))
        print("Warning, pls double check, this is untested!")
    elif header[2]=='DB':
        string2num = lambda db,a: (10**(float(db)/20))*np.exp(1j*np.deg2rad(float(a)))
        print("Warning, pls double check, this is untested!")
    else:
        string2num = lambda r,i: (float(r)+1j*float(i))
        print(f"Warning: could not reognized number format {header[2]}, defaulting to RI")


    # now we have to decide if it is 1, 2, 3 or nn ports, and, if it is a 2 port data, if there is noise
    # data attached at the bottom:
    # we look at the data in the first 4 data lines and see if we can match a pattern:
    # 1 port: always 3 columns
    # fre re im  
    # 2 port: always 9 columns
    # fre s11r s11i s21r s21i s12r s12i s22r s22i 
    # 3 port: 7-6-6 data per line
    # freq s11 s12 s13
    # s21 s22 s23
    # s31 s32 s33
    # 4 port 9-8-8-8
    # freq s11 s12 s13 s14
    # s21 s22 s23 s24
    # s31 s32 s33 s34
    # s41 s42 s43 s44

    dt = len(datalines)
    line0=len( datalines[0].strip().split() )
    line1=len( datalines[1].strip().split() )
    line2=len( datalines[2].strip().split() )
    line3=len( datalines[3].strip().split() )

    
    if line0==3 and line1==3: # 1 Port data: freq real imag
        scatt=np.zeros(dt, dtype=np.complex_)
        freq=np.zeros(dt, dtype=np.float_)
        for i, line in enumerate(datalines):
            f, a, b = line.strip().split()
            freq[i]=float(f)
            scatt[i]=string2num(a,b)
    elif line0==9 and line1==9 : # 2 port data
        scatt=np.zeros((dt,2,2), dtype=np.complex_)
        freq=np.zeros(dt, dtype=np.float_)
        for i, line in enumerate(datalines):
            f, s11a, s11b, s21a, s21b, s12a, s12b, s22a, s22b = line.strip().split() # notare la famosa inversione S21 - S12
            freq[i]=float(f)
            scatt[i,0,0]=string2num(s11a,s11b)
            scatt[i,0,1]=string2num(s12a,s12b)
            scatt[i,1,0]=string2num(s21a,s21b)
            scatt[i,1,1]=string2num(s22a,s22b)
    elif line0==7 and line1==6 and line2==6 and line3==7 : # 3 port data
        # test that the lines are a multiple of three
        nfreq=dt//3 # integre division
        if dt%3>0: # remainder of the division
            print("Fatal, suppoesed to process 3 port data file but data lines are not multiple of 3")
            return
        scatt=np.zeros((nfreq,3,3), dtype=np.complex_)
        freq=np.zeros(nfreq, dtype=np.float_)
        for i in range(nfreq):
            f, s11a, s11b, s12a, s12b, s13a, s13b = datalines[3*i].strip().split() 
            s21a, s21b, s22a, s22b, s23a, s23b = datalines[3*i+1].strip().split() 
            s31a, s31b, s32a, s32b, s33a, s33b = datalines[3*i+2].strip().split() 
            freq[i]=float(f)
            scatt[i,0,0]=string2num(s11a,s11b)
            scatt[i,0,1]=string2num(s12a,s12b)
            scatt[i,0,2]=string2num(s13a,s13b)
            scatt[i,1,0]=string2num(s21a,s21b)
            scatt[i,1,1]=string2num(s22a,s22b)
            scatt[i,1,2]=string2num(s23a,s23b)
            scatt[i,2,0]=string2num(s31a,s31b)
            scatt[i,2,1]=string2num(s32a,s32b)
            scatt[i,2,2]=string2num(s33a,s33b)
    elif line0==9 and line1==8 and line2==8 and line3==8 : # 4 port data
        # test that the lines are a multiple of three
        nfreq=dt//4 # integre division
        if dt%4>0: # remainder of the division
            print("Fatal, suppoesed to process 4 port data file but data lines are not multiple of 3")
            return
        scatt=np.zeros((nfreq,4,4), dtype=np.complex_)
        freq=np.zeros(nfreq, dtype=np.float_)
        for i in range(nfreq):
            f, s11a, s11b, s12a, s12b, s13a, s13b, s14a, s14b = datalines[4*i].strip().split() 
            s21a, s21b, s22a, s22b, s23a, s23b, s24a, s24b = datalines[4*i+1].strip().split() 
            s31a, s31b, s32a, s32b, s33a, s33b, s34a, s34b = datalines[4*i+2].strip().split() 
            s41a, s41b, s42a, s42b, s43a, s43b, s44a, s44b = datalines[4*i+3].strip().split() 
            freq[i]=float(f)
            scatt[i,0,0]=string2num(s11a,s11b)
            scatt[i,0,1]=string2num(s12a,s12b)
            scatt[i,0,2]=string2num(s13a,s13b)
            scatt[i,0,3]=string2num(s14a,s14b)            
            scatt[i,1,0]=string2num(s21a,s21b)
            scatt[i,1,1]=string2num(s22a,s22b)
            scatt[i,1,2]=string2num(s23a,s23b)
            scatt[i,1,3]=string2num(s24a,s24b)
            scatt[i,2,0]=string2num(s31a,s31b)
            scatt[i,2,1]=string2num(s32a,s32b)
            scatt[i,2,2]=string2num(s33a,s33b)            
            scatt[i,2,3]=string2num(s34a,s34b)                        
            scatt[i,3,0]=string2num(s41a,s41b)
            scatt[i,3,1]=string2num(s42a,s42b)
            scatt[i,3,2]=string2num(s43a,s43b)            
            scatt[i,3,3]=string2num(s44a,s44b)                        
    else:
        print("Fatal, only 1 port and 2 port are supported!")
        return
    freq=ff*freq # fix 2023/03   
    if SaveComments: 
        if ReadFrequency:
            return (freq, scatt, comments)
        else:
            return (scatt, comments)
    else:
        if Verbose:
             print(comments)
        if ReadFrequency:
            return freq, scatt
        else:
            return scatt

def WincalMDFparse(filename, Verbose=False, SaveComments=False, ReadFrequency=False):
    """
    quick routine to load files saved from wincal.
    There are many difference between wincal mdf and PNAX mdf.
    Better have two different parsers
    """
    # parse the file in one go
    with open(filename) as fin:
        lines=fin.read().splitlines()

    header=list()    # usually "# Hz S RI R 50" 
    comments=list()  # usually at the header of the file, line starting with '!'
    datalines=list() # pure data

    # find comments, start, stop and list of variables in the header
    for i,line in enumerate(lines): 
        #print(line)
        if len(line)<1: continue # drops empty lines
        if line[0]=='!':            # comment lines start with "!"
            stripl=line[1:].strip()  # cancella leading '!' and trailing spaces
            if len( stripl )<1: continue # drops empty comments lines
            comments.append( stripl )
            continue
        if line[0:3]=='VAR':            # comment lines start also with VAR
            stripl=line[3:].strip()  # cancella leading '!' and trailing spaces
            if len( stripl )<1: continue # drops empty comments lines
            comments.append( stripl )
            continue
        if line[0]=='%':            # not a comment but a useless header
            stripl=line[1:].strip()  # cancella leading '!' and trailing spaces
            if len( stripl )<1: continue # drops empty comments lines
            comments.append( stripl )
            continue
 
        if line[0:5]=='BEGIN': continue
        if line[0:3]=='END': continue
        b=re.search("^# (.*)", line) # header line starts with "#" 
        if b:  # found a header files
            header=b.group(1).upper().split() # Split the header into a list ['HZ', 'S', 'RI', 'R', '50']
            continue
        datalines.append(line)

    fact=header[0] # usually 'HZ', has been made uppercase!

    if header[1]!='S':
        print(f"Warning, expectiong 'S' for Sparams, but getting {header[0]}")
    if header[3]!='R':
        print(f"Warning, Header is not conforming to TOuchstone!")
    if header[4]!='50':
        print(f"Warning, parsing non-50-Ohms data? {Header[4]} Ohms!")

    if fact=='HZ':
        ff=1.0
    elif fact=='GHZ':
        ff=1e9
    elif fact=='MHZ':
        ff=1e6
    elif fact=='KHZ':
        ff=1e3
    elif fact=='THZ':
        ff=1e12
    else:
        print(f"Warning: cannot make sense of frequency units specified as {fact} !")
        ff=1

    if header[2]=='RI':
        string2num = lambda r,i: (float(r)+1j*float(i))
    elif header[2]=='MA':
        string2num = lambda m,a: float(m)*np.exp(1j*np.deg2rad(float(a)))
        print("Warning, pls double check, this is untested!")
    elif header[2]=='DB':
        string2num = lambda db,a: (10**(float(db)/20))*np.exp(1j*np.deg2rad(float(a)))
        print("Warning, pls double check, this is untested!")
    else:
        string2num = lambda r,i: (float(r)+1j*float(i))
        print(f"Warning: could not reognized number format {header[2]}, defaulting to RI")


    # now we have to decide if it is 1, 2, 3 or nn ports, and, if it is a 2 port data, if there is noise
    # data attached at the bottom:

    dt = len(datalines)
    line0=datalines[0].strip().split()

    freq=np.zeros(dt, dtype=np.float_)
    if len(line0)==3: # 1 Port data: freq real imag
        scatt=np.zeros(dt, dtype=np.complex_)
        for i, line in enumerate(datalines):
            f, a, b = line.strip().split()
            freq[i]=float(f)
            scatt[i]=string2num(a,b)
    elif len(line0)==9: # 2 port data
        scatt=np.zeros((dt,2,2), dtype=np.complex_)
        for i, line in enumerate(datalines):
            #f, s11a, s11b, s21a, s21b, s12a, s12b, s22a, s22b = line.strip().split() # notare la famosa inversione S21 - S12
            # se il file e' salvato in mdf, l'ordine tre S12 e S21 è  naturale
            f, s11a, s11b, s12a, s12b, s21a, s21b, s22a, s22b = line.strip().split() # notare la famosa inversione S21 - S12
            freq[i]=float(f)
            scatt[i,0,0]=string2num(s11a,s11b)
            scatt[i,0,1]=string2num(s12a,s12b)
            scatt[i,1,0]=string2num(s21a,s21b)
            scatt[i,1,1]=string2num(s22a,s22b)
    else:
        print("Fatal, only 1 port and 2 port are supported!")
        return
    

    if SaveComments: 
        if ReadFrequency:
            return (freq, scatt, comments)
        else:
            return (scatt, comments)
    else:
        if Verbose:
            print(comments)
        if ReadFrequency:
            return freq, scatt
        else:
            return scatt


def MDFparse(filename, trans={"Freq":"F","S1_1":"s11", "S1_2":"s12", "S2_1":"s21", "S2_2":"s22","A_1":"a", "B_2":"b", "R1_1":"r1", "R2_2":"r2", "A_2":"ap", "B_1":"bp"}, 
        forceRealNF=True, debug=False, SaveComments=False):
    """
    This routine is meant to load files saved by the PNAX
    Example usage:
       da = MDFparse(datadir+"/1stCal/sweep_Z60X_-10dBm.mdf", {"Freq":"F","S1_1":"s11", "S1_2":"s12", "S2_1":"s21", "S2_2":"s22"})
       comm, da = MDFparse(datadir+"/1stCal/sweep_Z60X_-10dBm.mdf", {"Freq":"F","S1_1":"s11", "S1_2":"s12", "S2_1":"s21", "S2_2":"s22"}, SaveComments=True)
    it is possible to implemnt 'on the fly' conversions passing a 'Trans' dictionary: 
    this is how the PNAX is saving the data in the MDFs

    MdifSp12={  "Freq"  :"F",
                  "S1_1":"s11", 
                  "S1_2":"s12", 
                  "S2_1":"s21", 
                  "S2_2":"s22",
                  "A,1" :"a", 
                  "B,2" :"b", 
                  "R1,1":"r1", 
                  "R2,2":"r2", 
                  "A,2" :"ap", 
                  "B,1" :"bp"
                }               # ap is a' prime. bp is b'
    Without noise most of the time is better to use port 3 and 4

    MdifSp34unr={ "Freq":"F",
                  "S3_3":"s11", 
                  "S3_4":"s12", 
                  "S4_3":"s21", 
                  "S4_4":"s22",
                  "C_3" :"a", 
                  "D,4" :"b", 
                  "R3,3":"r1", 
                  "R4,4":"r2",
                  "C,4" :"ap", 
                  "D,3" :"bp"
                } # tradotti per essere processati come dati port1 port2 
    

    MdifSp34   ={"Freq":"F","S3_3":"s11", "S3_4":"s12", "S4_3":"s21", "S4_4":"s22"} # from  and to variable names
    MdifSp34unr={"Freq":"F","S3_3":"s11", "S3_4":"s12", "S4_3":"s21", "S4_4":"s22","C_3":"a", "D_4":"b", "R3_3":"r1", "R4_4":"r2","C_4":"ap", "D_3":"bp"} # ap is a' prime. bp is b'
    MdifSp12   ={"Freq":"F","S1_1":"s11", "S1_2":"s12", "S2_1":"s21", "S2_2":"s22"} # from  and to variable names
    MdifSp12unr={"Freq":"F","S1_1":"s11", "S1_2":"s12", "S2_1":"s21", "S2_2":"s22","A_1":"a", "B_2":"b", "R1_1":"r1", "R2_2":"r2", "A_2":"ap", "B_1":"bp"} # ap is a' prime. bp is b'
    """

    datablocks=dict()
    comments=list()
    nocom=list()
    headers=list() # list of lists!
    starts=list()   # list of starting blocks
    stops=list()    # list of stoping blocks
    blocks=list()

    # parse the file in one go
    if debug: print("read file in one go...")
    with open(filename, encoding='Latin1') as fin:
        lines=fin.read().splitlines()
    if debug: print("done.\nparsing file")        
    # find comments, start, stop and list of variables in the header
    for i,line in enumerate(lines): 
        #print(line)
        if len(line)<1: continue # drops empty lines
        if line[0]=='!':            # comment lines start with "!"
            comments.append(line)
            continue
        b=re.search("^% (.*)", line) # header line starts with "%"
        if b:  # found a header files
            headers.append( b.group(1).split() ) # Trick adds a list to the end of a list of lists
        start=re.search("^BEGIN (.*)",line)
        stop =re.search("^END(.*)", line)
        if start:
            blockname=start.group(1)  # prende il (.*) della search 
            #if debug: print(f"found a block {blockname}")
            starts.append(i)
            blocks.append(blockname)

        if stop:
            stops.append(i)
    if len(starts) != len (stops) or len(starts) != len (blocks) or  len(stops) != len (blocks) : 
        print(f"problems in identifying blocks: starts:[{starts}], stops:[{stops}], blocks:[{blocks}]")
        return

    if debug:
        print(f"identified blocks: starts:[{starts}], stops:[{stops}], starts:[{blocks}]")
        
    # create the hash and parse the variables  

    for b, bname in enumerate(blocks): # for each block
        varnames, vartypes, data=MDFparseHeader(headers=headers[b], trans=trans, debug=debug)
        if debug: print(f"processing {bname}...")
        # actually parse the data
        for l in range(starts[b]+2,stops[b]):       # skip begin line and header line i is line number in original file
            dataline=lines[l].split()   # split that line
            j=0  # dummy counter
            for i,var in enumerate(varnames):
                if vartypes[i]=="I":
                    data[var].append(int(float(dataline[j]))) # will do a floor of that data 
                    j+=1
                if vartypes[i]=="F":
                    data[var].append(float(dataline[j]))
                    j+=1
                if vartypes[i]=="C":
                    data[var].append(complex(float(dataline[j]), float(dataline[j+1]) ) )
                    j+=2
        # last step, convert lists into np.arrays
        for varname in data.keys():
            data[varname]=np.array(data[varname])
            if varname=='NF' and forceRealNF: 
                data[varname]=np.real(data[varname])
        datablocks[bname]=data.copy()   
    if SaveComments:
        return comments, datablocks
    else:
        print(f"{comments}")
    return datablocks


################################################################################
################################################################################
################################################################################
### SMALL CIRCUIT EXTRACTION
### need fix
### asking a single value:
###        self.defs['Gm']=rft.findHz(self.F, 5e9, 5e9)  # for Tau
### should yield the value at that exact value.
### instead there is a division by zero ssince in that case
### fstart==fstop
### and anything[fstart:fstop]=[]
### and then also len()=0

class SS_trivial(object):
    def __init__(self, Freq, Ym, Vg, Vd=None):
        self.Y=Ym
        self.vg=np.array(Vg)
        if len(Vd)==len(Vg):
            self.vd=np.array(Vd)
        elif len(Vd==1):
            self.vd=Vd*np.ones(np.shape(self.vg))
        elif Vd==None:
            self.vd=None
        else:
            print("Warning, Vg and Vg are not matching dimensions!")
            self.Vd=None
        
        self.F=np.array(Freq)  # Hz
        self.FHz = self.F
        self.FGHz = self.F/1e9 # GHz
        self.omega=2*np.pi*self.F
        self.deftable=dict() # this stores only start and stop limits
        self.deftable['C']=(10e9, 20e9)  # most consistent with Ousmane 20GHz to 40GHz
        self.deftable['R']=( 5e9, 15e9)
        self.deftable['G']=(10e9, 20e9)  # most consistent with Ousmane 20GHz to 40GHz
        self.deftable['F']=( 3e9, 10e9)  # for Ft. Fmax
        self.deftable['T']=(20e9, 40e9)  # for Tau
        

        self.defs=dict()
        #self.defs['C']=rft.findHz(self.F, 10e9, 20e9)  # most consistent with Ousmane 20GHz to 40GHz
        #self.defs['R']=rft.findHz(self.F,  5e9, 15e9)
        #self.defs['G']=rft.findHz(self.F, 10e9, 20e9)  # most consistent with Ousmane 20GHz to 40GHz
        #self.defs['F']=rft.findHz(self.F,  3e9, 10e9) # for Ft. Fmax
        #self.defs['T']=rft.findHz(self.F, 20e9, 40e9)  # for Tau
        self.sync_defs()

        self.Y11=self.Y[:,:,0,0]
        self.Y22=self.Y[:,:,1,1]
        self.Y12=self.Y[:,:,0,1]
        self.Y21=self.Y[:,:,1,0]
        self.Z=uW.Y2Z(self.Y)
        self.S=uW.Z2S(self.Z)

        self.Z11=self.Z[:,:,0,0]
        self.Z22=self.Z[:,:,1,1]
        self.Z12=self.Z[:,:,0,1]
        self.Z21=self.Z[:,:,1,0]
        
        self.S11=self.S[:,:,0,0]
        self.S22=self.S[:,:,1,1]
        self.S12=self.S[:,:,0,1]
        self.S21=self.S[:,:,1,0]        

    def findbias(self, vg=0.0, vd=0.0): # will work with or without BG
        for i, vvg in enumerate(self.vg):    # for each bias point
            if vvg==vg and self.vd[i]==vd:
                return i
        raise NameError('No Bias point found')

    def sync_defs(self):
        self.defs=dict()
        for name,element in self.deftable.items():
            self.defs[name]=rft.findHz(self.F,element[0], element[1])

    def set_span(self, var, start, stop):
        """  
        this methods set the dafault span to extract the 
        property defined in "var" as a string (upper case converted)
        from the start and stop values expressed in Hz
        """
        self.deftable[var.upper()]=(start, stop)
        self.sync_defs()

    def get_span(self, var):
        testvar=var.upper()
        if testvar in self.defs:
            return self.defs[testvar]
        elif testvar[0] in self.defs:
            return self.defs[testvar[0]]
        else:
            print(f"Warning, could not identify {var} as a SS property")
            return (0, len(self.F)) # by default return all sweep

    @property
    def Gmv(self):
        return np.real(self.Y21-self.Y12)
    @property
    def Gdsv(self):
        return np.real(self.Y22)
    @property
    def Cgsv(self):
        return np.imag(self.Y11+self.Y12)/self.omega
    @property
    def Cgdv(self):
        return -np.imag(self.Y12)/(self.omega)
    @property
    def Cdsv(self):
        return np.imag(self.Y22+self.Y12)/self.omega
    @property
    def Ftv(self):
        num=self.Y21
        den=self.Y11
        h21=num/den
        ftv=(self.F*np.abs(h21))
        return ftv
    @property
    def Fmaxv(self):
        num=np.abs(self.Y21-self.Y12)**2
        den=4*(self.Y11.real*self.Y22.real-self.Y12.real*self.Y21.real)
        fmaxv=(self.FHz*np.sqrt(np.abs(num/den)))    
        return fmaxv
    @property
    def Rolletk(self):
        num=2.0*np.real(self.Y11)*np.real(self.Y22)
        num=num-np.real(self.Y12*self.Y21)
        Rolletk=num/np.abs(self.Y12*self.Y21)
        return Rolletk
    @property
    def kstability(self):
    # taken from rft.S2Pstability(S):
        S11=self.S11
        S22=self.S22
        S12=self.S12
        S21=self.S21
        Det=S11*S22-S12*S21
        kstb=1-np.abs(S11)**2-np.abs(S22)**2+np.abs(Det)**2
        kstb=kstb/(2*np.abs(S12*S21))
        return kstb
    @property
    def UMasondB(self):
        S11=self.S11
        S22=self.S22
        S12=self.S12
        S21=self.S21
        Ratio=S21/S12
        Det=S11*S22-S12*S21
        kstb=1-np.abs(S11)**2-np.abs(S22)**2+np.abs(Det)**2
        kstb=kstb/(2*np.abs(S12*S21))
        mason = (np.abs(Ratio-1)**2)/(2*kstb*np.abs(Ratio)-2*np.real(Ratio))
        return 10*np.log10(mason)   
    @property
    def UMason(self):
        S11=self.S11
        S22=self.S22
        S12=self.S12
        S21=self.S21
        Ratio=S21/S12
        Det=S11*S22-S12*S21
        kstb=1-np.abs(S11)**2-np.abs(S22)**2+np.abs(Det)**2
        kstb=kstb/(2*np.abs(S12*S21))
        mason = (np.abs(Ratio-1)**2)/(2*kstb*np.abs(Ratio)-2*np.real(Ratio))
        return mason   

    @property
    def stab_fact(self):
       # stransleted from ADS
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
       a1 = np.abs(self.S11)*np.abs(self.S11)+np.abs(self.S22)*np.abs(self.S22)
       a2 = np.abs(self.S11*self.S22-self.S12*self.S21)*np.abs(self.S11*self.S22-self.S12*self.S21)
       numer = 1.0-a1+a2
       denom = 2*np.abs(self.S12*self.S21)
       return numer/denom

    @property
    def stab_meas(self):
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
       a1 = np.abs(self.S11)*np.abs(self.S11) - np.abs(self.S22)*np.abs(self.S22)
       a2 = np.abs(self.S11*self.S22-self.S12*self.S21) * np.abs(self.S11*self.S22-self.S12*self.S21)
       return (1.0+a1-a2)
    
    @property
    def max_gain(self): 
       # this is translation of the function max_gain defined in ADS
       # in the ael file "C:\Program Files\Keysight\ADS2022_update1\expressions\ael\rf_system_fun.ael"
       # it is taken from rftools.py
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
        magS21=np.abs(self.S21)
        magS12=np.abs(self.S12)
        magSqS11=np.abs(self.S11)*np.abs(self.S11)
        magSqS22=np.abs(self.S22)*np.abs(self.S22)
        Bstab = self.stab_meas
        k = self.stab_fact
        k[k<=1]=1.0 # cancella i valori di k<1 e li sostituisce con 1
        gain0  = (magS21*magS21)/((1.0 - magSqS11)*(1.0 - magSqS22)) # 
        gain1  = (magS21/magS12)*(k - np.sqrt(k*k-1.0)) # MAG + 
        gain2  = (magS21/magS12)*(k + np.sqrt(k*k-1.0)) # MAG - 
        gainMSG =(magS21/magS12)
        maxGain = np.empty(np.shape(gain0), dtype=np.float_)
        for b in range(len(self.S11)): 
           for f in range(len(self.S11[b])):
              if magS12[b][f]==0: 
                 tmp = gain0[b][f]
              else: # magS12[f] > 0
                if Bstab[b][f] > 0 :
                    tmp = gain1[b][f]
                else:
                    tmp = gain2[b][f]
              if tmp < 1e-304:
                tmp=1e-304
              maxGain[b][f]=tmp
        return 10.0*np.log10(maxGain)

    def MSGdB(self):
    # taken from rft.S2PMSGdB(S):
        S11=self.S11
        S22=self.S22
        S12=self.S12
        S21=self.S21
        return 10*np.log10(S21/S12)
    def MSG(self):
    # taken from rft.S2PMSGdB(S):
        S11=self.S11
        S22=self.S22
        S12=self.S12
        S21=self.S21
        return S21/S12
    @property
    def Unilateral_FOM(self):
    # taken from rft.Unilateral_FOM(S):
        S11=self.S11
        S22=self.S22
        S12=self.S12
        S21=self.S21        
        S11a=np.abs(S11)
        S22a=np.abs(S22)
        S12a=np.abs(S12)
        S21a=np.abs(S21)
        U=S11a*S12a*S21a*S22a
        U=U/(1-S11a*S11a)
        U=U/(1-S22a*S22a)
        Uminus=10*np.log10(1/((1-U)*(1-U)))
        Uplus =10*np.log10(1/((1+U)*(1+U)))
        return (U, Uminus, Uplus)

    @property 
    def Gm(self):
        (fstart, fstop)=self.get_span('Gm')
        return np.array( [ (sum(tmp[fstart:fstop])/len(tmp[fstart:fstop])) for tmp in self.Gmv ] )
    @property
    def Gds(self):
        (fstart, fstop)=self.get_span('Gds')
        return np.array( [ (sum(tmp[fstart:fstop])/len(tmp[fstart:fstop])) for tmp in self.Gdsv ] )
 
    @property
    def Cgs(self):
        (fstart, fstop)=self.get_span('Cgs')
        return np.array( [ (sum(tmp[fstart:fstop])/len(tmp[fstart:fstop])) for tmp in self.Cgsv ] )
     
    @property
    def Cgd(self):
        (fstart, fstop)=self.get_span('Cgd')
        return np.array( [ (sum(tmp[fstart:fstop])/len(tmp[fstart:fstop])) for tmp in self.Cgdv ] )
 
    @property
    def Cds(self):
        (fstart, fstop)=self.get_span('Cds')
        return np.array( [ (sum(tmp[fstart:fstop])/len(tmp[fstart:fstop])) for tmp in self.Cdsv ] )
     
    @property
    def Ft(self):
        # Ft value is returned in GHz as usually expected
        (fstart, fstop)=self.get_span('Ft')
        return np.array( [ (sum(tmp[fstart:fstop])/1e9/len(tmp[fstart:fstop])) for tmp in self.Ftv ] )

    @property
    def Fmax(self):
        # Fmax value is returned in GHz as usually expected
        (fstart, fstop)=self.get_span('Fmax')
        return np.array( [ (sum(tmp[fstart:fstop])/1e9/len(tmp[fstart:fstop])) for tmp in self.Fmaxv ] )

class SS_STM(SS_trivial):
    """
    This class should reproduce exactly the extraction
    implemented in ICARE at ST
    """
    def __init__(self,Freq, Ym, Vg, Vd=None):
        SS_trivial.__init__(self, Freq, Ym, Vg, Vd) #call SS_trivial __init__ function, then adds to it :-)
        self.defs['A']=rft.findHz(self.F, 5e9, 5e9)  # pick Analogue gain at 5GHz by default
    @property
    def Cgdv(self):
        """
        perfect match, same as angelov
        """
        fact1=-np.imag(self.Y12)/(self.omega)
        fact2=np.real(self.Y12)/np.imag(self.Y12)
        return fact1*(1.0+fact2*fact2)    
    @property
    def Cggv(self): 
        """
        Cannot fully debug this one
        """
        return np.imag(self.Y11/self.omega)
    @property
    def Gdsv(self):
        """
        Perfect match. same as Angelov
        """
        return np.real(self.Y22)+np.real(self.Y12)
    @property
    def Coutv(self):
        """
        apparently, I get perfect fit, but no much sense, with Imag(Y22)/omega.
        Imag(Y22+Y12)/omega would be better
        """
        return np.imag(self.Y22)/self.omega
    @property
    def Rgdv(self): # what i cold find in mathematica, and it works (his is negative)
        """
        In this case STM agrees with me in correcting Angelov Formula 
        """
        return -np.real(1.0/self.Y12) 
    @property
    def Routv(self):
        """
        Cannot match this perfectly my best guess is Re[1/(Y22+Y21)]
        But 1/gds is a good match also: 1/Re[Y22] + 1/Re[Y21]
        but none is perfect :-(
        """
        return np.real(1.0/(self.Y22+self.Y12))
    @property
    def Rggv(self):
        """
        Called 'Enz' definition and valid at Vds=0V only! 
        """
        return np.real(self.Y11)/(np.imag(self.Y11)*np.imag(self.Y11))
    @property
    def AnGainv(self):
        """
        Analogue_gain gm/gds 
        """
        return self.Gmv/self.Gdsv
    @property
    def AnGain(self):
        (fstart, fstop)=self.get_span('AnGain')
        return np.array( [ (sum(tmp[fstart:fstop])/len(tmp[fstart:fstop])) for tmp in self.AnGainv ] )
    @property
    def Rout(self):
        (fstart, fstop)=self.get_span('Rout')
        return np.array( [ (sum(tmp[fstart:fstop])/len(tmp[fstart:fstop])) for tmp in self.Routv ] )
    @property
    def Rgg(self):
        (fstart, fstop)=self.get_span('Rgg')
        return np.array( [ (sum(tmp[fstart:fstop])/len(tmp[fstart:fstop])) for tmp in self.Rggv ] )
    @property
    def Cgg(self):
        (fstart, fstop)=self.get_span('Cgg')
        return np.array( [ (sum(tmp[fstart:fstop])/len(tmp[fstart:fstop])) for tmp in self.Cggv ] )
    @property
    def Cout(self):
        (fstart, fstop)=self.get_span('Cout')
        return np.array( [ (sum(tmp[fstart:fstop])/len(tmp[fstart:fstop])) for tmp in self.Coutv ] )
    @property
    def Rgd(self):
        (fstart, fstop)=self.get_span('Rgd')
        return np.array( [ (sum(tmp[fstart:fstop])/len(tmp[fstart:fstop])) for tmp in self.Rgdv ] )

class SS_Angelov(SS_trivial):
    @property
    def Gmv(self):
        arg=np.abs((self.Y21-self.Y12)/(self.Y11+self.Y12))
        return -arg/np.imag(1/(self.Y11+self.Y12))
    @property
    def Gdsv(self):
        return np.real(self.Y22)+np.real(self.Y12)
    @property
    def Cgsv(self):
        return -1.0/( np.imag( 1/(self.Y11+self.Y12) )*self.omega)
    @property
    def Cgdv(self):
        fact1=-np.imag(self.Y12)/(self.omega)
        fact2=np.real(self.Y12)/np.imag(self.Y12)
        return fact1*(1.0+fact2*fact2)
    @property
    def Cdsv(self):
        return np.imag(self.Y22+self.Y12)/self.omega
    @property
    def Riv(self):
        return np.real(1/(self.Y11+self.Y12))
    @property
    def Rgdv(self): # real definition from his book
        fact=-np.real(self.Y12)/np.imag(self.Y12)
        fact2=1+fact*fact
        return fact/fact2
    @property
    def Tauv(self):
        arg=(self.Y21-self.Y12)/(self.Y22+self.Y12)
        return -(np.deg2rad(np.angle(arg))+np.pi/2)/self.omega # cannot make it work!
  
    @property
    def Tau(self):
        (fstart, fstop)=self.get_span('Tau')
        return np.array( [ (sum(tmp[fstart:fstop])/len(tmp[fstart:fstop])) for tmp in self.Tauv ] )
    @property
    def Ri(self):
        (fstart, fstop)=self.get_span('Ri')
        return np.array( [ (sum(tmp[fstart:fstop])/len(tmp[fstart:fstop])) for tmp in self.Riv ] )
    @property
    def Rgd(self):
        (fstart, fstop)=self.get_span('Rgd')
        return np.array( [ (sum(tmp[fstart:fstop])/len(tmp[fstart:fstop])) for tmp in self.Rgdv ] )


#### Angelov definition of Rgd above is very difficult to understand and I could not derive
#### it in mathematica. The definition below is much simpler, but IT IS YIELDING
#### MUCH HIGHER VALUES!!!
#### ST also is using my definition of Rgdv

class SS_Angelov_Rgd(SS_Angelov):
    @property
    def Rgdv(self): # what i cold find in mathematica, and it works (his is negative)
        return -np.real(1.0/self.Y12)     



############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################

class PiMatrix(object):
    """ 
    it creates a Pi matrix from a set of S parameters (e.g. for an Open structure):
    
      p1 o-------||------o P2
            |   Yc    | 
            |         | 
        Ya ===       ===   Yb 
            |         |
            |         |
         o---------------o
          Where Y. (.=a, b or c) is either series
          Rs. + Cs.  i.e. Rsa + Csa  in case of Ya
    
               o---/\/\/\---||---o
                     Rs      Cs

          or parallel Rp. || Cp.

                        Cp
                   ----||----
              o____|         |___o
                   |         |
                   --/\/\/\---
                       Rp
    """
    def __init__(self, Freq, S=None, Y=None,Z=None):
        if S is not None:
            self.S=np.array(S)
            if Y is not None:
                print("WARNING, Init is given Sparams and Ymatrix. The Ymatrix is ignored!")
                Y=None # really forget the passed Yp matrix, since Sp is there
            self.Y=rft.StoY(self.S) # compute Ymatrix from Sp
            if Z is not None:
                print("WARNING, Init is given Sparams and Zmatrix. The Zmatrix is ignored!")
                Y=None # really forget the passed Yp matrix, since Sp is there
            self.Z=rft.YtoZ_2P(self.Y) # compute Ymatrix from Sp
        elif (Y is not None) and (Z is None):
            self.Y=np.array(Y)
            self.S=rft.YtoS(self.Y)
            self.Z=rft.YtoZ_2P(self.Y)
        elif (Y is None) and (Z is not None):
            self.Z=np.array(Z)
            self.Y=rft.ZtoY_2P(self.Z)
            self.S=rft.YtoS(self.Y)
        else:
            print("FATAL, cannon init the matrix!")
        self.SaveFig=False
        self.FigFmt='.svg'
        self.figsize=(5,4)
        d=np.shape(self.S)
        if len(Freq) != d[0]:
            print("Warning: leng of frequency sweep and length of matrix are not matching!")
        self.F=Freq
        self.omega=2*np.pi*self.F
        self.defs=dict()
        self.defs['C']=rft.findHz(self.F, 0.01e9, 150e9)
        self.defs['R']=rft.findHz(self.F, 0.01e9, 150e9)
        
        self.Y11=self.Y[:,0,0]
        self.Y22=self.Y[:,1,1]
        self.Y12=self.Y[:,0,1]
        self.Y21=self.Y[:,1,0]
        self.Yc=-0.5*(self.Y12+self.Y21)
        self.Ya=self.Y11-self.Yc
        self.Yb=self.Y22-self.Yc
        self.Za=1/self.Ya
        self.Zb=1/self.Yb
        self.Zc=1/self.Yc

    def set_span(self, var, start, stop):
        """  
        this methods set the dafault span to extract the 
        property defined in "var" as a string (upper case converted)
        from the start and stop values expressed in Hz
        """
        self.defs[var.upper()]=rft.findHz(self.F, start, stop)
    def get_span(self, var):
        testvar=var.upper()
        if testvar in self.defs:
            return self.defs[testvar]
        elif testvar[0] in self.defs:
            return self.defs[testvar[0]]
        else:
            print(f"Warning, could not identify {var} as a SS property")
            return (0, len(self.F)-1) # by default return all sweep    
    @property
    def Csav(self):
        return -1.0/(np.imag(self.Za)*(self.omega))
    @property    
    def Csbv(self):
        return -1.0/(np.imag(self.Zb)*(self.omega))
    @property
    def Cscv(self):
        return -1.0/(np.imag(self.Zc)*(self.omega))
    @property
    def Rsav(self):
        return np.real(self.Za)
    @property
    def Rsbv(self):
        return np.real(self.Zb)
    @property
    def Rscv(self):
        return np.real(self.Zc)
    @property
    def Cpav(self):
        return np.imag(self.Ya)/(self.omega)
    @property
    def Rpav(self):
        return 1.0/np.real(self.Ya)
    @property
    def Cpbv(self):
        return np.imag(self.Yb)/(self.omega)
    @property
    def Rpbv(self):
        return 1.0/np.real(self.Yb)
    @property
    def Cpcv(self):
        return np.imag(self.Yc)/(self.omega)
    @property
    def Rpcv(self):
        return 1.0/np.real(self.Yc)
    @property
    def Csa(self):
        (fstart, fstop)=self.get_span("Cs")
        return sum(self.Csav[fstart:fstop])/len(self.Csav[fstart:fstop])
    @property
    def Csb(self):
        (fstart, fstop)=self.get_span("Cs")
        return sum(self.Csbv[fstart:fstop])/len(self.Csbv[fstart:fstop])    
    @property
    def Csc(self):
        (fstart, fstop)=self.get_span("Cs")
        return sum(self.Cscv[fstart:fstop])/len(self.Cscv[fstart:fstop]) 
    @property
    def Cpa(self):
        (fstart, fstop)=self.get_span("Cp")
        return sum(self.Cpav[fstart:fstop])/len(self.Cpav[fstart:fstop])
    @property
    def Cpb(self):
        (fstart, fstop)=self.get_span("Cp")
        return sum(self.Cpbv[fstart:fstop])/len(self.Cpbv[fstart:fstop])    
    @property
    def Cpc(self):
        (fstart, fstop)=self.get_span("Cp")
        return sum(self.Cpcv[fstart:fstop])/len(self.Cpcv[fstart:fstop])
    @property
    def Rsa(self):
        (fstart, fstop)=self.get_span("Rs")
        return sum(self.Rsav[fstart:fstop])/len(self.Rsav[fstart:fstop])
    @property
    def Rsb(self):
        (fstart, fstop)=self.get_span("Rs")
        return sum(self.Rsbv[fstart:fstop])/len(self.Rsbv[fstart:fstop])    
    @property
    def Rsc(self):
        (fstart, fstop)=self.get_span("Rs")
        return sum(self.Rscv[fstart:fstop])/len(self.Rscv[fstart:fstop]) 
    @property
    def Rpa(self):
        (fstart, fstop)=self.get_span("Rp")
        return sum(self.Rpav[fstart:fstop])/len(self.Rpav[fstart:fstop])
    @property
    def Rpb(self):
        (fstart, fstop)=self.get_span("Rp")
        return sum(self.Rpbv[fstart:fstop])/len(self.Rpbv[fstart:fstop])    
    @property
    def Rpc(self):
        (fstart, fstop)=self.get_span("Rs")
        return sum(self.Rpcv[fstart:fstop])/len(self.Rpcv[fstart:fstop])     
    def plotC(self, what=['Csav', 'Csbv', 'Cscv']):
        if not what:
            what=['Csav', 'Csbv', 'Cscv']
        if type(what) is str:
            what=[what]
        fig=plt.figure(figsize=self.figsize)
        Cp = fig.add_subplot(111) # one row, two columns, second plot
        for w in what:          
            Cp.plot(self.F*1e-9,1e15*getattr(self, w), label=w)
        Cp.legend(loc='best'); Cp.set_xlabel('f [GHz]'); Cp.set_ylabel("C [fF]")  ; # Cp.set_ylim(0,10);
        Cp.minorticks_on()
        Cp.grid(which='major',ls='-')
        Cp.grid(which='minor',ls=':')
        return fig, Cp
    def plotR(self, what=['Rsav', 'Rsbv', 'Rscv']):
        if not what:
            what=['Rsav', 'Rsbv', 'Rscv']
        if type(what) is str:
            what=[what]
        fig=plt.figure(figsize=self.figsize)
        Cp = fig.add_subplot(111) # one row, two columns, second plot
        for w in what:          
            Cp.plot(self.F*1e-9,getattr(self, w), label=w)
        Cp.legend(loc='best'); Cp.set_xlabel('f [GHz]'); Cp.set_ylabel("R [Ohms]")  ; # Cp.set_ylim(0,10);
        Cp.minorticks_on()
        Cp.grid(which='major',ls='-')
        Cp.grid(which='minor',ls=':')
        return fig, Cp


class TeeMatrix(object):
    """ 
    it creates a Tee matrix from a set of S parameters (e.g. for an Short structure):
              __      __
      P1 o---|__|----|__|---o P2
              Za  |   Zb
                  _
               Zc| |
                 |_|
                  |
         o----------------o

          Where Z. (.=a, b or c) is either series
          Rs. + Ls.  i.e. Rsa + Csa  in case of Ya
    
               o---/\/\/\---UUUU---o
                     Rs      Ls

          or parallel Rp. || Lp.

                       Lp
                   ---UUUUU---
              o____|         |___o
                   |         |
                   --/\/\/\---
                       Rp
    """
    def __init__(self, Freq, S=None, Y=None, Z=None):
        if S is not None:
            self.S=np.array(S)
            if Y is not None:
                print("WARNING, Init is given Sparams and Ymatrix. The Ymatrix is ignored!")
                Y=None # really forget the passed Yp matrix, since Sp is there
            self.Y=rft.StoY(self.S) # compute Ymatrix from Sp
            if Z is not None:
                print("WARNING, Init is given Sparams and Zmatrix. The Zmatrix is ignored!")
                Y=None # really forget the passed Yp matrix, since Sp is there
            self.Z=rft.YtoZ_2P(self.Y) # compute Ymatrix from Sp
        elif (Y is not None) and (Z is None):
            self.Y=np.array(Y)
            self.S=rft.YtoS(self.Y)
            self.Z=rft.YtoZ_2P(self.Y)
        elif (Y is None) and (Z is not None):
            self.Z=np.array(Z)
            self.S=rft.YtoS(rft.ZtoY_2P(self.Z))
            self.Y=rft.ZtoY_2P(self.Z)
        else:
            print("FATAL, cannon init the matrix!")
        self.SaveFig=False
        self.FigFmt='.svg'
        self.figsize=(5,4)
        d=np.shape(self.S)
        if len(Freq) != d[0]:
            print("Warning: leng of frequency sweep and length of matrix are not matching!")
        self.F=Freq
        self.omega=2*np.pi*self.F
        self.defs=dict()
        self.defs['L']=rft.findHz(self.F, 0.01e9, 150e9)
        self.defs['R']=rft.findHz(self.F, 0.01e9, 150e9)
        
        self.Z11=self.Z[:,0,0]
        self.Z22=self.Z[:,1,1]
        self.Z12=self.Z[:,0,1]
        self.Z21=self.Z[:,1,0]
        self.Zc=0.5*(self.Z12+self.Z21)
        self.Za=self.Z11-self.Zc
        self.Zb=self.Z22-self.Zc
        self.Ya=1/self.Za
        self.Yb=1/self.Zb
        self.Yc=1/self.Zc

    def set_span(self, var, start, stop):
        """  
        this methods set the dafault span to extract the 
        property defined in "var" as a string (upper case converted)
        from the start and stop values expressed in Hz
        """
        self.defs[var.upper()]=rft.findHz(self.F, start, stop)
    def get_span(self, var):
        testvar=var.upper()
        if testvar in self.defs:
            return self.defs[testvar]
        elif testvar[0] in self.defs:
            return self.defs[testvar[0]]
        else:
            print(f"Warning, could not identify {var} as a SS property")
            return (0, len(self.F)-1) # by default return all sweep    
    @property
    def Lsav(self):
        return (np.imag(self.Za)/(self.omega))
    @property    
    def Lsbv(self):
        return (np.imag(self.Zb)/(self.omega))
    @property
    def Lscv(self):
        return (np.imag(self.Zc)/(self.omega))
    @property
    def Rsav(self):
        return np.real(self.Za)
    @property
    def Rsbv(self):
        return np.real(self.Zb)
    @property
    def Rscv(self):
        return np.real(self.Zc)
    @property
    def Lpav(self):
        return -1.0/(np.imag(self.Ya)*(self.omega))
    @property
    def Rpav(self):
        return 1.0/np.real(self.Ya)
    @property
    def Lpbv(self):
        return -1.0/(np.imag(self.Yb)*(self.omega))
    @property
    def Rpbv(self):
        return 1.0/np.real(self.Yb)
    @property
    def Lpcv(self):
        return -1.0/(np.imag(self.Yc)*(self.omega))
    @property
    def Rpcv(self):
        return 1.0/np.real(self.Yc)
    @property
    def Lsa(self):
        (fstart, fstop)=self.get_span("Ls")
        return sum(self.Lsav[fstart:fstop])/len(self.Lsav[fstart:fstop])
    @property
    def Lsb(self):
        (fstart, fstop)=self.get_span("Ls")
        return sum(self.Lsbv[fstart:fstop])/len(self.Lsbv[fstart:fstop])    
    @property
    def Lsc(self):
        (fstart, fstop)=self.get_span("Ls")
        return sum(self.Lscv[fstart:fstop])/len(self.Lscv[fstart:fstop]) 
    @property
    def Lpa(self):
        (fstart, fstop)=self.get_span("Lp")
        return sum(self.Lpav[fstart:fstop])/len(self.Lpav[fstart:fstop])
    @property
    def Lpb(self):
        (fstart, fstop)=self.get_span("Lp")
        return sum(self.Lpbv[fstart:fstop])/len(self.Lpbv[fstart:fstop])    
    @property
    def Lpc(self):
        (fstart, fstop)=self.get_span("Lp")
        return sum(self.Lpcv[fstart:fstop])/len(self.Lpcv[fstart:fstop])
    @property
    def Rsa(self):
        (fstart, fstop)=self.get_span("Rs")
        return sum(self.Rsav[fstart:fstop])/len(self.Rsav[fstart:fstop])
    @property
    def Rsb(self):
        (fstart, fstop)=self.get_span("Rs")
        return sum(self.Rsbv[fstart:fstop])/len(self.Rsbv[fstart:fstop])    
    @property
    def Rsc(self):
        (fstart, fstop)=self.get_span("Rs")
        return sum(self.Rscv[fstart:fstop])/len(self.Rscv[fstart:fstop]) 
    @property
    def Rpa(self):
        (fstart, fstop)=self.get_span("Rp")
        return sum(self.Rpav[fstart:fstop])/len(self.Rpav[fstart:fstop])
    @property
    def Rpb(self):
        (fstart, fstop)=self.get_span("Rp")
        return sum(self.Rpbv[fstart:fstop])/len(self.Rpbv[fstart:fstop])    
    @property
    def Rpc(self):
        (fstart, fstop)=self.get_span("Rs")
        return sum(self.Rpcv[fstart:fstop])/len(self.Rpcv[fstart:fstop])     
    def plotL(self, what=['Lsav', 'Lsbv', 'Lscv']):
        if type(what) is str:
            what=[what]
        fig=plt.figure(figsize=self.figsize)
        Lp = fig.add_subplot(111) # one row, two columns, second plot
        for w in what:          
            Lp.plot(self.F*1e-9,1e12*getattr(self, w), label=w)
        Lp.legend(loc='best'); Lp.set_xlabel('f [GHz]'); Lp.set_ylabel("L [pH]")  ; # Cp.set_ylim(0,10);
        Lp.minorticks_on()
        Lp.grid(which='major',ls='-')
        Lp.grid(which='minor',ls=':')
        return fig, Lp
    def plotR(self, what=['Rsav', 'Rsbv', 'Rscv']):
        if type(what) is str:
            what=[what]
        fig=plt.figure(figsize=self.figsize)
        Cp = fig.add_subplot(111) # one row, two columns, second plot
        for w in what:          
            Cp.plot(self.F*1e-9,getattr(self, w), label=w)
        Cp.legend(loc='best'); Cp.set_xlabel('f [GHz]'); Cp.set_ylabel("R [Ohms]")  ; # Cp.set_ylim(0,10);
        Cp.minorticks_on()
        Cp.grid(which='major',ls='-')
        Cp.grid(which='minor',ls=':')   
        return fig, Cp





#!/usr/bin/env python
# coding: utf-8

# # icare.py - IO with Crolles 
** DO NOT EDIT THIS .py FILE. EDIT THE JUPYTER NOTEBOOK **
# To actually create the module, export as an executable python script

# Change log:
# - v 0.1: inclusion of first macros

# In[ ]:


# this line is useless, will be added by export filter anyway. here for redundant checks
# -*- coding: utf-8 -*-

# import re   as re 
import numpy as np  # so far is the only dependable module to be included
# import pickle

################################################################################


# functions preparse, ParseIdVg, ParseValues, parseMatrix were developted to read in huge number of files sent by Joao in "ICARE" native format. 
# Format is quite similar, but different from what modelTK is reading in.

# ## Native icare read in (header with Hash)

# In[ ]:


def preparse(filename, readAllKeys=False): 
        """
         this routines opens the file given as input parameter then:
          - check that the file actually exists 
          - chack that files cames from icare (header has # as first char)
          - reads in the header variables and stores them in a hash
          - reads the column headers and store the column names in a hash
         file format accepted is quite strict:
           # key1 val1
           # key2 val2
           header1 header2 headern
           val1 val2  ...  valn
          no comments and no empty lines permitted!
          return header, column
                           list of strings
                   hash
          readAllKeys  =  reads all header, even if key has no associated value
                          by defaults ignores header keys that have no values associated
        """
        header=dict()
        columns=list()
        hline=list()
        with open(filename) as f:
            for line in f:
                if line[0] != '#': break
                hline = line.strip().split(None ,2) # 2 is the max numb of splits
                if len(hline)==2 :
                    if readAllKeys:  header[hline[1]]=None
                elif len(hline)==3 : 
                    header[hline[1]]=hline[2]
                else:
                    print(f"Warning: non conforming line {line}")
            columns=line.strip().split()
        return header, columns


# In[ ]:


def parseIdVg(filename, colVg, colVd, colId, colIg):
    """ 
    this routine parse a IdVg sweep in all lines where it makes sense, that is in all lines 
    where column[colVg], column[colVd], column[colId] and column[colId] are not "None"
    ex. 
    pd['IdVg']=icare.parseIdVg(datadir+"MOS_3.txt", colVg=5, colVd=4, colIg=9, colId=8)
    returns a dictionary
    e.g.
    pd['IdVg']['vd']=[0, 0.05, 0.8]
    pd['IdVg']['vg']= linspace 0.0 to 0.8 
    pd['IdVg']['id']= (len(vd), len(vg)) real value matrix
    pd['IdVg']['ig']= (len(vd), len(vg)) real value matrix
    """
    IdVg=dict()
    vVd=list()
    vVg=list()
    vId=list()
    vIg=list()
    data=list()
    first=True # skip one more line with headers
    with open(filename) as f:
        for line in f: 
            if line[0]=='#': continue # skips header
            if first:
                # parsing header line
                data = line.strip().split()
                print(f"Parsing {data[colVg]} as Vg, {data[colVd]} as Vd, {data[colIg]} as Ig, {data[colId]} as Id")
                first=False
                continue

            data=line.strip().split()
            Vd=data[colVd]
            Vg=data[colVg]
            Id=data[colId]
            Ig=data[colIg]
            if Vg == "None" : continue
            if Vd == "None" : continue
            if Id == "None" : continue
            if Ig == "None" : continue
            vVd.append(float(Vd))
            vId.append(float(Id))
            vVg.append(float(Vg))
            vIg.append(float(Ig))
            # file should have been parsed completely now
        IdVg['vvd']=np.array(vVd)
        IdVg['vvg']=np.array(vVg)
        IdVg['vid']=np.array(vId)
        IdVg['vig']=np.array(vIg)

        allVg=list(sorted(set(vVg)))
        allVd=list(sorted(set(vVd)))
        IdVg['vd']=np.array(allVd)
        IdVg['vg']=np.array(allVg)
        IdVg['id']=np.empty((len(allVd), len(allVg)), dtype=float)
        IdVg['ig']=np.empty((len(allVd), len(allVg)), dtype=float)

        # magic trick, suppose there are no redundant Vd, Vg data
        # and all lines are orderer in Vg ascending:
        for i, valVd in enumerate(allVd):   # for each Vd values
            vdmask=(IdVg['vvd']==valVd)  # mask for e.g. Vd=1.1 Values
            IdVg['id'][i]=IdVg['vid'][vdmask]
            IdVg['ig'][i]=IdVg['vig'][vdmask]
    return IdVg


# In[ ]:


def parseValues(filename, colvg, colvd, colf, tuplelist):
    """ 
    this routine parse a complex, int or float vector sweep
    ex.
    dat = icare.parseValues(datadir+"MOS_3.txt",[("vd", 4), ( "vg", 5), ("freq", 6)])
    returns a dictionary with vd, vg, freq array, when ALL the values specified
    are not None
    This is important: e.g. vd and vg may be specified for RF and DC measurements
    on the same file. so they are not None in basically all lines, while freq is None
    in all lines refering to DC measurements,
    """
    var=list()
    numpyvar=dict()
    varpos=list()
    vartype=list()
    varname=list()
    
    ## paranoid
    vvg  =list()
    vvd  =list()
    vfreq=list()

    for tvar in tuplelist:
        varname.append(tvar[0])
        varpos.append(int(tvar[1]))
        if len(tvar)==2:
            vartype.append(float)
        elif len(tvar)==3:
            vartype.append(tvar[2])
        var.append(list())
     
    ##  varname=[  "NF"  ,   "Gopt"]   # aka dict name
    ##  vattype=[ float,    complex ]  # aka datatype
    ##  varpos =[  11    ,   12  ]     # column  to read
    ##  var    =[ [ ] , [] , ]       # list of list ready to be filled

    first=True # skip one more line with headers
    with open(filename) as f:
        for line in f: 
            if line[0]=='#': continue # skips header
            if first:
                # parsing header line
                data = line.strip().split()
                print(f"Parsing", end="")
                for i,name in enumerate(varname):
                    print(f" {data[varpos[i]]} as {name},", end="")
                print(f" also {data[colvd]} as vd, {data[colvg]} as vg and {data[colf]} as freq")
                first=False
                continue

            data=line.strip().split()
            # test if none of the variables are "None"
            relevant=all([data[i] != "None" for i in varpos])
            if relevant:
                for i,name in enumerate(varname):
                    thisdata=data[varpos[i]]
                    var[i].append(thisdata)
                vvd.append(data[colvd])
                vvg.append(data[colvg])
                vfreq.append(int(float(data[colf])))
    allvg=list(sorted(set(vvg))) # ordered list of vg values
    allvd=list(sorted(set(vvd))) # ordered list of vd values
    allf =list(sorted(set(vfreq))) # ordered list of freq values
    fsweeplen=len(allf)
    allbias=len(allvg)*len(allvd)

    if len(vvd) != allbias*fsweeplen:
        print(allvg)
        print(allvd)
        print(allf)
        print(allbias)
        raise("inconsistent bias point data")

    for i,name in enumerate(varname):
        numpyvar[name]=np.empty((allbias,fsweeplen), dtype=vartype[i])
    numpyvar['vg']=np.empty(allbias, dtype=float)
    numpyvar['vd']=np.empty(allbias, dtype=float)
    numpyvar['F'] =np.array(allf)

    #numpyvar['vg']=np.array(allvg)
    #numpyvar['vd']=np.array(allvd)

    #numpyvar['vvg']=np.empty(allbias, dtype=float)
    #numpyvar['vvd']=np.empty(allbias, dtype=float)

    #numpyvar['vg']=np.array(allvg)
    #numpyvar['vd']=np.array(allvd)


    for i in range(allbias): # for each bias point
        # paranoid, this is a bias point, check all freq values are correct and all vd, vg 
        # values are constant
        start=i*fsweeplen
        vd0=vvd[start] # vd and vg at the first f point
        vg0=vvg[start]
        for j in range(fsweeplen):
            if vvd[start+j] != vd0: raise("inconsistent vd bias value")
            if vvg[start+j] != vg0: raise("inconsistent vd bias value")
            if vfreq[start+j] != allf[j]: raise("inconsistent frequency value")            
        numpyvar['vg'][i]=vg0
        numpyvar['vd'][i]=vd0
        for k,name in enumerate(varname): # for each reading variable
            for f in range(fsweeplen): # for each freq point
                numpyvar[name][i][f]=var[k][start+f]
    return numpyvar


# In[ ]:


def parseMatrix(filename, cfreq, M11, M12, M21, M22):
    """ 
    this routine parse a complex values 2x2 matrix and associated freq sweep 
    """
    vF=list()
    vS11=list()
    vS21=list()
    vS12=list()
    vS22=list()
 

    first=True # skip one more line with headers
    with open(filename) as f:
        for line in f: 
            if line[0]=='#': continue # skips header
            if first:
                # parsing header line
                data = line.strip().split()
                print(f"Parsing {data[cfreq]} as f, {data[M11]} as 11, {data[M12]} as 12, "
                      f"{data[M21]} as 21, {data[M22]} as 22")
                first=False
                continue

            data=line.strip().split()
            freq=data[cfreq]
            if freq == "None" : continue
            vF.append(int(float(freq)))
            vS11.append(complex(data[M11]))
            vS12.append(complex(data[M12]))
            vS21.append(complex(data[M21]))
            vS22.append(complex(data[M22]))

            # file should have been parsed completely now

    allF=list(sorted(set(vF)))
    freqpoints=len(allF)
    biaspoints=len(vS22)//freqpoints  # integer division
    if biaspoints*freqpoints != len(vS22):
        raise "consistency problem"
    if biaspoints>1:
        Matrix=np.empty((biaspoints, freqpoints, 2, 2), dtype=complex)
    else:
        Matrix=np.empty((freqpoints, 2, 2), dtype=complex)
    sweep=np.array(allF)

    if biaspoints>1:
        for i in range(biaspoints):
            Matrix[i,:,0,0]=vS11[i*freqpoints:(i+1)*freqpoints]
            Matrix[i,:,1,0]=vS21[i*freqpoints:(i+1)*freqpoints]
            Matrix[i,:,0,1]=vS12[i*freqpoints:(i+1)*freqpoints]
            Matrix[i,:,1,1]=vS22[i*freqpoints:(i+1)*freqpoints]
            # paranoid
            for j in range(len(allF)):
                if sweep[j]!=vF[i*freqpoints+j]:
                    raise("mimatch in f values!")
    else:
        Matrix[:,0,0]=vS11
        Matrix[:,1,0]=vS21
        Matrix[:,0,1]=vS12
        Matrix[:,1,1]=vS22

    return sweep, Matrix


#!/usr/bin/env python
# coding: utf-8

# # ICCAP.py - IO with mom files

# **DO NOT EDIT THIS .py FILE. EDIT THE JUPYTER NOTEBOOK**

# To actually create the module, export as an executable python script

# Change log:
# - v 0.2: fix and adds DC import
# - v 0.1: inclusion of first macros

# In[4]:


# this line is useless, will be added by export filter anyway. here for redundant checks
# -*- coding: utf-8 -*-

import numpy as np  # so far is the only dependable module to be included
import os # to check if file exists
import re 
################################################################################


# functions: parse
# 

# ## Parse  of sparam data from  ICCAP native mdm files

# In[5]:


def parseRF(filename, readAllKeys=False, trans={'Grille':'vvg','Vds':'vvd','vg':'vvg','vd':'vvd'}, fenc='Latin1', verbose=False): 
        """
         this routines opens the file given as input parameter then:
          - check that the file actually exists 
          - check that files cames with proper ICCAP headers
          - reads in the header variables and stores them in a hash
          - reads the column headers and store the column names in a hash
          
          by default file is encoded as 'Latin1' i.e. windows origin
          in this routine we tacitly aussme the primary sweep is frequency
          i.e. we are reading in some Sparams. Of a passive (simple sweep, no bias)
          or of an active (with Vg and Vd bias)
        """
        header=dict() # headders stuff
        DBs=list()    # Databases, i.e. blocks of data
        lines=list()  # input file
        with open(filename, encoding=fenc) as f:
            lines=f.readlines()
        lines=list(map(str.strip, lines))
        # check the first line for the [ VERSION=6.0 ]
        if not re.search("^! VERSION = (.*)", lines[0]):
            print(f"FATAL:missing first line check! {lines[0]}")
            return
        # finde header start and stop lines
        start=-1
        stop=-1
        for i,line in enumerate(lines):
            if re.search("^BEGIN_HEADER",line): start=i
            if re.search("^END_HEADER", line): stop=i
        if not start>0 and stop>0 : 
            print(f"FATAL: Could not find start and end of header! start {start} stop {stop} ")
        else:
            header=lines[start+1:stop]
        # print(header)
        # in case of ST files, pls add section here
        for i, line in enumerate(header):
            if re.search("^.*ICCAP_INPUT",line): inputvars=i
            if re.search("^.*ICCAP_OUTPUT", line): outputvars=i
        inlist=header[inputvars+1:outputvars] # line containing an input var directive, i.e. a sweep variable
        outlist=header[outputvars+1:]         # line containing an output, i.e. a measured or a simulated values
        # print(inlist)
        inputs=[line.split() for line in inlist]
        outputs=[line.split() for line in outlist]
        # print(outputs)
        varBias=list()
        ### end of common part. start to make decision about the sweep we are reading in
        # look for primary sweep in frequency
        for i, var in enumerate(inputs):
            if var[1]=="F":
                fsname=var[0]
                if not int(var[3])==1: 
                    print(f"WARNING:frequency is not the primary sweep but {var[3]}!?")
                fstype=var[2]
                if fstype=='LIST': fslen=int(var[4])
                if fstype=='LIN': fslen=int(var[6])
                if fstype=='LOG': fslen=int(var[6])                

            if var[1]=='V':
                varBias.append(var[0])
            if var[1]=='I':
                varBias.append(var[0])   
        if not fsname=='F' and verbose:
            print(f"WARNING: frequency sweep variable name is '{fsname}', will rename it to 'F'")
        if len(varBias)>0 and verbose:
            print(f"File contains bias sweeps on {varBias}")
  
        tosave=list()
        for i, var in enumerate(outputs):
            if var[1]=="S":
                spname=var[0]
            if var[1]=='V':
                tosave.append(var[0])
            if var[1]=='I':
                tosave.append(var[0])   
        if not spname=='Sp' and verbose:
            print(f"WARNING: Sparam var name is '{spname}' will rename it to 'Sp'")
        if len(tosave)>0 and verbose:
            print(f"found extra output variables named {tosave}")
        
        
        dbstart=list()
        dbstop=list()
        for i,line in enumerate(lines):
            if re.search("^BEGIN_DB",line): dbstart.append(i+1)
            if re.search("^END_DB", line): dbstop.append(i)
            
        DBsraw=list() # list of databases, there's at least one if no bias
        for i, start in enumerate(dbstart):
            DBsraw.append(lines[start:dbstop[i]])
        colstring=DBsraw[0][len(inputs)]
        cols=colstring.split()
        save=dict()

        for vname in varBias:
            save[vname]=list()

        for vname in tosave:
            save[vname]=dict()
            save[vname]['col']=cols.index(vname)
            save[vname]['data']=list()
        s11r=cols.index(f"R:{spname}(1,1)")
        s11i=cols.index(f"I:{spname}(1,1)")
        s22r=cols.index(f"R:{spname}(2,2)")
        s22i=cols.index(f"I:{spname}(2,2)")
        s21r=cols.index(f"R:{spname}(2,1)")
        s21i=cols.index(f"I:{spname}(2,1)") 
        s12r=cols.index(f"R:{spname}(1,2)")
        s12i=cols.index(f"I:{spname}(1,2)")         
        save['F']=list()
        
        #print(cols)
        #print(f"save   is {save} ")
        #print(f"tosave is {tosave}")
        DBs=list()
        
        # parsebias points
        for DB in DBsraw:
            if not DB[len(inputs)]==colstring:
                print("FATAL: sweep header are not matching!")
                return
            DBhead=DB[0:len(inputs)]
            #print(DBhead)
            for line in DBhead:
                b=re.search('^\s*ICCAP_VAR (\S+)\s+(\S+)',line)
                if b:
                    # gropu 0 is the string
                    # group (1) is the first parentheses
                    # group (2) is the second parentheses
                    save[b.group(1)].append(float(b.group(2))+0.0)
                    #print(f"{b.group(1)} is {b.group(2)}")
        
            DBdata=DB[len(inputs)+1:] # acconts fot an empty lne and a column header
            DBsplit=[line.split() for line in DBdata]
            DBs.append(DBsplit)
        # guess number of bias points:
        if len(varBias)>0:
            bp=len(save[varBias[0]])
            save['Sp']=np.zeros((bp,fslen, 2,2), dtype=np.complex_)
        else:
            bp=0
            save['Sp']=np.zeros((fslen, 2,2), dtype=np.complex_)    
        DB=DBs[0]
        save['F']=np.array([float(p[0]) for p in DB])
        for i, DB in enumerate(DBs):
            if not len(DB)==fslen:  
                print("FATAL: inconsistent DB length")
                print(DBdata)
                return
            for j, DBl in enumerate(DB):                 
                if bp>0: 
                    save['Sp'][i,j,0,0]=complex( float(DBl[s11r]), float(DBl[s11i]))
                    save['Sp'][i,j,1,1]=complex( float(DBl[s22r]), float(DBl[s22i]))
                    save['Sp'][i,j,1,0]=complex( float(DBl[s21r]), float(DBl[s21i]))
                    save['Sp'][i,j,0,1]=complex( float(DBl[s12r]), float(DBl[s12i]))
                else:
                    save['Sp'][j,0,0]=complex( float(DBl[s11r]), float(DBl[s11i]))
                    save['Sp'][j,1,1]=complex( float(DBl[s22r]), float(DBl[s22i]))
                    save['Sp'][j,1,0]=complex( float(DBl[s21r]), float(DBl[s21i]))
                    save['Sp'][j,0,1]=complex( float(DBl[s12r]), float(DBl[s12i]))
            # trick only last value is saved!!
            for var in tosave:
                save[var]['data'].append(float(DBl[save[var]['col']]))
        for var in tosave:
            save[var]=np.array(save[var]['data'])
        for var in varBias:
            if var in trans:
                save[trans[var]]=np.array(save[var])
                del save[var]
            else:            
                save[var]=np.array(save[var])
        return save            
            
     


# In[6]:


def parseDC(filename, trans={"R:gd(1,1)":'gds', 'R:gm(1,1)':'gm'}, fenc='Latin1', verbose=False): 
    """
         this routines opens the file given as input parameter then:
          - check that the file actually exists 
          - check that files cames with proper ICCAP headers
          - reads in the header variables and stores them in a hash
          - reads the column headers and store the column names in a hash
          
          by default file is encoded as 'Latin1' i.e. windows origin
          in this routine we tacitly aussme the primary sweep is a voltage, i.e. 
    """
    #filename="C:/home/LABLUCCI/230707_STC65_ALTAIR_IEMN/Q201122_wf09_120723_yannick_mos_thomas_st/Q1BHF17_1_ret32/IDEV.mdm"
    #filename="C:/home/LABLUCCI/230707_STC65_ALTAIR_IEMN/Q201122_wf09_120723_yannick_mos_thomas_st/Q1BHF17_1_ret32/vth.mdm"
    #verbose=True
    #fenc='Latin1'

    header=dict() # headders stuff
    DBs=list()    # Databases, i.e. blocks of data
    lines=list()  # input file
    with open(filename, encoding=fenc) as f:
        lines=f.readlines()
    lines=list(map(str.strip, lines))
    # check the first line for the [ VERSION=6.0 ]
    if not re.search("^! VERSION = (.*)", lines[0]):
        print(f"FATAL:missing first line check! {lines[0]}")
        return
    # finde header start and stop lines
    start=-1
    stop=-1
    for i,line in enumerate(lines):
        if re.search("^BEGIN_HEADER",line): start=i
        if re.search("^END_HEADER", line): stop=i
    if not start>0 and stop>0 : 
        print(f"FATAL: Could not find start and end of header! start {start} stop {stop} ")
        return
    else:
        header=lines[start+1:stop]
    # print(header)
    # in case of ST files, pls add section here
    for i, line in enumerate(header):
        if re.search("^.*ICCAP_INPUT",line): inputvars=i
        if re.search("^.*ICCAP_OUTPUT", line): outputvars=i
    inlist=header[inputvars+1:outputvars] # line containing an input var directive, i.e. a sweep variable
    outlist=header[outputvars+1:]         # line containing an output, i.e. a measured or a simulated values
    # print(inlist)
    inputs=[line.split() for line in inlist]
    outputs=[line.split() for line in outlist]
    # print(outputs)
    varBias=list()
    ### end of common part. start to make decision about the sweep we are reading in

    # look for primary sweep and all bias sweep
    # will fill a varBias list
    primary=None
    for i, var in enumerate(inputs):
        if var[1]=='V':
            varBias.append(var[0])
        elif var[1]=='I':
            varBias.append(var[0])         
        else:
            print(f"FATAL: vartype {var[1]} for {var[0]} is a type not supported") 
            return
        if var[6] in ("LIN", "LOG", "LIST"):
            if int(var[7])==1:
                primary=var[0]
    if not primary:
        print(f"FATAL: could not identify primary sweep")
        return
    if len(varBias)>0 and verbose:
        print(f"File contains bias sweeps on {varBias}")
    if verbose: print(f"primary sweep is {primary}")

    ## start with outputs
    # will fill a tosave dictionary
    tosave=list()
    for i, var in enumerate(outputs):
        if var[1]=='V':
            tosave.append(var[0])
        elif var[1]=='I':
            tosave.append(var[0])
        elif var[1]=='U':
            tosave.append(f"R:{var[0]}(1,1)")         
    if len(tosave)>0 and verbose:
        print(f"found output variables named {tosave}")

    # find all available databases
    dbstart=list()
    dbstop=list()
    for i,line in enumerate(lines):
        if re.search("^BEGIN_DB",line): dbstart.append(i+1)
        if re.search("^END_DB", line): dbstop.append(i)

    DBsraw=list() # list of databases, there's at least one if no bias
    for i, start in enumerate(dbstart):
        DBsraw.append(lines[start:dbstop[i]])
    colstring=DBsraw[0][len(inputs)]
    cols=colstring.split()
    save=dict()

    storeBiasInfo=dict()
    for vname in varBias:
        save[vname]=list()
        save[f'v{vname}']=list()
        if vname==primary: continue
        storeBiasInfo[vname]=list()

    for vname in tosave:
        save[vname]=dict()
        save[vname]['col']=cols.index(vname)
        save[vname]['data']=list() 

    if verbose:
        print(f"cols are {cols}")
        print(f"save   is {save} ")
        print(f"tosave is {tosave}")
    DBs=list()

    # parsebias points

    for DB in DBsraw:
        if not DB[len(inputs)]==colstring: # test for same header for DB columns
            print("FATAL: sweep header are not matching!")
            return
        DBhead=DB[0:len(inputs)]
        #print(DBhead)
        # separates the header, with secondary bias points
        # from the bulk of the database with the outputs
        for line in DBhead:
            b=re.search('^\s*ICCAP_VAR (\S+)\s+(\S+)',line)
            if b:
                # gropu 0 is the string
                # group (1) is the first parentheses
                # group (2) is the second parentheses
                storeBiasInfo[b.group(1)].append(float(b.group(2))+0.0)
                #print(f"{b.group(1)} is {b.group(2)}")

        DBdata=DB[len(inputs)+1:] # acconts fot an empty lne and a column header
        DBsplit=[line.split() for line in DBdata]
        DBs.append(DBsplit)
    # print(f"store is {storeBiasInfo}")
    save[primary]=np.array([float(p[0]) for p in DBs[0]])
    for var in storeBiasInfo:
        save[var]=np.array([float(v) for v in storeBiasInfo[var]])
    for i, DB in enumerate(DBs):
        save[f'v{primary}'].append([float(p[0]) for p in DB])
        for var in storeBiasInfo.keys():
            save[f'v{var}'].append([storeBiasInfo[var][i] for dummy in DB]) 
        for var in tosave:
            save[var]['data'].append([float(p[save[var]['col']]) for p in DB])

    for var in tosave:
        save[var]=np.array(save[var]['data'])
    #for var in set(varBias) | set(tosave):
    #    if var in trans:
    #        save[trans[var]]=np.array(save[var])
    #        del save[var]
    #    else:            
    #        save[var]=np.array(save[var])
    return save


# test the sparamenter imput

# In[7]:


#d1=parseRF("C:/home/LABLUCCI/230707_STC65_ALTAIR_IEMN/Q201122_wf09_120723_yannick_mos_thomas_st/Q1BHF17_1_ret32/Sparameters.mdm", verbose=True)
#d2=parseRF("C:/home/LABLUCCI/230707_STC65_ALTAIR_IEMN/Q201122_wf09_120723_yannick_mos_thomas_st/Q1BHF17_1_ret32/Self.mdm", verbose=True)
#d3=parseRF("C:/home/LABLUCCI/230707_STC65_ALTAIR_IEMN/Q201122_wf09_120723_yannick_mos_thomas_st/Q1BHF17_1_ret32/S_ft.mdm", verbose=True)

#d4=parseDC("C:/home/LABLUCCI/230707_STC65_ALTAIR_IEMN/Q201122_wf09_120723_yannick_mos_thomas_st/Q1BHF17_1_ret32/IDEV.mdm",verbose=True)
#d5=parseDC("C:/home/LABLUCCI/230707_STC65_ALTAIR_IEMN/Q201122_wf09_120723_yannick_mos_thomas_st/Q1BHF17_1_ret32/vth.mdm",verbose=True)


# In[8]:


#d4.keys()


# In[9]:


# tg=d4
# for k, dat in tg.items():
#     print(f"{k} is {np.shape(dat)}")


# In[187]:


# d1.keys()
# for k, dat in d1.items():
#     print(f"{k} is {np.shape(dat)}")


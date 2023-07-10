#!/usr/bin/env python
# coding: utf-8

# # ICCAP.py - IO with mom files

# **DO NOT EDIT THIS .py FILE. EDIT THE JUPYTER NOTEBOOK**

# To actually create the module, export as an executable python script

# Change log:
# - v 0.1: inclusion of first macros

# In[1]:


# this line is useless, will be added by export filter anyway. here for redundant checks
# -*- coding: utf-8 -*-

import numpy as np  # so far is the only dependable module to be included
import os # to check if file exists
import re 
################################################################################


# functions: parse
# 

# ## Parse  of sparam data from  ICCAP native mdm files

# In[303]:


def parse(filename, readAllKeys=False): 
        """
         this routines opens the file given as input parameter then:
          - check that the file actually exists 
          - chack that files cames with proper ICCAP headers
          - reads in the header variables and stores them in a hash
          - reads the column headers and store the column names in a hash
        """
        header=dict()
        DBs=list()
        lines=list()
        with open(filename, encoding='Latin1') as f:
            lines=f.readlines()
        lines=list(map(str.strip, lines))

        if not re.search("^! VERSION = (.*)", lines[0]):
            print(f"FATAL:missing first line check! {lines[0]}")
            return
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
        
        for i, line in enumerate(header):
            if re.search("^.*ICCAP_INPUT",line): inputvars=i
            if re.search("^.*ICCAP_OUTPUT", line): outputvars=i
        inlist=header[inputvars+1:outputvars]
        outlist=header[outputvars+1:]
        # print(inlist)
        inputs=[line.split() for line in inlist]
        outputs=[line.split() for line in outlist]
        # print(outputs)
        varBias=list()
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
        if not fsname=='F':
            print(f"WARNING: frequency sweep variable name is '{fsname}', will rename it to 'F'")
        if len(varBias)>0:
            print(f"File contains bias sweeps on {varBias}")
        tosave=list()
        for i, var in enumerate(outputs):
            if var[1]=="S":
                spname=var[0]
            if var[1]=='V':
                tosave.append(var[0])
            if var[1]=='I':
                tosave.append(var[0])   
        if not spname=='Sp':
            print(f"WARNING: Sparam var name is '{spname}' will rename it to 'Sp'")
        if len(tosave)>0:
            print(f"found extra output variables named {tosave}")
        
        
        dbstart=list()
        dbstop=list()
        for i,line in enumerate(lines):
            if re.search("^BEGIN_DB",line): dbstart.append(i+1)
            if re.search("^END_DB", line): dbstop.append(i)
            
        DBsraw=list()
        for i, start in enumerate(dbstart):
            DBsraw.append(lines[start:dbstop[i]])
        colstring=DBsraw[0][len(tosave)+1]
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
        #print(save)
        
        DBs=list()
        
        # parsebias points
        for DB in DBsraw:
            if not DB[len(tosave)+1]==colstring:
                print("FATAL: sweep header are not matching!")
                return
            DBhead=DB[0:len(tosave)]
            #print(DBhead)
            for line in DBhead:
                b=re.search('^\s*ICCAP_VAR (\S+)\s+(\S+)',line)
                if b:
                    # gropu 0 is the string
                    # group (1) is the first parentheses
                    # group (2) is the second parentheses
                    save[b.group(1)].append(float(b.group(2))+0.0)
        
            DBdata=DB[len(tosave)+2:] # acconts fot an empty lne and a column header
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
            save[var]=np.array(save[var])
        return save            
            
     


# In[304]:


# d=preparse("C:/home/LABLUCCI/230707_STC65_ALTAIR_IEMN/Q1FHF26_6_ret23/S_parameters.mdm")


# In[ ]:





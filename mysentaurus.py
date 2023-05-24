#!/usr/bin/env python
# coding: utf-8
# original code taken from 
# "C:\home\Investigations\220105_PiN_GaN_Tours\01_TCAD\01_parsing data.ipynb"

import re
import shlex

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def parse_plt_in_csv(filename="", trans={"cathode OuterVoltage":"Vc",
                                         "anode OuterVoltage":"Va",
                                         'anode TotalCurrent':'Ia', 
                                         'cathode TotalCurrent':'Ic',
                                         'anode eCurrent':'Iae',
                                         'anode hCurrent':'Iah',
                                         'cathode eCurrent':'Ice',
                                         'cathode hCurrent':'Ich'}, debug=False, skip=2):
    """
    parse a plt output file. original plt are also plain ascii easy to parse. will do.
    by default first line is datafiled and second line is dtatype. 
    but if you have some step '0' simulations (Initial condition for bias),
    You may want to skip more lines ...
    """
    return parse_csv(filename=filename, trans=trans, debug=debug, skip=skip)

def parse_cut_in_csv(filename="", trans={"X":"x", 
       'ElectrostaticPotential':'V', 
       'DonorConcentration':'Nd', 
       'AcceptorConcentration':'Na',
       'ValenceBandEnergy':'Ev',
       'ConductionBandEnergy':'Ec',
       'hDensity':'h',
       'eDensity':'e',
       'BandGap':'Eg',
       'hQuasiFermiPotential':'Efh',
       'eQuasiFermiPotential':'Efe',
       'eMobility':'mue',
       'hMobility':'muh'}, debug=False, skip=2):
    """
    parse a cut saved from a tdr outputfile. 
    by default first line is datafiled and second line is dtatype. 
    but if you have some step '0' simulations (Initial condition for bias),
    You may want to skip more lines ...
    """
    return parse_csv(filename=filename, trans=trans, debug=debug, skip=skip)

def parse_csv(filename="", trans={}, debug=False, skip=2):
    """
    parse a plt output file. original plt are also plain ascii easy to parse. will do.
    by default first line is datafiled and second line is dtatype. 
    but if you have some step '0' simulations (Initial condition for bias),
    You may want to skip more lines ...
    """
    with open(filename, encoding='Latin1') as fin:
        lines=fin.read().splitlines()
    datafields=lines[0].split(",")
    data=dict()
    for field in datafields: data[field]=list()
           
    if debug: 
        print(f"datafields={datafields}")
        if skip!=2:  print(f"New skip!")
            
    rest=lines[skip:]
    for i,line in enumerate(rest):
        varlist=line.split(",")
        if len(varlist)==1: continue # skip empty or white space only lines
        if len(varlist)!=len(datafields):
            print(f"soemthing wrong in line: '{line}'")
            continue
        for idx,field in enumerate(datafields): 
            data[field].append(float(varlist[idx]))
        # print(f"{i}, {line}")
    
    for key, newk in trans.items():
        if key in data:
            if debug: print(f"adding {newk} as link to {key}")
            data[newk]=data[key]
    for key, item in data.items():
        data[key]=np.array(item)
    return data



def parse_plt(filename="", trans={}, debug=False):
    """
    parse a plt output file. ...
    """
    lines=list()
    with open(filename, encoding='Latin1') as fin:
        lines=fin.read().splitlines()
    lines=[line.strip() for line in lines]

    if lines[0] != "DF-ISE text":
        print(f"mismatch in HEADER:{lines[0]}")
        return 
    
    head=list()
    data=list()
    
    # try to have header in head and data in data: very bad parsing
    stp=set()
    for i, line in enumerate(lines):
        if line=="Info {":
            hst=i
        if line=="}":
            stp.add(i)
        if line=="Data {":
            dst=i
    # print(f"{hst} {dst} {stp}")
     
    head=lines[hst+1:sorted(stp)[0]]
    data=lines[dst+1:sorted(stp)[-1]]
    # return head, data
    
    ### parse the header
    stp=set()
    for i, line in enumerate(head):
        if re.match("version\s*=\s*1.0", line):
            # print(line) OK
            continue
        if re.match("type\s*=\s*xyplot", line):
            # print(line) OK
            continue
        if re.match("datasets\s*=\s*\[", line):
            dst=i+1
        if re.match("functions\s*=\s*\[", line):
            fst=i+1
        if re.match(".* ]\s*", line):            
            stp.add(i+1)
    # print(f"{dst} {fst} {stp}")
    datasets=head[dst:sorted(stp)[0]]
    functions=head[fst:sorted(stp)[-1]]
    
    datasets[-1]=datasets[-1][:-2] # drops ' ]' a tend of last string
    functions[-1]=functions[-1][:-2]
    # print(datasets)
    # print(functions)
    dsets=[  shlex.split(line) for line in datasets  ]
    func =[  shlex.split(line) for line in functions ]
    
    rawdata=dict()
    
    dsetsl=sum(dsets, []) # list of lists to list
    for name in dsetsl:
        rawdata[name]=list()
    i=0
    l=len(data)
    while i<l:
        for dline in dsets:
            locdata=data[i].split()
            if len(locdata)!=len(dline):
                print(f"dataset in header an data not consistent: {data[i]}\n{locdata}\n{dline}")
                return dsets,data
            i=i+1
            for j,name in enumerate(dline):
                rawdata[name].append(float(locdata[j])) 
    for name in rawdata.keys():
        rawdata[name]=np.array(rawdata[name])
    rawdata['datasets']=dsets
    rawdata['functions']=func
    return rawdata 

def dump_plt(filename, data):
    """
    data must be a dictionary. It must contain a datasets and functions addressing all contained nparray
    use the datastructure loaded from parse_plt as example
    """
    if not type(data)==dict:
        print(f"routins is expecting data to be a dictionary")
        return
    if not 'datasets' in data:
        print(f"could not find datasets structure")
        return
    if not 'functions' in data:
        print(f"could not find functions structure")
        return
    
    with open(filename, "w", newline='\n') as f:
        print("DF-ISE text\n", file=f)
        
        print("Info {\n  version   = 1.0\n  type      = xyplot\n  datasets  = [", file=f)
        for j, line in enumerate(data['datasets']):
            print(f'    "{line[0]}"', end='', file=f)
            for i, name in enumerate(line):
                if i==0: continue ## contorted to be compliant for spaces
                print(f' "{name}"', end='', file=f)
            if j==len(data['datasets'])-1:
                print(f' ]', file=f)
            else:
                print("", file=f)

        print("  functions = [", file=f)
        for j,line in enumerate(data['functions']):
            print(f'    {line[0]}', end='', file=f)
            for i, name in enumerate(line):
                if i==0: continue ## contorted to be compliant for spaces
                print(f' {name}', end='', file=f)
            if j==len(data['functions'])-1:
                print(f' ]', file=f)
            else:
                print("", file=f)
        print("}\n", file=f)
        
        print("Data {", file=f)
        
        ps=data['datasets'][0][0]
        #print(ps)
        for ln in range(len(data[ps])):
            for j, line in enumerate(data['datasets']):
                print(f'    {data[line[0]][ln]:22.14E}', end='', file=f)
                for i, name in enumerate(line):
                    if i==0: continue
                    print(f' {data[name][ln]:22.14E}', end='', file=f)
                print('', file=f)
    
        print("}",file=f)


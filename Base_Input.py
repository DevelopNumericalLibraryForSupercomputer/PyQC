import sys

def get_default_values():
    # key names starting with 'N****' should have integer number values
    # key names starting with 'V****' should have real number values
    # all others should have characters values
    EnvVal={}
    EnvVal['NIRREP']=1
    EnvVal['NBAS']=0
    EnvVal['NOCC']=0
    EnvVal['NVRT']=0
    EnvVal['REF_ORB']='RHF'
    EnvVal['NMAX_ITER']=50
    EnvVal['NROOT']=1
    EnvVal['VTOL_ENG']=1.0E-7    #Energy Tolerence
    EnvVal['NDIM_SUBSP']=10      #Maximum subspace dimension
    EnvVal['GUESS_TYPE']='HDIAG' # Hdiag/CIS_FILE
    EnvVal['NDIM_GUESS']=1
    EnvVal['HBAR_OUT']='FILE'
    EnvVal['HBAR_SUBMO']='DIRECT'
    EnvVal['HBAR_TYPE']='FILE'
    EnvVal['HBAR_DEBUG']='FALSE'
    EnvVal['DATA_DIR']='NA'
    return EnvVal

def read_input_file(EnvVal,argv):
    # this will be moved to a separte file
    finp=argv[1]
    f=open(finp,'r')
    for line in f:
        EnvVal=read_input_line(EnvVal,line)
    f.close()
    
    #print(EnvVal)
    return EnvVal


def read_input_line(EnvVal,line):
    #print(line) 
    if (line[0]!='#'):
       key=line.split('=')[0].upper()
       val=line.split('=')[1].strip()
       #print('key='+key)
       #print('val='+val)
       if key=='DATA_DIR':
          EnvVal[key]=val
       elif key[0]=='N':
          EnvVal[key]=int(val)
       elif key[0]=='V':
          EnvVal[key]=float(val)
       else:
          EnvVal[key]=val.upper()
    return EnvVal


def print_inpkey_values(EnvVal):
    print(' - Input keys:') 
    print('   '+'-'*45)
    for key, value in EnvVal.items():
        print('   %-20s : %-20s' % (key, str(value)))
    print('   '+'-'*45)

def check_inpkey_values(EnvVal):
    for key, value in EnvVal.items():
        if (key=='NOCC' or key=='NVRT'):
           if (EnvVal[key]==0):
              print(' - Error: '+str(key)+' is not set.')
              exit(1)

def get_inpkey_values(EnvVal):
    if (EnvVal['NBAS']==0): 
       EnvVal['NBAS']=EnvVal['NOCC']+EnvVal['NVRT']
    return EnvVal

def driver(argv):
    print(' * Input module') 
    #print(argv)
    if (len(argv)==2):
       print(' - Input file: '+str(argv[1]))

       EnvVal=get_default_values()
       EnvVal=read_input_file(EnvVal,argv)

       check_inpkey_values(EnvVal)
       EnvVal=get_inpkey_values(EnvVal)
       print_inpkey_values(EnvVal)
    else:
       print(' - Error: One input file should be specified.')
       print('   [usage] python EOMCCSD.py input_file')
       exit(1)

    return EnvVal


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
    EnvVal['NMAXITER']=20
    EnvVal['NROOT']=1
    EnvVal['VTOL_ENG']=1.0E-7    #Energy Tolerence
    EnvVal['NDIM_SUBSP']=5       #Maximum subspace dimension
    EnvVal['GUESS_TYPE']='HDIAG' # Hdiag/CIS_FILE
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
       inpkey=line.split('=')[0].upper()
       inpval=line.split('=')[1].strip()
       #print('key='+inpkey)
       #print('val='+inpval)
       if inpkey=='DATA_DIR':
          EnvVal[inpkey]=inpval
       elif inpkey[0]=='N':
          EnvVal[inpkey]=int(inpval)
       elif inpkey[0]=='V':
          EnvVal[inpkey]=float(inpval)
       else:
          EnvVal[inpkey]=inpval.upper()
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


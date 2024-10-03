
def get_default_values():
    # this will be moved to a separte file
    EnvVal={}
    EnvVal['Irrep']=1
    EnvVal['Nbas']=11
    EnvVal['Nocc']=2
    EnvVal['Nvrt']=9
    EnvVal['RefOrb']='RHF'
    EnvVal['MaxIter']=20
    EnvVal['Nroot']=1
    EnvVal['EngTol']=1.0E-5    #Energy Tolerence
    EnvVal['DavSubSpDim']=5    #Maximum subspace dimension
#   EnvVal['GuessType']='CIS_FILE'
    EnvVal['GuessType']='Hdiag'
    EnvVal['HbarType']='FILE'
    return EnvVal

def get_input_file(EnvVal,finp):
    f=open(finp,'r')
    for line in f:
        print(line) 
    f.close()
    return EnvVal

def driver():
    EnvVal=get_default_values()
    return EnvVal


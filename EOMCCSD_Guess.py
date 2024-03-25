import sys
import numpy as np
import EOMCCSD_Util as util

def make_Hdiag(EnvVal,F):
    # diagonal approximation
    EngOcc=F['oo'].diagonal()
    EngVrt=F['vv'].diagonal()
    D1a = EngOcc.reshape(-1,1)-EngVrt
    D2aa= EngOcc.reshape(-1,1,1,1)+EngOcc.reshape(-1,1,1) \
         -EngVrt.reshape(-1,1)-EngVrt

    D1a =D1a.flatten()
    D2aa=D2aa.flatten()
    D=np.concatenate((D1a,D2aa),axis=0)
    return D

def sort_and_get_indices(Eval,Evec):
    indexed_Eval = list(enumerate(Eval))
    sorted_Eval = sorted(indexed_Eval, key=lambda x: x[1])
    sorted_indices = [index for index, _ in sorted_Eval]
    #print(sorted_Eval)
    #print(sorted_indices)

    dim=len(Eval)
    Eval_new=np.zeros(dim)
    Evec_new=np.zeros((dim,dim))
    for i in range(dim):
        isort=sorted_indices[i]
        Eval_new[i]=Eval[isort]           #reorder eig. values
        for j in range(dim):
            jsort=sorted_indices[j]
            Evec_new[i][j]=Evec[i][jsort] #reorder eig. vectors

    lprint=False
    if (lprint):
       print('Val, original')
       print(Eval)
       print('Val, sorted')
       print(Eval_new)
       print('Vec, original')
       print(Evec)
       print('Vec, sorted')
       print(Evec_new)

    Evec_new=Evec_new.transpose()  #transpose for convenience

    return Eval_new,Evec_new

def get_CIS(EnvVal,lprint):
    Nocc=EnvVal['Nocc']
    Nvrt=EnvVal['Nvrt']
    RefOrb=EnvVal['RefOrb']
    GuessType=EnvVal['GuessType']
    Nov=Nocc*Nvrt

    if (GuessType=='CIS_ACES2'):
       print(' - Reading CIS matrix from ACES2')
       #if aces2py not in sys.modules:
       import EOMCCSD_AcesPy as ap 
       EigVal,EigVec = ap.get_CIS_matrix(Nocc,Nvrt)
       EigVec=EigVec.reshape(Nov,Nov)
       Val,Vec = sort_and_get_indices(EigVal,EigVec)
    elif (GuessType=='CIS_FILE'):
       print(' - Reading CIS matrix from file')
       finp='CIS-Matrix' 
       CISMat=util.read_data(finp,True)
       #print('CIS Matrix')
       #print(CISMat)
       EigVal,EigVec = np.linalg.eig(CISMat)
       EigVal=np.real(EigVal)
       EigVec=np.real(EigVec)
       if RefOrb=='RHF': 
          EigVec=EigVec.reshape(Nov,Nov)
       elif RefOrb=='UHF': 
          EigVec=EigVec.reshape(Nov*2,Nov*2)

    if (lprint):
       print('\n * CIS Eigenvalue')
       print(EigVal)
       print('\n * CIS Eigemvector, len='+str(len(Vec)))
       print(EigVec)

    return EigVal,EigVec 

def get_CISvec(EnvVal): 
    Nocc=EnvVal['Nocc']
    Nvrt=EnvVal['Nvrt']
    RefOrb=EnvVal['RefOrb']
    GuessType=EnvVal['GuessType']
    Nov=Nocc*Nvrt
    Nroot=1
    
    finp='ACES2-CISvec' 
    CISvec=util.read_data(finp,True)

    print('CISvec')
    print(CISvec)
 
    Gvec=np.zeros((Nroot,Nov))
    Gvec[0,:]=CISvec

    return Gvec


def driver(EnvVal,F,W):
    Nocc=EnvVal['Nocc']
    Nvrt=EnvVal['Nvrt']
    Nroot=EnvVal['Nroot']
    GuessType=EnvVal['GuessType']
    Nov=Nocc*Nvrt
    Rdim=Nov+Nov*Nov
    NGuessSp=2
    lprint=False   

    #initial R
    if (GuessType=='Hdiag'):
       print('\n * Guess : diagonal approxiation')
       Hdiag=make_Hdiag(EnvVal,F)
       idx=Hdiag[:Nov].argsort()[::-1][:Nroot*NGuessSp]
       R = np.eye(Rdim)[:,idx]

    elif (GuessType[0:3]=='CIS'):
       print('\n * Guess : from CIS vectors')
       Val,Vec=get_CIS(EnvVal,lprint)
       idx=Val.argsort()[:Nroot*NGuessSp]
       R = np.zeros((Rdim,Nroot*NGuessSp))
       for i in range(Nroot*NGuessSp):
           R[:Nov,i] = Vec[:Nov,idx[i]]

    else:
       print('\n * Guess : Error, Unknown Guess type = '+GuessType)

    return R
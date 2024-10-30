import sys
import numpy as np
import Base_Util as util

def make_Hdiag(EnvVal,F):
    #print('Hdiag : diagonal approximation')
    EngOcc=F['oo'].diagonal()
    EngVrt=F['vv'].diagonal()

    D1a = EngOcc.reshape(-1,1)-EngVrt.reshape(1,-1)  #e(i)-e(a)
    D2oo= EngOcc.reshape(-1,1,1,1)+EngOcc.reshape(-1,1,1)
    D2vv= EngVrt.reshape(-1,1)+EngVrt.reshape(1,-1)
    D2ab= D2oo-D2vv  #e(i)+e(j)-e(a)-e(b)

    D1a =D1a.flatten()
    D2ab=D2ab.flatten()
    D=np.concatenate((D1a,D2ab),axis=0)
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
    Nocc=EnvVal['NOCC']
    Nvrt=EnvVal['NVRT']
    RefOrb=EnvVal['REF_ORB']
    GuessType=EnvVal['GUESS_TYPE']
    Nov=Nocc*Nvrt

    if (GuessType=='CIS_ACES2'):
       print(' - Reading CIS matrix from ACES2')
       #if aces2py not in sys.modules:
       import Base_AcesPy as ap 
       EigVal,EigVec = ap.get_CIS_matrix(Nocc,Nvrt)
       EigVec=EigVec.reshape(Nov,Nov)
       Val,Vec = sort_and_get_indices(EigVal,EigVec)
    elif (GuessType=='CIS_FILE'):
       print(' - Reading CIS matrix from file')
       finp='CIS-Matrix' 
       CISMat=util.read_data(finp,EnvVal)
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
    Nocc=EnvVal['NOCC']
    Nvrt=EnvVal['NVRT']
    RefOrb=EnvVal['REF_ORB']
    GuessType=EnvVal['GUESS_TYPE']
    Nov=Nocc*Nvrt
    Nroot=1
    
    finp='ACES2-CISvec' 
    CISvec=util.read_data(finp,EnvVal)

    print('CISvec')
    print(CISvec)
 
    Gvec=np.zeros((Nroot,Nov))
    Gvec[0,:]=CISvec

    return Gvec


def driver(EnvVal,F,W):
    Nocc=EnvVal['NOCC']
    Nvrt=EnvVal['NVRT']
    Nroot=EnvVal['NROOT']
    GuessType=EnvVal['GUESS_TYPE']
    NdimGuess=EnvVal['NDIM_GUESS']
    Nov=Nocc*Nvrt
    Rrow=Nov+Nov*Nov
    Rcol=Nroot*NdimGuess
    lprint=False   

    #initial R
    if (GuessType=='HDIAG'):
       print('\n * Guess : diagonal approximation')
       Hdiag=make_Hdiag(EnvVal,F)
       D1a=Hdiag[:Nov]
       idx=D1a.argsort()[::-1]  #D1a(ordered)
       idx=idx[:Rcol] #Choose the M lowest elements (M=Nroor*NdimGuess)

       R = np.zeros([Rrow,Rcol])
       for i in range(Rcol):
           R[idx[i],i]=1.0

    elif (GuessType[0:3]=='CIS'):
       print('\n * Guess : from CIS vectors')
       Val,Vec=get_CIS(EnvVal,lprint)
       idx=Val.argsort()[:Rcol]

       R = np.zeros((Rrow,Rcol))
       for i in range(Rcol):
           R[:Nov,i] = Vec[:Nov,idx[i]]

    else:
       print('\n * Guess : Error, Unknown Guess type = '+GuessType)
       exit(1)
       
    print(' - Size of the guess vector (Rrow,Rcol) = '+str(R.shape))
    print('   Rrow: Rov+Roovv       = '+str(Nov)+'+'+str(Nov*Nov))
    print('   Rcol: Nroot*NdimGuess = '+str(Nroot*NdimGuess)) 
    return R

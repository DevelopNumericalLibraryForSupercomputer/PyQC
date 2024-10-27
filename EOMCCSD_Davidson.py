import os
import numpy as np
import Base_Util as util
import EOMCCSD_Sigma as sig

def Diag_Davidson(EnvVal,F,W,T,R):
    Nocc=EnvVal['NOCC']
    Nvrt=EnvVal['NVRT']
    RefOrb=EnvVal['REF_ORB']
    MaxIter=EnvVal['NMAXITER']
    Nroot=EnvVal['NROOT']  #for now, Nroot=1
    EngTol=EnvVal['VTOL_ENG']
    DavSubSpDim=EnvVal['NDIM_SUBSP'] #Maximum subspace dimension
    Nov=Nocc*Nvrt
    DimR=Nov+Nov*Nov
    NGuessSp=2

    Hdiag=make_Hdiag(EnvVal,F,W)

    Eng=[0.0]*Nroot*NGuessSp

    print('\n * EOM-CCSD iteration, Max='+str(MaxIter))
    for Iter in range(MaxIter): 

        # make each R vector orthogonal with QR decomp. 
        Qval,Rval=np.linalg.qr(R)
        R=Qval

        # save energies
        Eng_old=Eng[:Nroot] 

        # make the subspace and diagonalize
        DimG=R.shape[1]
        G=np.zeros((DimG,DimG))
        Z=np.zeros((DimR,DimG))
        for i in range(DimG):
            Z[:,i]=sig.make_Sigma(EnvVal,F,W,T,R[:,i])
            for j in range(i+1):
                G[j,i]=np.dot(R[:,j],Z[:,i])
                G[i,j]=G[j,i]
        Eng,Evec=np.linalg.eig(G)

        # select root (by energy) 
        idx=Eng.argsort()[:Nroot]
        Eng=np.real(Eng[idx])
        Evec=np.real(Evec[:,idx])

        # get new vector
        for i in range(Nroot):
            ResVec = np.dot(Z-Eng[i]*R, Evec[:,i])
            CorrVec=ResVec/(Eng[i]-Hdiag)
            dE=abs(Eng[i]-Eng_old[i])
            print(' - #%d :  E[%d] = %.10f    dE = %.10f ' % (Iter,i+1,Eng[i],dE))

        # check conv. / if not, update R
        if (dE < EngTol) and (Iter>1): 
           print('\n * EOM-CCSD Converged')
           for i in range(Nroot):
               print(' - Energy for the state %d = %.10f ' % (i+1,Eng[i]))
           return
        else:   
           if DimG >= DavSubSpDim:
              R=np.dot(R,Evec)    # start the new subspace with the last estimate
              Eng_old=Eng
           else: 
              R=np.c_[R,CorrVec]  # update R 

    return


def make_Hdiag(EnvVal,F,W):
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


